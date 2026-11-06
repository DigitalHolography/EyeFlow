"""Resolved execution-plan exports, independent of scientific execution and UI."""

import json
from pathlib import Path

from pipeline_engine.planning.dag import PipelineDAG


def resolve_execution_graph(spec, *, observed_output=None):
    """Describe selected graphs; variant templates remain explicitly dynamic."""
    dag = PipelineDAG(spec.plan.descriptors)
    varying = set()
    producers = [item.name for item in spec.plan.descriptors if item.produces_execution_variants]
    for name in producers:
        varying.update(
            dag.dependents_of(name, transitive=True, pipeline_options=spec.pipeline_options)
        )
    pipelines, edges, plans = [], [], {}
    for descriptor in spec.plan.descriptors:
        name = descriptor.name
        graph = descriptor.build_graph(
            spec.pipeline_options.get(name, ()), selection=spec.selection
        )
        plan = (
            graph.plan(external_products=(item.alias for item in descriptor.product_requires))
            if graph
            else None
        )
        plans[name] = plan
        stages = []
        if plan is not None:
            for stage in graph.stages:
                stages.append(
                    dict(
                        id=stage.id,
                        artifact=stage.artifact,
                        wait_for_artifacts=stage.wait_for_artifacts,
                        products=[
                            dict(id=p.id, optional=p.optional, exposed=p.exposed)
                            for p in stage.products
                        ],
                    )
                )
                for kind, requirements in (
                    ("required_product", stage.requires),
                    ("optional_product", stage.optional_requires),
                ):
                    for product in requirements:
                        parent = plan.producer_by_product[product]
                        if parent is not None:
                            edges.append(
                                dict(
                                    source=f"{name}/{parent}",
                                    target=f"{name}/{stage.id}",
                                    kind=kind,
                                    product=product,
                                )
                            )
                edges.extend(
                    dict(
                        source=f"{name}/{parent}",
                        target=f"{name}/{stage.id}",
                        kind="successful_stage",
                    )
                    for parent in stage.after
                )
        pipelines.append(
            dict(
                id=name,
                execution_mode=descriptor.execution_mode,
                scope=(
                    "per_successful_variant"
                    if name in varying and not descriptor.run_once
                    else "shared"
                ),
                wait_for_artifacts=descriptor.wait_for_artifacts,
                options=list(spec.pipeline_options.get(name, ())),
                stages=stages,
            )
        )
        for parent in dag.dependencies_of(
            name, pipeline_options=spec.pipeline_options, include_products=False
        ):
            edges.append(
                dict(source=f"{parent}/completion", target=f"{name}/entry", kind="whole_pipeline")
            )
    for descriptor in spec.plan.descriptors:
        consumer = plans[descriptor.name]
        for requirement in descriptor.product_requires:
            producer = plans[requirement.producer_pipeline]
            source = f"{requirement.producer_pipeline}/{producer.producer_by_product[requirement.product]}"
            for stage in consumer.stage_by_id.values():
                if requirement.alias in (*stage.requires, *stage.optional_requires):
                    edges.append(
                        dict(
                            source=source,
                            target=f"{descriptor.name}/{stage.id}",
                            kind="cross_pipeline_product",
                            product=requirement.product,
                            optional=requirement.alias in stage.optional_requires,
                            variant=requirement.variant,
                        )
                    )
        if descriptor.wait_for_artifacts:
            edges.extend(
                dict(
                    source=f'{item["id"]}/artifacts',
                    target=f"{descriptor.name}/entry",
                    kind="artifact_barrier",
                )
                for item in pipelines
                if item["id"] != descriptor.name
                and any(stage["artifact"] for stage in item["stages"])
            )
    result = dict(
        format="eyeflow.execution_graph.v1",
        targets=list(spec.plan.targets),
        pipelines=pipelines,
        variant_discovery=dict(
            producers=producers, requested=spec.velocity_estimation_method, successful=None
        ),
        edges=edges,
        observations=None,
    )
    if observed_output is not None:
        import h5py
        from app_settings import selected_velocity_methods

        with h5py.File(observed_output, "r") as h5:
            observed_options = json.loads(h5.attrs.get("pipeline_options", "{}"))
            expected_options = {
                name: list(options) for name, options in spec.pipeline_options.items()
            }
            if (
                observed_options != expected_options
                or list(h5.attrs.get("pipeline_targets", ())) != list(spec.plan.targets)
                or list(h5.attrs.get("velocity_estimation_methods", ()))
                != list(selected_velocity_methods(spec.velocity_estimation_method))
            ):
                raise ValueError("Observed run has different selected targets or options.")
            observed = {
                key: json.loads(h5.attrs.get(key, "{}"))
                for key in (
                    "stage_status",
                    "stage_failures",
                    "stage_skips",
                    "stage_products",
                    "stage_timings",
                    "artifact_status",
                    "artifact_failures",
                    "artifact_timings",
                )
            }
            for scope, statuses in observed["stage_status"].items():
                name = scope.split(" [", 1)[0]
                if (
                    name not in plans
                    or plans[name] is None
                    or set(statuses) != set(plans[name].stage_by_id)
                ):
                    raise ValueError(f"Observed graph differs from resolved settings: {scope}")
            result["variant_discovery"]["successful"] = list(
                h5.attrs.get("execution_variants_completed", ())
            )
            result["observations"] = observed
    return result


def execution_graph_dot(graph):
    """Render dependency kinds; boundary diamonds are ports, never runtime stages."""
    quote = lambda value: json.dumps(str(value))
    lines = [
        "digraph EyeFlow {",
        "  rankdir=LR;",
        "  compound=true;",
        '  node [shape=box, fontname="Segoe UI"];',
    ]
    observations = graph["observations"] or {}
    for pipeline in graph["pipelines"]:
        name = pipeline["id"]
        label = name + (" (per successful variant)" if pipeline["scope"] != "shared" else "")
        lines.extend([f'  subgraph {quote("cluster_"+name)} {{', f"    label={quote(label)};"])
        for port in ("entry", "completion", "artifacts"):
            lines.append(
                f'    {quote(name+"/"+port)} [shape=diamond, fontsize=8, label={quote(port)}];'
            )
        for stage in pipeline["stages"]:
            label = stage["id"] + (" [artifact]" if stage["artifact"] else "")
            statuses = []
            for scope, records in observations.get("stage_status", {}).items():
                if scope.split(" [", 1)[0] == name and stage["id"] in records:
                    timing = (
                        observations.get("stage_timings", {}).get(scope, {}).get(stage["id"], {})
                    )
                    label += f'\n{scope}: {records[stage["id"]]}'
                    if "wall_seconds" in timing:
                        label += f' {timing["wall_seconds"]:.3f}s'
                    statuses.append(records[stage["id"]])
            color = (
                "red"
                if "failed" in statuses
                else ("darkgreen" if statuses and set(statuses) == {"succeeded"} else "black")
            )
            lines.append(
                f'    {quote(name+"/"+stage["id"])} [label={quote(label)}, color={quote(color)}, penwidth={2 if any(p["exposed"] for p in stage["products"]) else 1}];'
            )
        lines.append("  }")
        nodes = {f'{name}/{stage["id"]}' for stage in pipeline["stages"]}
        incoming = {edge["target"] for edge in graph["edges"] if edge["source"] in nodes}
        outgoing = {edge["source"] for edge in graph["edges"] if edge["target"] in nodes}
        for node in sorted(nodes - incoming):
            lines.append(
                f'  {quote(name+"/entry")} -> {quote(node)} [color=gray, label="boundary"];'
            )
        for node in sorted(nodes - outgoing):
            lines.append(
                f'  {quote(node)} -> {quote(name+"/completion")} [color=gray, label="boundary"];'
            )
        for stage in pipeline["stages"]:
            if stage["artifact"]:
                lines.append(
                    f'  {quote(name+"/"+stage["id"])} -> {quote(name+"/artifacts")} [color=purple, label="export completion"];'
                )
    styles = dict(
        required_product=("black", "solid"),
        optional_product=("gray", "dashed"),
        successful_stage=("black", "dotted"),
        whole_pipeline=("black", "bold"),
        cross_pipeline_product=("blue", "solid"),
        artifact_barrier=("purple", "dashed"),
    )
    for edge in graph["edges"]:
        color, style = styles[edge["kind"]]
        label = edge.get("product", edge["kind"])
        lines.append(
            f'  {quote(edge["source"])} -> {quote(edge["target"])} [color={color}, style={style}, label={quote(label)}];'
        )
    return "\n".join([*lines, "}"]) + "\n"


def export_execution_graph(spec, destination, *, observed_output=None):
    """Write sibling JSON/DOT files; never calculate data or replace run outputs."""
    graph = resolve_execution_graph(spec, observed_output=observed_output)
    base = Path(destination)
    base.parent.mkdir(parents=True, exist_ok=True)
    json_path, dot_path = base.with_suffix(".json"), base.with_suffix(".dot")
    json_path.write_text(json.dumps(graph, indent=2), encoding="utf-8")
    dot_path.write_text(execution_graph_dot(graph), encoding="utf-8")
    return json_path, dot_path
