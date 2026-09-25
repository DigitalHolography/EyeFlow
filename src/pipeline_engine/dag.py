from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from .base import PipelineDescriptor


@dataclass(frozen=True)
class PipelineExecutionPlan:
    """Resolved pipeline targets in the order they must run."""

    targets: tuple[str, ...]
    descriptors: tuple[PipelineDescriptor, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(pipeline.name for pipeline in self.descriptors)


class PipelineDAG:
    """
    Resolve pipeline dependencies from declared data keys.

    Each pipeline implicitly produces its own pipeline name. Developers can add
    business-level keys with `dag_produces` and consume them with `dag_requires`.
    Required keys with no producer are treated as external inputs.
    """

    def __init__(self, pipelines: Iterable[PipelineDescriptor]) -> None:
        self._pipelines = tuple(pipelines)
        self._pipelines_by_name = self._build_pipeline_index()
        self._original_index = {
            pipeline.name: index for index, pipeline in enumerate(self._pipelines)
        }
        self._key_producers = self._build_key_producers()
        self.graph, self._dependencies = self._build_graph()
        self.execution_order = tuple(self._topological_sort())

    @property
    def ordered_descriptors(self) -> tuple[PipelineDescriptor, ...]:
        return tuple(
            self._pipelines_by_name[name] for name in self.execution_order
        )

    def resolve_targets(
        self,
        targets: Sequence[str],
        *,
        pipeline_options: Mapping[str, Sequence[str]] | None = None,
    ) -> PipelineExecutionPlan:
        target_names = tuple(dict.fromkeys(name for name in targets if name))
        unknown = [
            name for name in target_names if name not in self._pipelines_by_name
        ]
        if unknown:
            available = ", ".join(self.execution_order)
            raise ValueError(
                f"Unknown pipeline target(s): {', '.join(unknown)}. "
                f"Available: {available}"
            )
        if not target_names:
            return PipelineExecutionPlan(targets=(), descriptors=())

        required: set[str] = set()

        def collect(pipeline_name: str) -> None:
            if pipeline_name in required:
                return
            required.add(pipeline_name)
            for dependency in self._active_dependencies(
                pipeline_name,
                pipeline_options,
            ):
                collect(dependency)

        for target in target_names:
            collect(target)

        descriptors = tuple(
            self._pipelines_by_name[name]
            for name in self.execution_order
            if name in required
        )
        return PipelineExecutionPlan(targets=target_names, descriptors=descriptors)

    def dependencies_of(
        self,
        pipeline_name: str,
        *,
        transitive: bool = False,
        pipeline_options: Mapping[str, Sequence[str]] | None = None,
    ) -> tuple[str, ...]:
        """Return direct or transitive upstream pipeline dependencies."""
        dependencies = self._active_dependency_graph(pipeline_options)
        return self._related_pipelines(
            pipeline_name,
            dependencies,
            transitive=transitive,
        )

    def dependents_of(
        self,
        pipeline_name: str,
        *,
        transitive: bool = False,
        pipeline_options: Mapping[str, Sequence[str]] | None = None,
    ) -> tuple[str, ...]:
        """Return direct or transitive downstream pipeline dependents."""
        dependencies = self._active_dependency_graph(pipeline_options)
        graph = {name: set() for name in self._pipelines_by_name}
        for child, parents in dependencies.items():
            for parent in parents:
                graph[parent].add(child)
        return self._related_pipelines(
            pipeline_name,
            graph,
            transitive=transitive,
        )

    def _active_dependency_graph(
        self,
        pipeline_options: Mapping[str, Sequence[str]] | None,
    ) -> dict[str, set[str]]:
        return {
            name: self._active_dependencies(name, pipeline_options)
            for name in self._pipelines_by_name
        }

    def _active_dependencies(
        self,
        pipeline_name: str,
        pipeline_options: Mapping[str, Sequence[str]] | None,
    ) -> set[str]:
        pipeline = self._pipelines_by_name[pipeline_name]
        required_keys = list(pipeline.dag_requires)
        selected_options = self._selected_option_names(
            pipeline,
            pipeline_options,
        )
        for option in pipeline.options:
            if option.name in selected_options:
                required_keys.extend(option.dag_requires)
        return {
            producer
            for key in required_keys
            if (producer := self._key_producers.get(key)) is not None
            and producer != pipeline_name
        }

    @staticmethod
    def _selected_option_names(
        pipeline: PipelineDescriptor,
        pipeline_options: Mapping[str, Sequence[str]] | None,
    ) -> set[str]:
        if pipeline_options is None or pipeline.name not in pipeline_options:
            return {
                option.name
                for option in pipeline.options
                if option.default_enabled
            }
        return {str(name) for name in pipeline_options[pipeline.name]}

    def _related_pipelines(
        self,
        pipeline_name: str,
        relationships: dict[str, set[str]],
        *,
        transitive: bool,
    ) -> tuple[str, ...]:
        if pipeline_name not in self._pipelines_by_name:
            raise ValueError(f"Unknown pipeline: '{pipeline_name}'")
        related = set(relationships[pipeline_name])
        if transitive:
            pending = list(related)
            while pending:
                current = pending.pop()
                for item in relationships[current]:
                    if item not in related:
                        related.add(item)
                        pending.append(item)
        return tuple(name for name in self.execution_order if name in related)

    def _build_pipeline_index(self) -> dict[str, PipelineDescriptor]:
        pipelines_by_name: dict[str, PipelineDescriptor] = {}
        for pipeline in self._pipelines:
            if not pipeline.name:
                raise ValueError("Pipeline DAG contains a pipeline without a name.")
            if pipeline.name in pipelines_by_name:
                raise ValueError(f"Duplicate pipeline name: '{pipeline.name}'")
            pipelines_by_name[pipeline.name] = pipeline
        return pipelines_by_name

    def _build_key_producers(self) -> dict[str, str]:
        key_producers: dict[str, str] = {}
        for pipeline in self._pipelines:
            for key in self._produced_keys(pipeline):
                producer = key_producers.get(key)
                if producer is not None and producer != pipeline.name:
                    raise ValueError(
                        f"Multiple pipelines produce DAG key '{key}': "
                        f"{producer}, {pipeline.name}"
                    )
                key_producers[key] = pipeline.name
        return key_producers

    def _build_graph(self) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
        graph: dict[str, set[str]] = defaultdict(set)
        dependencies: dict[str, set[str]] = defaultdict(set)
        for pipeline_name in self._pipelines_by_name:
            graph[pipeline_name]
            dependencies[pipeline_name]

        for pipeline in self._pipelines:
            required_keys = [
                *pipeline.dag_requires,
                *(
                    key
                    for option in pipeline.options
                    for key in option.dag_requires
                ),
            ]
            for required_key in required_keys:
                producer = self._key_producers.get(required_key)
                if producer is None or producer == pipeline.name:
                    continue
                graph[producer].add(pipeline.name)
                dependencies[pipeline.name].add(producer)

        return dict(graph), dict(dependencies)

    def _topological_sort(self) -> list[str]:
        in_degree = {name: 0 for name in self._pipelines_by_name}
        for children in self.graph.values():
            for child in children:
                in_degree[child] += 1

        queue = deque(
            sorted(
                (name for name, degree in in_degree.items() if degree == 0),
                key=self._original_index.__getitem__,
            )
        )
        order: list[str] = []

        while queue:
            name = queue.popleft()
            order.append(name)
            for child in sorted(
                self.graph[name], key=self._original_index.__getitem__
            ):
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    queue.append(child)

        if len(order) != len(self._pipelines_by_name):
            raise RuntimeError("Cycle detected in pipeline DAG.")
        return order

    @staticmethod
    def _produced_keys(pipeline: PipelineDescriptor) -> tuple[str, ...]:
        return tuple(dict.fromkeys((pipeline.name, *pipeline.dag_produces)))
