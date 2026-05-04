from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Iterable, Sequence
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

    def resolve_targets(self, targets: Sequence[str]) -> PipelineExecutionPlan:
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
            for dependency in self._dependencies[pipeline_name]:
                collect(dependency)

        for target in target_names:
            collect(target)

        descriptors = tuple(
            self._pipelines_by_name[name]
            for name in self.execution_order
            if name in required
        )
        return PipelineExecutionPlan(targets=target_names, descriptors=descriptors)

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
            for required_key in pipeline.dag_requires:
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
