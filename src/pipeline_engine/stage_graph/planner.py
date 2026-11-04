"""Deterministic local DAG planning, independent of global pipeline planning."""

from dataclasses import dataclass
from heapq import heappop, heappush
from types import MappingProxyType


@dataclass(frozen=True, slots=True)
class StagePlan:
    order: tuple
    stage_by_id: object
    producer_by_product: object
    required: object
    optional: object
    successful_stages: object
    _parents: object
    _children: object

    def dependencies(self, stage_id):
        return self._parents[stage_id]

    def dependents(self, stage_id):
        return self._children[stage_id]

    def ancestors(self, stage_id):
        found = set()
        pending = list(self.dependencies(stage_id))
        while pending:
            parent = pending.pop()
            if parent not in found:
                found.add(parent)
                pending.extend(self.dependencies(parent))
        return tuple(name for name in self.order if name in found)


class StageGraphPlanner:
    def plan(self, graph, *, external_products=()):
        stages, producers, output_paths = {}, {name: None for name in external_products}, []
        for stage in graph.stages:
            if not isinstance(stage.id, str) or not stage.id.strip() or stage.id in stages:
                raise ValueError(f"Duplicate or empty stage ID: {stage.id!r}")
            if not callable(stage.run):
                raise ValueError(f"Non-callable stage: {stage.id}")
            stages[stage.id] = stage
            for product in stage.products:
                if not isinstance(product.id, str) or not product.id.strip() or product.id in producers:
                    raise ValueError(f"Duplicate or empty product ID: {product.id!r}")
                if product.path is not None and product.pack is not None:
                    raise ValueError(f"Product {product.id} declares both path and pack.")
                if any(callback is not None and not callable(callback)
                       for callback in (product.pack, product.validate)):
                    raise ValueError(f"Invalid product callback: {product.id}")
                if product.path is not None:
                    path = product.path.strip("/")
                    if not path:
                        raise ValueError("Product output path cannot be empty.")
                    if any(path == other or path.startswith(other + "/") or other.startswith(path + "/")
                           for other in output_paths):
                        raise ValueError("Product output specifications must be non-overlapping.")
                    output_paths.append(path)
                producers[product.id] = stage.id
        parents, children = {key: set() for key in stages}, {key: set() for key in stages}
        for stage in stages.values():
            if set(stage.requires) & set(stage.optional_requires):
                raise ValueError(f"Contradictory required/optional inputs: {stage.id}")
            for group in (stage.requires, stage.optional_requires, stage.after):
                if len(group) != len(set(group)):
                    raise ValueError(f"Duplicate dependency declarations: {stage.id}")
            for name in (*stage.requires, *stage.optional_requires):
                if name not in producers:
                    raise ValueError(f"Stage {stage.id} has undeclared product {name!r}")
                if producers[name] is not None:
                    parents[stage.id].add(producers[name])
            for name in stage.after:
                if name not in stages:
                    raise ValueError(f"Unknown stage-completion dependency: {name}")
                parents[stage.id].add(name)
            if stage.id in parents[stage.id]:
                raise ValueError(f"Stage {stage.id} depends on itself.")
            for parent in parents[stage.id]:
                children[parent].add(stage.id)
        indices = {name: index for index, name in enumerate(stages)}
        degree = {name: len(items) for name, items in parents.items()}
        ready, order = [], []
        for name, count in degree.items():
            if not count:
                heappush(ready, (indices[name], name))
        while ready:
            _, name = heappop(ready)
            order.append(name)
            for child in children[name]:
                degree[child] -= 1
                if not degree[child]:
                    heappush(ready, (indices[child], child))
        if len(order) != len(stages):
            raise ValueError("Cycle detected in stage graph.")
        ordered = lambda items: tuple(name for name in order if name in items)
        immutable = lambda values: MappingProxyType(values)
        return StagePlan(tuple(order), immutable(stages), immutable(producers),
                         immutable({k: v.requires for k, v in stages.items()}),
                         immutable({k: v.optional_requires for k, v in stages.items()}),
                         immutable({k: v.after for k, v in stages.items()}),
                         immutable({k: ordered(v) for k, v in parents.items()}),
                         immutable({k: ordered(v) for k, v in children.items()}))
