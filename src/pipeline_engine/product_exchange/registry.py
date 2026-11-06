"""Coordinator-owned product availability, independent of physical data reads."""

from dataclasses import dataclass
from enum import Enum


class ProductStatus(str, Enum):
    PENDING = "pending"
    AVAILABLE = "available"
    FAILED = "failed"


@dataclass
class ProductEntry:
    status: ProductStatus = ProductStatus.PENDING
    handle: object = None
    commit_id: str | None = None
    reason: str | None = None


class ProductExchangeRegistry:
    def __init__(self):
        self._entries = {}
        self._dependents = {}
        self.revision = 0
        self._invocations = {}

    def register_invocation(self, invocation_id):
        self._invocations.setdefault(invocation_id.pipeline, set()).add(invocation_id)

    def invocations_for(self, pipeline):
        return tuple(sorted(self._invocations.get(pipeline, ()), key=lambda item: item.variant or ""))

    def remove_placeholder(self, invocation_id):
        """Replace a pending shared placeholder once workflow identities are known."""
        owned = [key for key in self._entries if key.producer == invocation_id]
        if any(self._entries[key].status != ProductStatus.PENDING for key in owned):
            raise ValueError("Cannot replace an invocation that already resolved products.")
        for key in owned:
            del self._entries[key]
            self._dependents.pop(key, None)
        self._invocations.get(invocation_id.pipeline, set()).discard(invocation_id)
        for dependents in self._dependents.values():
            dependents.discard(invocation_id)

    def declare(self, product_id):
        self.register_invocation(product_id.producer)
        self._entries.setdefault(product_id, ProductEntry())

    def add_dependent(self, product_id, invocation_id):
        if product_id not in self._entries:
            raise ValueError(f"Undeclared producer product: {product_id}")
        self._dependents.setdefault(product_id, set()).add(invocation_id)

    def validate_accept(self, product_id, handle, commit_id):
        if product_id not in self._entries:
            raise ValueError(f"Undeclared exposed product: {product_id}")
        entry = self._entries[product_id]
        if entry.status == ProductStatus.FAILED:
            raise ValueError("Cannot accept a finalized failed product.")
        if entry.status == ProductStatus.AVAILABLE and (entry.commit_id != commit_id or entry.handle != handle):
            raise ValueError("Conflicting duplicate product delivery.")

    def accept(self, product_id, handle, commit_id):
        self.validate_accept(product_id, handle, commit_id)
        self._entries[product_id] = ProductEntry(ProductStatus.AVAILABLE, handle, commit_id)
        self.revision += 1

    def fail(self, product_id, reason):
        entry = self._entries[product_id]
        if entry.status == ProductStatus.PENDING:
            self._entries[product_id] = ProductEntry(ProductStatus.FAILED, reason=str(reason))
            self.revision += 1

    def finalize_producer(self, invocation_id, reason="Producer finished without supplying product."):
        for product_id in self._entries:
            if product_id.producer == invocation_id:
                self.fail(product_id, reason)

    def status(self, product_id):
        return self._entries[product_id].status

    def get(self, product_id):
        entry = self._entries[product_id]
        return entry.handle if entry.status == ProductStatus.AVAILABLE else None

    def dependents(self, product_id):
        return frozenset(self._dependents.get(product_id, ()))
