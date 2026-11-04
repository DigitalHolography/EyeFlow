"""Invocation-local logical availability and declared consumer access only."""


class ProductStore:
    def __init__(self, plan):
        self.plan = plan
        self._handles = {}
        self.failures = {}

    def has(self, product_id):
        return product_id in self._handles

    def require(self, product_id):
        if not self.has(product_id):
            raise KeyError(f"Product {product_id!r} has not been handed off.")
        return self._handles[product_id]

    def get(self, product_id, default=None):
        return self._handles.get(product_id, default)

    def validate_read(self, stage_id, product_id):
        stage = self.plan.stage_by_id[stage_id]
        if product_id not in (*stage.requires, *stage.optional_requires):
            raise ValueError(f"Stage {stage_id} did not declare input {product_id!r}.")

    def validate_handoff(self, stage_id, product_id):
        if self.plan.producer_by_product.get(product_id) != stage_id or self.has(product_id):
            raise ValueError(f"Undeclared or duplicate product handoff: {product_id!r}")

    def record_accepted(self, product_id, handle):
        if product_id not in self.plan.producer_by_product or self.has(product_id):
            raise ValueError(f"Undeclared or duplicate product: {product_id!r}")
        self._handles[product_id] = handle

    def record_failed(self, product_id, reason):
        if not self.has(product_id):
            self.failures[product_id] = str(reason)
