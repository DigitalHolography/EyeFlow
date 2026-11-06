"""Logical identities do not depend on HDF5 paths or artifact filenames."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PipelineInvocationId:
    run_id: str
    pipeline: str
    variant: str | None = None

    @property
    def label(self):
        return self.pipeline if self.variant is None else f"{self.pipeline} [{self.variant}]"


@dataclass(frozen=True, slots=True)
class ProductId:
    producer: PipelineInvocationId
    name: str


@dataclass(frozen=True, slots=True, kw_only=True)
class ProductRequirement:
    alias: str
    producer_pipeline: str
    product: str
    variant: str = "shared"

    def __post_init__(self):
        if any(not isinstance(value, str) or not value.strip()
               for value in (self.alias, self.producer_pipeline, self.product, self.variant)):
            raise ValueError("Product requirements need nonempty alias, producer, product and variant scope.")

    def resolve(self, consumer):
        variant = (consumer.variant if self.variant == "same" else
                   None if self.variant == "shared" else self.variant)
        return ProductId(PipelineInvocationId(consumer.run_id, self.producer_pipeline, variant), self.product)
