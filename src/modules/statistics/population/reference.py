"""Optional local reference annotation contract; no network implementation."""

from typing import Protocol

from src.core.models import Variant
from src.modules.statistics.population.models import ReferenceVariantAnnotation


class ReferenceProvider(Protocol):
    """Supply per-database evidence aligned to the specified reference and regions."""

    def annotate(
        self,
        variant: Variant,
        *,
        reference_sha256: str,
        region_definitions: dict[str, list[int]],
    ) -> dict[str, ReferenceVariantAnnotation]:
        """Return evidence; only observed=False establishes an explicit negative."""
        ...
