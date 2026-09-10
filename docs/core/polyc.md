# Shared polyC Policy

`src/core/polyc.py` defines the shared, tool-agnostic directional polyC policy.
It classifies already-called candidates; it does not trim sequence, mutate
variants, parse trace names, or change serialized `Sample` / `Variant` data.

## HV2 directional coverage

The policy remains in effect when a trace extends beyond standard HV2 coverage:

| Direction | Excluded candidate anchors |
|---|---|
| `forward` | `304` and every later position |
| `reverse` | `315` and every earlier position |

`HV2_POLYC_REGION` remains `303–315` for canonical insertion normalization.
Reverse exclusion includes that tract so raw calls cannot compete with its
canonical insertions. Normalization is Tracy-specific; the shared policy decides
only which already-called candidates enter trusted per-trace evidence.

## HV1 16189-created polyC

A called `C` at position anchor `16189` creates the shared HV1 polyC condition.
The policy deliberately uses the anchor and called base only; callers do not add
a reference-base requirement. When that condition is present:

| Direction | Suppressed candidates |
|---|---|
| `forward` | Positions strictly above `16189` |
| `reverse` | Positions strictly below `16189` |

Position `16189` itself and unknown direction remain usable. Tracy's separate
16189-deletion conversion and synthetic 16193 deletion are caller-specific
nomenclature transforms, not shared policy.

## API

```python
from src.core.polyc import (
    directional_hv1_polyc_suppression_reason,
    directional_hv2_polyc_exclusion_reason,
    is_hv1_polyc_created,
)
```

The suppression helpers return a centralized reason or `None`. Decimal
insertion positions are classified by their integer anchor through the shared
`pos_base()` helper. Tools map known primer provenance to `"forward"` or
`"reverse"`, then apply the result to their own candidate lifecycle. This keeps
all policy boundaries, triggers, and reasons authoritative in one module while
leaving tool-specific diagnostics and future quality-aware recovery independent.
