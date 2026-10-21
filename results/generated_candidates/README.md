# Candidates from global normalizing-flow training

[Complete original candidate list](generated_ic_candidates_old_version.csv): all 97 saved model outputs from the global training run, reproduced without changing any decimal strings. Columns include `candidate_id`, `x`, `y`, `T` and `T_star`; trailing fields are the original verification-tracking fields.

Of these 97 entries, 94 passed the geometry filter used for the subsequent candidate plots and RK4 selection. IC_000002, IC_000054 and IC_000071 were excluded by that filter. The complete file retains them for provenance. Thirty candidates were selected for RK4 testing from the 94 eligible entries.

`T_star` denotes the scale-invariant period, while `T` is the physical period in the adopted simulation units. These are model predictions, not confirmed periodic orbits.

[LaTeX table of the 30 selected candidates](norm_flow_original_candidates.tex), retaining every digit from the original source.
