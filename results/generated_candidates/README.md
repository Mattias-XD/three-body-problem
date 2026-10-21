# Candidates from global normalizing-flow training

[Complete original candidate list](generated_ic_candidates_old_version.csv): all 97 saved model outputs from the global training run, reproduced without changing any decimal strings. Columns include `candidate_id`, `x`, `y`, `T` and `T_star`; trailing fields are the original verification-tracking fields.

Of these 97 entries, 94 passed the geometry filter used for the subsequent candidate plots and RK4 selection. IC_000002, IC_000054 and IC_000071 were excluded by that filter. The complete file retains them for provenance. Thirty candidates were selected for RK4 testing from the 94 eligible entries.

`T_star` is the CSV notation for $T^{*}=T|E|^{3/2}$, the scale-invariant period. Here, `T` is the period in the adopted simulation units and `E` is the total Newtonian energy. For fixed masses and gravitational constant, this quantity is unchanged under Newtonian rescaling of length and time. The table uses the original `T_star` predictions without rounding; these are not periods of confirmed periodic orbits.

[LaTeX table of the 30 selected candidates](norm_flow_original_candidates.tex), retaining every digit from the original source.
