# Numerical verification results

- [Original initial conditions](original_initial_conditions/verification_original_ic_completed.csv): 31 completed catalogue-IC verifications, one final row per IC.
- [New candidates](new_candidates/new_candidates_final_results.csv): 30 candidates, consisting of 27 completed integrations and the latest partial result for three unfinished candidates (IC_000040, IC_000079 and IC_000094).

For each IC, the latest completed result is retained, including the higher-precision period-confirmation run where available. Intermediate pause/resume records are omitted. For the three unfinished candidates, the latest partial record is retained; its return error is not a full-period result.

The CSV files omit monitored conserved-quantity diagnostics (energy, momentum, angular momentum and centre-of-mass diagnostics, where present) and the accumulated local-error diagnostic. All retained values preserve their original precision, column names and recorded values, including historical machine paths and Swedish classification labels (`periodisk` = periodic; `icke-periodisk` = non-periodic). Those paths refer to the machines that produced the results and are not repository links. `classification_2` requires all four return errors to be strictly below `threshold`.

[Original-resolution catalogue verification plots](../verification/original_initial_conditions/).

## Report tables

Candidate-table state errors are displayed with six decimal places, without scientific notation.

- [Catalogue verification (LaTeX)](original_initial_conditions/catalogue_verification.tex): all 31 completed catalogue cases, in catalogue-row order.
- [Period-refined candidates (LaTeX)](new_candidates/best_candidates_period_refinement.tex): only the 20 candidates with a refined return time and a separate higher-precision confirmation integration, ranked by confirmed state error. Boundary-minimum cases and incomplete original integrations are excluded.

Candidate period refinement keeps initial positions fixed and searches 33 times in `[0.9 T0, 1.1 T0]`, followed by a golden-section search around the best interior grid point until the bracket width is at most `1e-12 T0`. The improved return time is checked by a new integration from the initial state, with 70 decimal digits and `atol = rtol = 1e-26` (original runs: 50 digits and `1e-24`). The displayed refined errors are from these confirmation integrations.

- [All raw candidates (LaTeX)](new_candidates/best_raw_candidates.tex): all 30 candidates before refinement; completed cases ranked by original-period state error and incomplete cases shown last in red.
