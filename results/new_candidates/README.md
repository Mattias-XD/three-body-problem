# New-candidate results

[Download the current top ten from Models 1–4](top10_first_four_models.csv).

Updated 25 November 2026. Candidates are ranked together by the full-state return error from the completed, independent one-period confirmation after period refinement. This includes the earlier Region 5 RAW generation. The later T* model and the 100-coupling-layer model are excluded.

All listed integrations use dps = 50. Original integration and period search: atol = rtol = 1e-15; independent confirmation at the refined period: atol = rtol = 1e-18. Paused segments and two-period follow-ups are not ranked. CSV values retain their source precision; displayed values below use six decimal places.

| Rank | Model | Candidate | Before refinement | Period search | Confirmation |
|---:|---:|---|---:|---:|---:|
| 1 | 3 | NEW 000505 | 21.532203 | 0.110358 | 0.110358 |
| 2 | 4 | NEW 000942 | 3.129159 | 0.325972 | 0.325788 |
| 3 | 4 | NEW 000930 | 0.836169 | 0.327060 | 0.327060 |
| 4 | 3 | NEW 000508 | 2.813713 | 0.340353 | 0.340354 |
| 5 | 4 | NEW 000317 | 1.436239 | 0.344333 | 0.344333 |
| 6 | 2 | RAW 000006 | 0.674369 | 0.370753 | 0.370753 |
| 7 | 4 | NEW 000486 | 1.955769 | 0.372454 | 0.372451 |
| 8 | 2 | RAW 000002 | 0.527926 | 0.386030 | 0.386030 |
| 9 | 3 | NEW 000164 | 2.949080 | 0.397562 | 0.397454 |
| 10 | 2 | RAW 000051 | 2.115184 | 0.408841 | 0.408840 |

| Rank | Shape error (confirmation) | Position error (confirmation) | Velocity error (confirmation) | State error (confirmation) |
|---:|---:|---:|---:|---:|
| 1 | 0.038608 | 0.058560 | 0.093540 | 0.110358 |
| 2 | 0.200260 | 0.129547 | 0.298924 | 0.325788 |
| 3 | 0.206019 | 0.325206 | 0.034770 | 0.327060 |
| 4 | 0.049528 | 0.193802 | 0.279789 | 0.340354 |
| 5 | 0.357612 | 0.329802 | 0.098974 | 0.344333 |
| 6 | 0.251444 | 0.335518 | 0.157751 | 0.370753 |
| 7 | 0.303888 | 0.247053 | 0.278719 | 0.372451 |
| 8 | 0.093206 | 0.056206 | 0.381917 | 0.386030 |
| 9 | 0.335177 | 0.312329 | 0.245805 | 0.397454 |
| 10 | 0.466821 | 0.405058 | 0.055478 | 0.408840 |

- Model 1: global, float32, four coupling layers; no entry among the overall top ten.
- Model 2: Region 5, dps 100, four coupling layers.
- Model 3: global, dps 20, eight coupling layers.
- Model 4: global, dps 20, sixteen coupling layers.

The CSV includes full candidate identifiers, original coordinates and periods, refined periods, tolerances, error components and original/confirmation run IDs. Training dps above is distinct from the verification dps of 50. None of these ten candidates meets the periodicity threshold of 0.001.

## Historical results

[new_candidates_final_results.csv](new_candidates_final_results.csv) is the historical result list for the earlier 30-candidate experiment. It is retained separately and is not the current Models 1–4 ranking.
