# candidate_summary

| category | group | parameter | value | count | minimum | q25 | median | q75 | maximum | unique_count | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| growth_interval | all | interval_start | all_parameters | 5420 | 7 | 5420 | 10705 | 17011.75 | 21039 | 1840 |  |
| growth_interval | all | interval_end | all_parameters | 5420 | 1038 | 5444.25 | 10764 | 17066 | 21083 | 1872 |  |
| growth_interval | by_window | window | 11 | 1765 | 7 | 6964 | 11938 | 17458 | 21030 | 606 |  |
| growth_interval | by_window | window | 21 | 1545 | 12 | 6692 | 11337 | 17169 | 21037 | 603 |  |
| growth_interval | by_window | window | 41 | 1175 | 22 | 5648 | 11182 | 17249 | 21036 | 505 |  |
| growth_interval | by_window | window | 81 | 642 | 41 | 4269.75 | 9922 | 16722.25 | 21039 | 301 |  |
| growth_interval | by_window | window | 161 | 177 | 80 | 1591 | 4221 | 7325 | 10563 | 56 |  |
| growth_interval | by_window | window | 321 | 116 | 160 | 208 | 5022 | 7366 | 10590 | 29 |  |
| growth_end | all_found | candidate_step | all_parameters | 59623 | 1339 | 5337 | 10774 | 16305 | 21090 | 679 | not_found rows retained separately: 5417 |
| growth_end | by_condition | end_condition | A | 19400 | 1416 | 5337 | 10777 | 16150 | 21090 | 243 |  |
| growth_end | by_condition | end_condition | B | 20823 | 1339 | 5334 | 10774 | 16614 | 21087 | 426 |  |
| growth_end | by_condition | end_condition | C | 19400 | 1407 | 5336 | 10651 | 16150 | 21089 | 283 |  |
| rank4_onset | all_found | candidate_step | all_parameters | 49 | 15 | 650 | 840 | 940 | 985 | 8 | persistence counts consecutive saved q records; no interpolation |
| rank4_onset | by_threshold | relative_threshold | 9.9999999999999995e-07 | 7 | 985 | 985 | 985 | 985 | 985 | 1 |  |
| rank4_onset | by_threshold | relative_threshold | 9.9999999999999995e-08 | 7 | 940 | 940 | 940 | 940 | 940 | 1 |  |
| rank4_onset | by_threshold | relative_threshold | 1e-08 | 7 | 265 | 885 | 885 | 885 | 885 | 2 |  |
| rank4_onset | by_threshold | relative_threshold | 1.0000000000000001e-09 | 7 | 15 | 562.5 | 840 | 840 | 840 | 4 |  |
| rank4_onset | by_threshold | relative_threshold | 1e-10 | 7 | 15 | 562.5 | 840 | 840 | 840 | 4 |  |
| rank4_onset | by_threshold | relative_threshold | 9.9999999999999994e-12 | 7 | 15 | 562.5 | 840 | 840 | 840 | 4 |  |
| rank4_onset | by_threshold | relative_threshold | 9.9999999999999998e-13 | 7 | 15 | 562.5 | 840 | 840 | 840 | 4 |  |
| time_difference | all_pairs | rank4_minus_growth_end | all_parameters | 2921527 | -21075 | -15620 | -9892 | -4497 | -354 | 4782 |  |
| time_difference | by_end_condition | end_condition | A | 950600 | -21075 | -15465 | -9837 | -4497 | -431 | 1849 |  |
| time_difference | by_end_condition | end_condition | B | 1020327 | -21072 | -15929 | -9973 | -4573 | -354 | 3155 |  |
| time_difference | by_end_condition | end_condition | C | 950600 | -21074 | -15465 | -9836 | -4496 | -422 | 2144 |  |
| growth_interval_start | recurring_exact_step | frequency_rank | 1 | 12 | 1560 | 1560 | 1560 | 1560 | 1560 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 2 | 10 | 14 | 14 | 14 | 14 | 14 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 3 | 10 | 25 | 25 | 25 | 25 | 25 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 4 | 10 | 160 | 160 | 160 | 160 | 160 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 5 | 10 | 1554 | 1554 | 1554 | 1554 | 1554 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 6 | 10 | 4353 | 4353 | 4353 | 4353 | 4353 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 7 | 9 | 4075 | 4075 | 4075 | 4075 | 4075 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 8 | 9 | 6350 | 6350 | 6350 | 6350 | 6350 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 9 | 9 | 9922 | 9922 | 9922 | 9922 | 9922 | 1 | 頻度順位であり採用順位ではない |
| growth_interval_start | recurring_exact_step | frequency_rank | 10 | 8 | 2051 | 2051 | 2051 | 2051 | 2051 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 1 | 1253 | 8161 | 8161 | 8161 | 8161 | 8161 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 2 | 1085 | 4591 | 4591 | 4591 | 4591 | 4591 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 3 | 812 | 1769 | 1769 | 1769 | 1769 | 1769 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 4 | 760 | 4592 | 4592 | 4592 | 4592 | 4592 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 5 | 747 | 3201 | 3201 | 3201 | 3201 | 3201 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 6 | 625 | 10916 | 10916 | 10916 | 10916 | 10916 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 7 | 560 | 2249 | 2249 | 2249 | 2249 | 2249 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 8 | 509 | 10777 | 10777 | 10777 | 10777 | 10777 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 9 | 483 | 6556 | 6556 | 6556 | 6556 | 6556 | 1 | 頻度順位であり採用順位ではない |
| growth_end | recurring_exact_step | frequency_rank | 10 | 480 | 8160 | 8160 | 8160 | 8160 | 8160 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 1 | 16 | 840 | 840 | 840 | 840 | 840 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 2 | 7 | 940 | 940 | 940 | 940 | 940 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 3 | 7 | 985 | 985 | 985 | 985 | 985 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 4 | 6 | 885 | 885 | 885 | 885 | 885 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 5 | 4 | 15 | 15 | 15 | 15 | 15 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 6 | 4 | 475 | 475 | 475 | 475 | 475 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 7 | 4 | 650 | 650 | 650 | 650 | 650 | 1 | 頻度順位であり採用順位ではない |
| rank4_onset | recurring_exact_step | frequency_rank | 8 | 1 | 265 | 265 | 265 | 265 | 265 | 1 | 頻度順位であり採用順位ではない |
