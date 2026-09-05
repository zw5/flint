# Flint phenotype atlas

**291 recorded phenotype fields · 1,164 label-by-dimension evaluations · five domains**

This is the complete retained domain-recovery result table, rendered for browsing. Every label and all four tested coordinate counts are included. Fields can be individual questionnaire items, repeated measurements, or alternate units of the same assay; 291 fields are not 291 independent biological discoveries.

These values are Pearson correlations between pooled leave-one-out standardized predictions and standardized outcomes. The broader identity-coordinate pipeline refits the EEG reliability basis and phenotype map within each fold. It does not apply the nuisance adjustment used in the primary fluid-reasoning result.

**Reading the columns:** `Observed n` counts original nonmissing labels in the 111-person cohort. Missing standardized targets are mean-imputed in the original scoring; the correlation therefore includes imputed outcomes. K is the number of retained identity coordinates. A negative prediction correlation is an inverted prediction relationship, not evidence that a particular EEG feature is a protective or adverse biological factor. The values have no per-label significance or multiplicity correction attached.

This catalog is generated from retained aggregate outputs. Their underlying participant predictions have not been freshly recomputed. The statistical verifier for the primary fluid-reasoning tables covers a separate experiment.

[Research note](paper.md) · [Full original label table](results/context/phenotypes/domain_label_reconstruction.csv) · [Full original domain table](results/context/phenotypes/domain_recovery_summary.csv)

## Domain-level performance

Pooled standardized leave-one-out R² across every label in each domain. These values assess the full domain reconstruction, rather than its most favorable individual row. Negative R² means the prediction errors exceed the globally centered mean baseline used by the original calculation.

| Domain | Fields | K = 3 | K = 6 | K = 10 | K = 20 |
|---|---:|---:|---:|---:|---:|
| Blood chemistry | 35 | +0.0098 | -0.0105 | -0.0399 | -0.1217 |
| Blood pressure | 11 | +0.1263 | +0.1008 | +0.0741 | +0.0393 |
| Body measurements | 4 | +0.1192 | +0.0871 | +0.0425 | +0.0080 |
| Cognition | 101 | -0.0074 | -0.0535 | -0.0796 | -0.1282 |
| Personality and affect | 140 | -0.0089 | -0.0500 | -0.0688 | -0.1394 |

## All label-level results

Labels appear in source-name order within each domain. Display names retain the final source group and exact field key; the linked CSV preserves complete original paths. No ordering or filtering uses the correlation magnitude.

### Blood chemistry · 35 fields

| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |
|---|---:|---:|---:|---:|---:|
| Blood_Results_LEMON / ALAT_in_ µkat_l | 108 | +0.0957 | +0.2851 | +0.1908 | +0.0681 |
| Blood_Results_LEMON / ASAT_in_µkat_l | 108 | +0.0431 | +0.3177 | +0.2844 | +0.2021 |
| Blood_Results_LEMON / CHOL_in_mmol_l | 108 | +0.1741 | +0.1001 | +0.1539 | +0.1768 |
| Blood_Results_LEMON / CKDEPI_in_ml_min_1.73m | 90 | +0.4152 | +0.4344 | +0.5090 | +0.4833 |
| Blood_Results_LEMON / CL_in_mmol_l | 108 | -0.0286 | +0.1247 | +0.1516 | +0.0121 |
| Blood_Results_LEMON / CRE_in_µmol_l | 108 | -0.0151 | +0.0196 | +0.0951 | +0.0948 |
| Blood_Results_LEMON / CRP_in_mg_l | 87 | +0.1881 | +0.2464 | +0.1316 | +0.1296 |
| Blood_Results_LEMON / GGT_in_µkat_l | 108 | +0.2068 | +0.2427 | +0.3042 | +0.3532 |
| Blood_Results_LEMON / GLU_in_mmol_l | 108 | +0.2245 | +0.1559 | +0.1422 | +0.1612 |
| Blood_Results_LEMON / HBA1CI_in_mmol_mol | 109 | +0.3690 | +0.3449 | +0.3191 | +0.2857 |
| Blood_Results_LEMON / HBA1C_in_% | 109 | +0.3691 | +0.3450 | +0.3193 | +0.2859 |
| Blood_Results_LEMON / HBA1G_in_g_dl | 109 | +0.2214 | +0.1903 | +0.2178 | +0.2586 |
| Blood_Results_LEMON / HB_in_g_dl | 109 | -0.0441 | -0.0379 | +0.0233 | +0.0512 |
| Blood_Results_LEMON / HCT _in_l_l | 107 | +0.0479 | -0.0606 | +0.0181 | -0.0286 |
| Blood_Results_LEMON / HDLC_in_mmol_l | 108 | +0.1924 | +0.0201 | +0.0221 | -0.1198 |
| Blood_Results_LEMON / HGBK_in_g_dl | 107 | +0.0644 | +0.0113 | +0.0175 | +0.0459 |
| Blood_Results_LEMON / HGB_in_mmol_l | 107 | +0.0624 | +0.0082 | +0.0168 | +0.0504 |
| Blood_Results_LEMON / INR | 104 | +0.1427 | +0.2347 | +0.2048 | +0.1845 |
| Blood_Results_LEMON / INRVOR | 104 | +0.1191 | +0.2315 | +0.2253 | +0.2099 |
| Blood_Results_LEMON / K_in_mmol_l | 108 | -0.0271 | -0.1708 | -0.1843 | -0.0848 |
| Blood_Results_LEMON / LDLC_in_mmol_l | 108 | +0.2108 | +0.1966 | +0.2162 | +0.2512 |
| Blood_Results_LEMON / MCHC _in_mmol_l | 107 | -0.1424 | +0.0912 | +0.0241 | +0.0178 |
| Blood_Results_LEMON / MCHCK_in_g_dl | 107 | -0.1461 | +0.0947 | +0.0259 | +0.0174 |
| Blood_Results_LEMON / MCHK_in_pg | 107 | +0.0986 | +0.1379 | +0.1128 | -0.0336 |
| Blood_Results_LEMON / MCH_in_fmol | 107 | +0.0972 | +0.1347 | +0.1091 | -0.0359 |
| Blood_Results_LEMON / MCV_in_fl | 107 | +0.0908 | +0.0933 | +0.0581 | -0.0913 |
| Blood_Results_LEMON / MPV_in_fl | 105 | +0.0257 | +0.0158 | +0.0848 | -0.0115 |
| Blood_Results_LEMON / NA_in_mmol_l | 108 | +0.1910 | +0.1453 | +0.0907 | +0.0097 |
| Blood_Results_LEMON / PLT_in_exp_9 _l | 107 | -0.2030 | -0.2298 | -0.2384 | -0.1197 |
| Blood_Results_LEMON / PT_in_% | 104 | +0.0922 | +0.2566 | +0.2658 | +0.2270 |
| Blood_Results_LEMON / RBC_in_exp_12_l | 107 | -0.0452 | -0.2029 | -0.0803 | -0.1079 |
| Blood_Results_LEMON / RDW_in_% | 107 | +0.1291 | +0.1242 | +0.1001 | +0.1169 |
| Blood_Results_LEMON / TRIG_in_mmol_l | 108 | +0.0455 | -0.0832 | -0.2160 | -0.0728 |
| Blood_Results_LEMON / TSH_in_mU_l | 108 | -0.4082 | -0.4539 | -0.0588 | -0.1195 |
| Blood_Results_LEMON / WBC_in_exp_9_l | 107 | -0.0577 | -0.1188 | -0.1210 | +0.0215 |

### Blood pressure · 11 fields

| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |
|---|---:|---:|---:|---:|---:|
| Blood_Pressure_LEMON / BP1_left_diastole | 110 | +0.3807 | +0.3918 | +0.3756 | +0.3763 |
| Blood_Pressure_LEMON / BP1_left_systole | 110 | +0.3904 | +0.3637 | +0.3982 | +0.4069 |
| Blood_Pressure_LEMON / BP1_right_diastole | 107 | +0.4195 | +0.3819 | +0.2970 | +0.4014 |
| Blood_Pressure_LEMON / BP1_right_systole | 107 | +0.4379 | +0.4205 | +0.3720 | +0.3455 |
| Blood_Pressure_LEMON / BP2_left_diastole | 109 | +0.3935 | +0.3728 | +0.3541 | +0.3048 |
| Blood_Pressure_LEMON / BP2_left_systole | 109 | +0.5027 | +0.4721 | +0.4537 | +0.4546 |
| Blood_Pressure_LEMON / BP3_left_diastole | 80 | +0.3919 | +0.3539 | +0.3337 | +0.2961 |
| Blood_Pressure_LEMON / BP3_left_systole | 80 | +0.2835 | +0.2840 | +0.2665 | +0.1592 |
| Blood_Pressure_LEMON / pulse1_left | 108 | +0.2344 | +0.1957 | +0.1586 | +0.2120 |
| Blood_Pressure_LEMON / pulse1_right | 106 | +0.2013 | +0.1221 | +0.1392 | +0.2383 |
| Blood_Pressure_LEMON / pulse2_left | 109 | +0.0883 | +0.0825 | +0.1251 | +0.1265 |

### Body measurements · 4 fields

| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |
|---|---:|---:|---:|---:|---:|
| Anthropometry_LEMON / Height_cm | 110 | -0.0655 | -0.0005 | -0.0587 | +0.1867 |
| Anthropometry_LEMON / Hip_cm | 110 | +0.3804 | +0.3079 | +0.2856 | +0.2456 |
| Anthropometry_LEMON / Waist_cm | 110 | +0.5133 | +0.4941 | +0.4786 | +0.5031 |
| Anthropometry_LEMON / Weight_kg | 110 | +0.3144 | +0.3186 | +0.2746 | +0.1694 |

### Cognition · 101 fields

| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |
|---|---:|---:|---:|---:|---:|
| _CVLT / CVLT_1 | 111 | +0.0028 | -0.0032 | -0.0032 | +0.0066 |
| _CVLT / CVLT_10 | 110 | +0.0164 | -0.0150 | -0.0940 | -0.0173 |
| _CVLT / CVLT_11 | 110 | +0.1094 | +0.1096 | +0.0209 | +0.0502 |
| _CVLT / CVLT_12 | 109 | +0.0686 | +0.1206 | +0.0457 | +0.0720 |
| _CVLT / CVLT_13 | 110 | -0.1129 | -0.2135 | -0.2596 | -0.1335 |
| _CVLT / CVLT_14 | 110 | -0.2236 | +0.0391 | -0.0555 | -0.1127 |
| _CVLT / CVLT_15 | 110 | +0.0697 | +0.0384 | +0.0352 | +0.1053 |
| _CVLT / CVLT_2 | 111 | +0.0672 | -0.0453 | -0.0377 | -0.0085 |
| _CVLT / CVLT_3 | 110 | +0.0222 | -0.0253 | +0.0113 | +0.0479 |
| _CVLT / CVLT_4 | 111 | +0.0235 | +0.0791 | +0.0459 | -0.0142 |
| _CVLT / CVLT_5 | 110 | +0.0800 | +0.0608 | +0.0393 | -0.1071 |
| _CVLT / CVLT_6 | 110 | +0.1963 | +0.1153 | +0.1194 | +0.0884 |
| _CVLT / CVLT_7 | 78 | -0.1157 | -0.0778 | -0.1056 | -0.1527 |
| _CVLT / CVLT_8 | 111 | +0.1955 | +0.0608 | +0.0034 | +0.0633 |
| _CVLT / CVLT_9 | 111 | +0.1535 | +0.1218 | +0.1196 | +0.1899 |
| LPS / LPS_1 | 111 | +0.4125 | +0.4426 | +0.4224 | +0.4589 |
| RWT / RWT_1 | 110 | +0.1117 | -0.0091 | -0.0683 | -0.1559 |
| RWT / RWT_10 | 110 | +0.0478 | +0.1177 | +0.0884 | +0.0877 |
| RWT / RWT_11 | 110 | -0.0492 | -0.0640 | -0.0452 | -0.0037 |
| RWT / RWT_13 | 110 | +0.1194 | +0.0596 | +0.1685 | +0.1757 |
| RWT / RWT_14 | 103 | +0.0243 | +0.0282 | +0.0968 | -0.0212 |
| RWT / RWT_15 | 110 | +0.2071 | +0.2344 | +0.1810 | +0.2459 |
| RWT / RWT_16 | 110 | -0.5656 | -0.2268 | -0.0037 | -0.0232 |
| RWT / RWT_17 | 110 | -0.0102 | +0.0717 | +0.1185 | +0.2431 |
| RWT / RWT_18 | 110 | -0.1069 | -0.0002 | -0.0943 | -0.0286 |
| RWT / RWT_19 | 110 | -0.0269 | -0.0117 | +0.1122 | +0.0895 |
| RWT / RWT_2 | 110 | -0.2075 | -0.1725 | -0.1851 | -0.1752 |
| RWT / RWT_20 | 110 | +0.0859 | +0.0457 | +0.1182 | +0.2238 |
| RWT / RWT_21 | 104 | +0.1077 | +0.0713 | +0.1100 | +0.0926 |
| RWT / RWT_22 | 110 | +0.1536 | +0.2050 | +0.1787 | +0.1533 |
| RWT / RWT_23 | 110 | -0.3821 | -0.1389 | +0.0533 | +0.0264 |
| RWT / RWT_3 | 110 | -0.0796 | +0.0540 | +0.0531 | -0.1385 |
| RWT / RWT_4 | 110 | +0.0399 | -0.0400 | -0.1352 | +0.0366 |
| RWT / RWT_5 | 110 | +0.0786 | -0.1100 | -0.1035 | -0.1016 |
| RWT / RWT_6 | 110 | +0.1454 | +0.0826 | +0.0256 | +0.1056 |
| RWT / RWT_7 | 110 | -0.1886 | -0.0899 | -0.0143 | -0.0277 |
| RWT / RWT_8 | 110 | +0.1244 | -0.0304 | -0.1025 | -0.1387 |
| RWT / RWT_9 | 110 | -0.0375 | -0.1269 | -0.1798 | -0.1776 |
| TAP-Alertness / TAP_A_1 | 111 | +0.1402 | +0.0061 | +0.0643 | +0.0885 |
| TAP-Alertness / TAP_A_10 | 111 | +0.1440 | +0.0252 | +0.0732 | +0.1317 |
| TAP-Alertness / TAP_A_11 | 111 | +0.0939 | -0.0401 | +0.0040 | +0.0618 |
| TAP-Alertness / TAP_A_12 | 110 | -0.0095 | -0.0269 | +0.0064 | -0.0535 |
| TAP-Alertness / TAP_A_13 | 111 | +0.1562 | +0.0770 | +0.0592 | +0.1837 |
| TAP-Alertness / TAP_A_14 | 109 | +0.0631 | +0.0961 | +0.0108 | +0.0705 |
| TAP-Alertness / TAP_A_15 | 111 | -0.0649 | -0.1575 | +0.0087 | -0.1124 |
| TAP-Alertness / TAP_A_16 | 110 | -0.1354 | -0.2108 | -0.0623 | -0.1469 |
| TAP-Alertness / TAP_A_2 | 111 | +0.1224 | +0.0160 | +0.0886 | +0.0798 |
| TAP-Alertness / TAP_A_3 | 111 | +0.0913 | -0.0299 | -0.0701 | +0.0559 |
| TAP-Alertness / TAP_A_4 | 111 | +0.1101 | +0.0530 | +0.0416 | +0.1358 |
| TAP-Alertness / TAP_A_5 | 111 | +0.1475 | +0.0138 | +0.0499 | +0.1477 |
| TAP-Alertness / TAP_A_6 | 111 | +0.1711 | +0.0493 | +0.0737 | +0.1344 |
| TAP-Alertness / TAP_A_7 | 110 | -0.0381 | -0.1040 | -0.0729 | -0.1210 |
| TAP-Alertness / TAP_A_8 | 111 | -0.0332 | -0.1302 | -0.1380 | +0.1132 |
| TAP-Alertness / TAP_A_9 | 107 | -0.1625 | -0.0679 | -0.0956 | -0.0624 |
| TAP-Incompatibility / TAP_I_1 | 111 | +0.1207 | +0.0256 | +0.2832 | +0.2598 |
| TAP-Incompatibility / TAP_I_10 | 110 | -0.1719 | -0.3055 | -0.0841 | -0.2640 |
| TAP-Incompatibility / TAP_I_11 | 111 | +0.1721 | +0.0347 | +0.0669 | +0.0479 |
| TAP-Incompatibility / TAP_I_12 | 109 | -0.6088 | -0.2979 | -0.2306 | -0.1695 |
| TAP-Incompatibility / TAP_I_13 | 111 | -0.6568 | -0.3167 | -0.4099 | -0.1230 |
| TAP-Incompatibility / TAP_I_14 | 108 | +0.0165 | -0.0639 | -0.1644 | -0.3153 |
| TAP-Incompatibility / TAP_I_15 | 111 | +0.2601 | +0.1701 | +0.3487 | +0.3442 |
| TAP-Incompatibility / TAP_I_16 | 111 | +0.2563 | +0.1750 | +0.3686 | +0.3632 |
| TAP-Incompatibility / TAP_I_17 | 110 | -0.2683 | -0.1250 | +0.0894 | -0.1462 |
| TAP-Incompatibility / TAP_I_18 | 111 | +0.2037 | +0.0645 | +0.1347 | +0.1532 |
| TAP-Incompatibility / TAP_I_19 | 110 | -0.2144 | -0.1061 | -0.1287 | -0.1878 |
| TAP-Incompatibility / TAP_I_2 | 111 | +0.1380 | +0.0169 | +0.2442 | +0.2314 |
| TAP-Incompatibility / TAP_I_20 | 111 | -0.4913 | -0.2098 | -0.2909 | -0.1418 |
| TAP-Incompatibility / TAP_I_21 | 107 | +0.0745 | -0.0094 | -0.0612 | -0.2369 |
| TAP-Incompatibility / TAP_I_22 | 109 | -0.1233 | -0.1141 | -0.1320 | -0.3137 |
| TAP-Incompatibility / TAP_I_23 | 107 | -0.1768 | -0.0954 | -0.1581 | -0.1878 |
| TAP-Incompatibility / TAP_I_24 | 110 | +0.2590 | +0.2327 | +0.2186 | +0.1245 |
| TAP-Incompatibility / TAP_I_25 | 107 | +0.0550 | +0.0863 | +0.0635 | +0.1532 |
| TAP-Incompatibility / TAP_I_26 | 110 | +0.1076 | +0.1766 | +0.1704 | +0.1265 |
| TAP-Incompatibility / TAP_I_27 | 109 | -0.0707 | -0.0343 | +0.0099 | -0.0766 |
| TAP-Incompatibility / TAP_I_3 | 109 | -0.1344 | +0.0922 | +0.1247 | -0.0699 |
| TAP-Incompatibility / TAP_I_4 | 111 | +0.1260 | +0.0204 | +0.1366 | +0.0761 |
| TAP-Incompatibility / TAP_I_5 | 110 | -0.0381 | +0.0664 | +0.0390 | -0.1513 |
| TAP-Incompatibility / TAP_I_6 | 111 | -0.3193 | -0.0770 | -0.0978 | -0.1474 |
| TAP-Incompatibility / TAP_I_8 | 111 | +0.2659 | +0.2007 | +0.3125 | +0.3569 |
| TAP-Incompatibility / TAP_I_9 | 111 | +0.2658 | +0.2212 | +0.2892 | +0.3752 |
| TAP-Working_Memory / TAP_WM_1 | 111 | -0.2798 | +0.0275 | +0.0947 | +0.0240 |
| TAP-Working_Memory / TAP_WM_10 | 110 | -0.5898 | +0.0062 | +0.0343 | +0.1080 |
| TAP-Working_Memory / TAP_WM_11 | 111 | -0.2492 | -0.1141 | -0.0187 | -0.0428 |
| TAP-Working_Memory / TAP_WM_2 | 111 | -0.1420 | +0.0332 | +0.0415 | +0.0189 |
| TAP-Working_Memory / TAP_WM_3 | 109 | -0.0191 | +0.0290 | +0.0121 | -0.2447 |
| TAP-Working_Memory / TAP_WM_4 | 111 | -0.0539 | -0.1951 | -0.0700 | -0.1668 |
| TAP-Working_Memory / TAP_WM_5 | 109 | -0.1758 | -0.2187 | -0.0373 | -0.1766 |
| TAP-Working_Memory / TAP_WM_6 | 111 | +0.1312 | +0.1456 | +0.1457 | +0.2308 |
| TAP-Working_Memory / TAP_WM_7 | 111 | -0.0999 | -0.0777 | -0.0368 | +0.1579 |
| TAP-Working_Memory / TAP_WM_8 | 65 | -0.0375 | -0.0844 | -0.0757 | -0.0202 |
| TAP-Working_Memory / TAP_WM_9 | 111 | +0.1324 | +0.1475 | +0.1429 | +0.2272 |
| TMT / TMT_1 | 111 | +0.3383 | +0.2763 | +0.3162 | +0.3414 |
| TMT / TMT_2 | 111 | +0.3095 | +0.2584 | +0.2577 | +0.3315 |
| TMT / TMT_3 | 111 | +0.1260 | +0.0123 | -0.0658 | -0.0393 |
| TMT / TMT_5 | 111 | +0.2407 | +0.2006 | +0.2028 | +0.3387 |
| TMT / TMT_6 | 111 | +0.2059 | +0.1167 | +0.2053 | +0.3374 |
| TMT / TMT_7 | 111 | +0.0739 | -0.0613 | +0.0755 | +0.1021 |
| WST / WST_1 | 111 | -0.0056 | +0.0269 | +0.1071 | +0.1164 |
| WST / WST_2 | 111 | -0.0280 | +0.0244 | +0.1339 | +0.1257 |
| WST / WST_3 | 111 | -0.0269 | +0.0202 | +0.1303 | +0.1271 |
| WST / WST_4 | 111 | -0.0342 | +0.0216 | +0.1367 | +0.1259 |

### Personality and affect · 140 fields

| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |
|---|---:|---:|---:|---:|---:|
| BISBAS / BAS_Drive | 111 | +0.1947 | +0.1263 | +0.0922 | -0.1164 |
| BISBAS / BAS_Fun | 111 | +0.0439 | +0.0211 | -0.0093 | +0.0322 |
| BISBAS / BAS_Reward | 111 | +0.1653 | +0.1889 | +0.1561 | -0.0537 |
| BISBAS / BIS | 111 | +0.1483 | +0.0477 | +0.0538 | +0.0128 |
| CERQ / CERQ_Acceptance | 111 | -0.0309 | +0.0955 | +0.0604 | +0.0129 |
| CERQ / CERQ_BlamingOthers | 111 | -0.3028 | +0.0146 | -0.0174 | -0.1114 |
| CERQ / CERQ_Catastrophizing | 111 | -0.3495 | -0.0803 | -0.1532 | -0.1353 |
| CERQ / CERQ_PositiveReappraisal | 111 | -0.1213 | -0.3144 | -0.0414 | -0.0881 |
| CERQ / CERQ_PuttingIntoPerspective | 111 | -0.0398 | -0.0379 | -0.0617 | -0.1662 |
| CERQ / CERQ_RefocusOnPlanning | 111 | -0.1788 | -0.1040 | -0.0832 | -0.1862 |
| CERQ / CERQ_Rumination | 111 | -0.2119 | -0.2421 | -0.0184 | -0.1632 |
| CERQ / CERQ_SelfBlame | 111 | -0.0012 | -0.0057 | +0.0229 | -0.0544 |
| CERQ / CERQ_positiveRefocusing | 111 | -0.0308 | -0.1733 | -0.1556 | +0.0842 |
| COPE / COPE_Acceptance | 111 | +0.0074 | -0.1275 | +0.0394 | +0.0032 |
| COPE / COPE_Alkohol‎_Drogen | 111 | -0.1397 | -0.2414 | -0.2408 | -0.3304 |
| COPE / COPE_BehavioralDisengagement | 111 | +0.0234 | -0.0737 | -0.1119 | -0.1042 |
| COPE / COPE_Denial | 111 | -0.0970 | +0.0018 | -0.0372 | -0.1082 |
| COPE / COPE_Humor | 111 | +0.1017 | +0.0859 | +0.0721 | -0.0616 |
| COPE / COPE_Planning | 111 | -0.1160 | -0.2818 | -0.3108 | -0.3424 |
| COPE / COPE_Religion | 111 | +0.0208 | +0.0351 | +0.0141 | +0.0597 |
| COPE / COPE_SelfBlame | 111 | -0.3022 | +0.0856 | -0.0239 | -0.0302 |
| COPE / COPE_SelfDistraction | 111 | -0.1275 | +0.0087 | +0.0871 | +0.1322 |
| COPE / COPE_UseOfEmotionalSupport | 111 | +0.2409 | +0.1914 | +0.1798 | +0.1905 |
| COPE / COPE_UseOfInstrumentalSupport | 111 | +0.2037 | +0.1586 | +0.1721 | +0.1185 |
| COPE / COPE_Venting | 111 | -0.2122 | +0.0148 | +0.0074 | +0.1128 |
| COPE / COPE_activeCoping | 111 | -0.0094 | -0.0024 | -0.0645 | -0.1766 |
| COPE / COPE_positiveReframing | 111 | -0.0835 | -0.1265 | +0.0745 | +0.0149 |
| ERQ / ERQ_reappraisal | 111 | +0.0421 | +0.0607 | -0.0482 | +0.1907 |
| ERQ / ERQ_suppression | 111 | +0.0586 | -0.0084 | -0.1652 | -0.1293 |
| F-SozU_K-22 / FSoZu_EU | 111 | -0.1036 | -0.2224 | -0.1575 | -0.2952 |
| F-SozU_K-22 / FSoZu_PU | 111 | +0.0382 | -0.1250 | +0.0358 | -0.0962 |
| F-SozU_K-22 / FSoZu_SI | 111 | -0.0138 | -0.1412 | +0.0112 | -0.1129 |
| F-SozU_K-22 / FSoZu_Vert | 111 | -0.1749 | -0.1625 | -0.0895 | -0.1233 |
| F-SozU_K-22 / FSoZu_Zuf | 111 | +0.0447 | +0.0136 | -0.0225 | -0.0579 |
| FEV / FEV_HUNGER | 111 | +0.2293 | +0.2194 | +0.2705 | +0.1902 |
| FEV / FEV_KK | 111 | +0.1477 | +0.1292 | +0.1603 | +0.0698 |
| FEV / FEV_STOER | 111 | +0.0297 | +0.0297 | +0.0517 | +0.0299 |
| FTP / FTP_SUM | 70 | +0.2987 | +0.2726 | +0.3148 | +0.2885 |
| LOT-R / LOT_Optimism | 111 | +0.3176 | +0.3015 | +0.3782 | +0.3586 |
| LOT-R / LOT_Pessimism | 111 | -0.0031 | +0.0082 | +0.1039 | +0.1204 |
| LOT-R / LOT_sumscore | 111 | +0.1723 | +0.1796 | +0.2160 | +0.1569 |
| MARS / MARS_Affect_focused | 103 | -0.0407 | +0.1805 | +0.2716 | +0.0892 |
| MARS / MARS_Avoidance | 103 | +0.0009 | -0.0853 | +0.0088 | -0.0191 |
| MARS / MARS_Behavioral_Distraction | 103 | +0.1749 | +0.2396 | +0.2623 | +0.2465 |
| MARS / MARS_Cognitive_Distraction | 103 | -0.0376 | -0.0450 | -0.0759 | +0.0135 |
| MARS / MARS_Disengagement | 103 | -0.0226 | -0.0878 | -0.0298 | -0.0503 |
| MARS / MARS_Situation_focused | 103 | +0.1110 | +0.0392 | +0.0109 | +0.0175 |
| MDBF_Day1 / MDBF_Day1_GS_Scale | 109 | +0.1826 | +0.2227 | +0.1507 | +0.0652 |
| MDBF_Day1 / MDBF_Day1_RU_Scale | 108 | +0.1196 | +0.0385 | +0.0213 | -0.0381 |
| MDBF_Day1 / MDBF_Day1_WM_Scale | 110 | -0.0340 | -0.0274 | -0.1197 | -0.2047 |
| MDBF_Day2 / MDBF_Day2_GS_Scale | 109 | -0.0474 | -0.0035 | +0.0613 | -0.0021 |
| MDBF_Day2 / MDBF_Day2_RU_Scale | 109 | -0.1328 | +0.0187 | -0.0071 | +0.0311 |
| MDBF_Day2 / MDBF_Day2_WM_Scale | 108 | +0.0723 | +0.0725 | -0.0631 | -0.0707 |
| MDBF_Day3 / MDBF_Day3_GS_Scale | 80 | +0.1771 | +0.1585 | +0.0933 | -0.0050 |
| MDBF_Day3 / MDBF_Day3_RU_Scale | 80 | -0.0327 | -0.0089 | -0.0589 | -0.1394 |
| MDBF_Day3 / MDBF_Day3_WM_Scale | 79 | +0.2176 | +0.2085 | +0.2293 | +0.0746 |
| MSPSS / MSPSS_Family | 111 | -0.0645 | -0.0095 | +0.0840 | -0.0296 |
| MSPSS / MSPSS_Friends | 111 | +0.1546 | +0.0892 | +0.1845 | +0.0394 |
| MSPSS / MSPSS_SignificantOthers | 111 | -0.2395 | -0.2672 | -0.2546 | -0.2243 |
| MSPSS / MSPSS_total | 111 | -0.2675 | -0.1232 | +0.0665 | -0.0915 |
| NEO_FFI / NEOFFI_Agreeableness | 111 | -0.0239 | +0.0845 | +0.0617 | -0.0484 |
| NEO_FFI / NEOFFI_Conscientiousness | 111 | +0.2484 | +0.1874 | +0.1316 | +0.0441 |
| NEO_FFI / NEOFFI_Extraversion | 111 | +0.0639 | +0.0865 | +0.0418 | -0.0891 |
| NEO_FFI / NEOFFI_Neuroticism | 111 | +0.1345 | +0.0516 | -0.0004 | +0.0666 |
| NEO_FFI / NEOFFI_OpennessForExperiences | 111 | +0.0304 | -0.0310 | -0.0476 | -0.0262 |
| NYC_Q_lemon / NYC-Q_lemon_1 | 110 | +0.1546 | +0.1383 | +0.0679 | +0.0942 |
| NYC_Q_lemon / NYC-Q_lemon_10 | 111 | +0.0072 | -0.0532 | -0.1218 | -0.0919 |
| NYC_Q_lemon / NYC-Q_lemon_11 | 111 | -0.2384 | +0.0792 | +0.0856 | +0.0707 |
| NYC_Q_lemon / NYC-Q_lemon_12 | 108 | -0.0864 | -0.1779 | +0.0543 | -0.1095 |
| NYC_Q_lemon / NYC-Q_lemon_13 | 110 | -0.0820 | -0.2288 | -0.0043 | -0.1423 |
| NYC_Q_lemon / NYC-Q_lemon_14 | 111 | +0.0529 | +0.0520 | +0.0322 | +0.0188 |
| NYC_Q_lemon / NYC-Q_lemon_15 | 110 | +0.0550 | +0.1336 | +0.0850 | +0.1473 |
| NYC_Q_lemon / NYC-Q_lemon_16 | 111 | +0.0235 | -0.0185 | -0.0316 | +0.0117 |
| NYC_Q_lemon / NYC-Q_lemon_17 | 108 | +0.1578 | +0.1382 | +0.1878 | -0.0004 |
| NYC_Q_lemon / NYC-Q_lemon_18 | 109 | -0.2376 | -0.1978 | -0.2047 | -0.2646 |
| NYC_Q_lemon / NYC-Q_lemon_19 | 111 | -0.0420 | -0.0307 | +0.1185 | +0.1409 |
| NYC_Q_lemon / NYC-Q_lemon_2 | 111 | +0.0333 | -0.0517 | +0.0147 | +0.1068 |
| NYC_Q_lemon / NYC-Q_lemon_20 | 111 | +0.0126 | -0.0072 | +0.2029 | +0.1563 |
| NYC_Q_lemon / NYC-Q_lemon_21 | 110 | -0.0204 | -0.0068 | +0.1100 | +0.0591 |
| NYC_Q_lemon / NYC-Q_lemon_22 | 111 | +0.0791 | +0.0012 | +0.1486 | +0.0379 |
| NYC_Q_lemon / NYC-Q_lemon_23 | 111 | +0.1458 | +0.1216 | +0.0989 | +0.1165 |
| NYC_Q_lemon / NYC-Q_lemon_24 | 110 | +0.1522 | +0.0843 | +0.1928 | +0.2425 |
| NYC_Q_lemon / NYC-Q_lemon_25 | 110 | -0.1960 | +0.0022 | -0.0529 | +0.0149 |
| NYC_Q_lemon / NYC-Q_lemon_26 | 110 | -0.1054 | -0.1458 | -0.2536 | -0.2659 |
| NYC_Q_lemon / NYC-Q_lemon_27 | 109 | +0.1282 | +0.0559 | +0.1017 | +0.0604 |
| NYC_Q_lemon / NYC-Q_lemon_28 | 110 | -0.4609 | -0.3682 | -0.3931 | -0.3683 |
| NYC_Q_lemon / NYC-Q_lemon_29 | 111 | +0.0907 | +0.0097 | +0.0112 | -0.1046 |
| NYC_Q_lemon / NYC-Q_lemon_3 | 110 | +0.1215 | +0.0550 | +0.1266 | +0.0642 |
| NYC_Q_lemon / NYC-Q_lemon_30 | 111 | +0.1408 | +0.0690 | +0.0917 | +0.1127 |
| NYC_Q_lemon / NYC-Q_lemon_31 | 111 | +0.2147 | +0.1363 | +0.2144 | +0.1516 |
| NYC_Q_lemon / NYC-Q_lemon_4 | 111 | -0.0341 | -0.0819 | -0.0886 | -0.2649 |
| NYC_Q_lemon / NYC-Q_lemon_5 | 110 | +0.0042 | -0.0188 | +0.0769 | -0.0919 |
| NYC_Q_lemon / NYC-Q_lemon_6 | 111 | +0.0191 | +0.0587 | +0.1879 | +0.0017 |
| NYC_Q_lemon / NYC-Q_lemon_7 | 110 | +0.1954 | +0.1754 | +0.1297 | +0.0695 |
| NYC_Q_lemon / NYC-Q_lemon_8 | 111 | +0.1251 | +0.1605 | +0.2178 | -0.0143 |
| NYC_Q_lemon / NYC-Q_lemon_9 | 111 | +0.1641 | +0.0742 | -0.0347 | -0.0537 |
| PSQ / PSQ_Demands | 111 | +0.0674 | -0.0466 | +0.0415 | -0.0857 |
| PSQ / PSQ_Joy | 111 | +0.0753 | +0.0195 | -0.0102 | -0.0765 |
| PSQ / PSQ_OverallScore | 111 | +0.1393 | +0.0500 | +0.1024 | +0.0327 |
| PSQ / PSQ_Tension | 111 | +0.1176 | +0.0462 | +0.0904 | +0.1227 |
| PSQ / PSQ_Worries | 111 | +0.1567 | +0.1435 | +0.1639 | +0.0668 |
| STAI_G_X2 / STAI_Trait_Anxiety | 111 | +0.1684 | +0.0884 | +0.1251 | +0.0207 |
| STAXI / STAXI_AC | 111 | -0.0428 | -0.1799 | -0.1582 | -0.0109 |
| STAXI / STAXI_AI | 111 | +0.0297 | +0.0455 | -0.0403 | -0.0590 |
| STAXI / STAXI_AO | 111 | -0.1885 | -0.2582 | +0.0471 | +0.0114 |
| STAXI / STAXI_State_Anger | 111 | -0.0817 | +0.0394 | +0.0254 | +0.0021 |
| STAXI / STAXI_TAR | 111 | +0.1995 | +0.0603 | -0.0305 | +0.0657 |
| STAXI / STAXI_TAT | 111 | +0.0004 | -0.1770 | -0.0832 | -0.2247 |
| STAXI / STAXI_Trait_Anger | 111 | +0.1615 | -0.0145 | -0.0331 | -0.0551 |
| TAS / TAS_Describing | 111 | -0.0136 | -0.1529 | +0.1026 | -0.0089 |
| TAS / TAS_ExternalThinking | 111 | +0.2398 | +0.1673 | +0.1995 | +0.0550 |
| TAS / TAS_Identification | 111 | -0.0422 | -0.0618 | -0.0678 | -0.1773 |
| TAS / TAS_OverallScore | 111 | +0.1595 | +0.0714 | +0.1646 | +0.0780 |
| TEIQue-SF / TeiQueSF_emotionality | 111 | -0.1720 | -0.2375 | +0.0105 | -0.1209 |
| TEIQue-SF / TeiQueSF_self_control | 111 | +0.0925 | +0.0298 | -0.0206 | +0.0953 |
| TEIQue-SF / TeiQueSF_sociability | 111 | -0.0966 | -0.0161 | -0.1192 | -0.1202 |
| TEIQue-SF / TeiQueSF_total | 111 | -0.0355 | -0.0694 | -0.1179 | -0.1768 |
| TEIQue-SF / TeiQueSF_well_being | 111 | +0.1215 | +0.0322 | -0.0221 | -0.1071 |
| TICS / TICS_ChronicWorrying | 111 | +0.0154 | -0.0759 | +0.0457 | +0.0442 |
| TICS / TICS_LackSocialRecognition | 111 | -0.0421 | -0.1885 | -0.1019 | -0.0296 |
| TICS / TICS_PressuereToPerform | 111 | +0.1236 | -0.0124 | +0.1098 | +0.1378 |
| TICS / TICS_ScreeningScale | 111 | +0.0184 | -0.0765 | +0.0375 | -0.0176 |
| TICS / TICS_SocialIsolation | 111 | -0.0211 | -0.0075 | +0.0725 | +0.0170 |
| TICS / TICS_SocialOverload | 111 | -0.0339 | -0.0727 | +0.0750 | +0.0844 |
| TICS / TICS_SocialTension | 111 | +0.1745 | +0.1323 | +0.2385 | +0.0748 |
| TICS / TICS_WorkDemands | 111 | -0.0261 | -0.0929 | -0.0237 | -0.0883 |
| TICS / TICS_WorkDiscontent | 111 | +0.1775 | +0.1435 | +0.2089 | +0.0257 |
| TICS / TICS_WorkOverload | 111 | -0.0619 | -0.1421 | +0.0769 | -0.0051 |
| UPPS / UPPS_lack_perseverance | 111 | +0.3465 | +0.3026 | +0.2624 | +0.1781 |
| UPPS / UPPS_lack_premeditation | 111 | +0.0948 | +0.0239 | +0.0379 | +0.1375 |
| UPPS / UPPS_sens_seek | 111 | +0.2296 | +0.2670 | +0.3043 | +0.1182 |
| UPPS / UPPS_urgency | 111 | +0.1015 | +0.0312 | -0.0189 | +0.1427 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / AUDIT | 110 | -0.1025 | -0.1051 | +0.0221 | -0.1937 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / BSL23_behavior | 79 | +0.0374 | -0.1410 | -0.0549 | -0.2489 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / BSL23_sumscore | 79 | -0.0290 | -0.1352 | +0.0622 | -0.0515 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / DRUG_0=negative_1=Positive | 111 | +0.1320 | +0.1855 | +0.1376 | +0.0279 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / Gender_ 1=female_2=male | 111 | -0.1248 | -0.0001 | +0.0335 | +0.2163 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / Hamilton_Scale | 110 | -0.1439 | -0.0964 | -0.1229 | -0.1417 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / Smoking_num_(Non-smoker=1, Occasional Smoker=2, Smoker=3) | 111 | -0.2165 | -0.0952 | -0.2008 | -0.0719 |
| META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON / Standard_Alcoholunits_Last_28days | 111 | -0.3499 | -0.0539 | +0.0214 | -0.2655 |

## Rebuild this catalog

```sh
python scripts/build_phenotype_catalog.py
```

The builder checks unique label/dimension keys and complete coverage of all 291 fields at K = 3, 6, 10, and 20. It formats the original values without changing their signs or selecting a best K.
