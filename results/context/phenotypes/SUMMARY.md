# LEMON Identity Yreg Domain Recovery

Question: are shared physiology/metabolic domains more recoverable from broad identity coordinates than personality/affect labels?
Subjects: 111

## Domain Summary

- blood_chemistry, K=3, labels=35: LOOCV_R2=+0.0098, mean_cv_corr=+0.080, mean_abs_cv_corr=+0.144, in_R2=+0.0545
- blood_chemistry, K=6, labels=35: LOOCV_R2=-0.0105, mean_cv_corr=+0.093, mean_abs_cv_corr=+0.170, in_R2=+0.0926
- blood_chemistry, K=10, labels=35: LOOCV_R2=-0.0399, mean_cv_corr=+0.101, mean_abs_cv_corr=+0.152, in_R2=+0.1286
- blood_chemistry, K=20, labels=35: LOOCV_R2=-0.1217, mean_cv_corr=+0.083, mean_abs_cv_corr=+0.130, in_R2=+0.1951
- blood_pressure, K=3, labels=11: LOOCV_R2=+0.1263, mean_cv_corr=+0.339, mean_abs_cv_corr=+0.339, in_R2=+0.1639
- blood_pressure, K=6, labels=11: LOOCV_R2=+0.1008, mean_cv_corr=+0.313, mean_abs_cv_corr=+0.313, in_R2=+0.1849
- blood_pressure, K=10, labels=11: LOOCV_R2=+0.0741, mean_cv_corr=+0.298, mean_abs_cv_corr=+0.298, in_R2=+0.2205
- blood_pressure, K=20, labels=11: LOOCV_R2=+0.0393, mean_cv_corr=+0.302, mean_abs_cv_corr=+0.302, in_R2=+0.3204
- body_size, K=3, labels=4: LOOCV_R2=+0.1192, mean_cv_corr=+0.286, mean_abs_cv_corr=+0.318, in_R2=+0.1673
- body_size, K=6, labels=4: LOOCV_R2=+0.0871, mean_cv_corr=+0.280, mean_abs_cv_corr=+0.280, in_R2=+0.1900
- body_size, K=10, labels=4: LOOCV_R2=+0.0425, mean_cv_corr=+0.245, mean_abs_cv_corr=+0.274, in_R2=+0.2172
- body_size, K=20, labels=4: LOOCV_R2=+0.0080, mean_cv_corr=+0.276, mean_abs_cv_corr=+0.276, in_R2=+0.3401
- cognition, K=3, labels=101: LOOCV_R2=-0.0074, mean_cv_corr=+0.002, mean_abs_cv_corr=+0.152, in_R2=+0.0420
- cognition, K=6, labels=101: LOOCV_R2=-0.0535, mean_cv_corr=+0.007, mean_abs_cv_corr=+0.100, in_R2=+0.0694
- cognition, K=10, labels=101: LOOCV_R2=-0.0796, mean_cv_corr=+0.036, mean_abs_cv_corr=+0.115, in_R2=+0.1227
- cognition, K=20, labels=101: LOOCV_R2=-0.1282, mean_cv_corr=+0.038, mean_abs_cv_corr=+0.142, in_R2=+0.1955
- personality_affect, K=3, labels=140: LOOCV_R2=-0.0089, mean_cv_corr=+0.018, mean_abs_cv_corr=+0.117, in_R2=+0.0381
- personality_affect, K=6, labels=140: LOOCV_R2=-0.0500, mean_cv_corr=+0.002, mean_abs_cv_corr=+0.103, in_R2=+0.0610
- personality_affect, K=10, labels=140: LOOCV_R2=-0.0688, mean_cv_corr=+0.034, mean_abs_cv_corr=+0.104, in_R2=+0.1062
- personality_affect, K=20, labels=140: LOOCV_R2=-0.1394, mean_cv_corr=-0.014, mean_abs_cv_corr=+0.104, in_R2=+0.1410

## Best Labels Per Domain At K=6

### blood_chemistry
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__TSH_in_mU_l: cv_corr=-0.454, in_sample_corr=+0.057, n=108
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__CKDEPI_in_ml_min_1.73m: cv_corr=+0.434, in_sample_corr=+0.512, n=90
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__HBA1C_in_%: cv_corr=+0.345, in_sample_corr=+0.443, n=109
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__HBA1CI_in_mmol_mol: cv_corr=+0.345, in_sample_corr=+0.443, n=109
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__ASAT_in_µkat_l: cv_corr=+0.318, in_sample_corr=+0.438, n=108
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__ALAT_in_ µkat_l: cv_corr=+0.285, in_sample_corr=+0.489, n=108
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__PT_in_%: cv_corr=+0.257, in_sample_corr=+0.373, n=104
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__CRP_in_mg_l: cv_corr=+0.246, in_sample_corr=+0.231, n=87
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__GGT_in_µkat_l: cv_corr=+0.243, in_sample_corr=+0.339, n=108
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__INR: cv_corr=+0.235, in_sample_corr=+0.352, n=104
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__INRVOR: cv_corr=+0.231, in_sample_corr=+0.353, n=104
- Medical_LEMON__Blood_Sample__Blood_Results_LEMON__PLT_in_exp_9 _l: cv_corr=-0.230, in_sample_corr=+0.166, n=107
### blood_pressure
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP2_left_systole: cv_corr=+0.472, in_sample_corr=+0.551, n=109
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP1_right_systole: cv_corr=+0.421, in_sample_corr=+0.498, n=107
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP1_left_diastole: cv_corr=+0.392, in_sample_corr=+0.457, n=110
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP1_right_diastole: cv_corr=+0.382, in_sample_corr=+0.475, n=107
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP2_left_diastole: cv_corr=+0.373, in_sample_corr=+0.459, n=109
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP1_left_systole: cv_corr=+0.364, in_sample_corr=+0.459, n=110
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP3_left_diastole: cv_corr=+0.354, in_sample_corr=+0.439, n=80
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__BP3_left_systole: cv_corr=+0.284, in_sample_corr=+0.417, n=80
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__pulse1_left: cv_corr=+0.196, in_sample_corr=+0.327, n=108
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__pulse1_right: cv_corr=+0.122, in_sample_corr=+0.290, n=106
- Medical_LEMON__Blood_Pressure__Blood_Pressure_LEMON__pulse2_left: cv_corr=+0.083, in_sample_corr=+0.256, n=109
### body_size
- Medical_LEMON__Anthropometry__Anthropometry_LEMON__Waist_cm: cv_corr=+0.494, in_sample_corr=+0.566, n=110
- Medical_LEMON__Anthropometry__Anthropometry_LEMON__Weight_kg: cv_corr=+0.319, in_sample_corr=+0.425, n=110
- Medical_LEMON__Anthropometry__Anthropometry_LEMON__Hip_cm: cv_corr=+0.308, in_sample_corr=+0.451, n=110
- Medical_LEMON__Anthropometry__Anthropometry_LEMON__Height_cm: cv_corr=-0.001, in_sample_corr=+0.237, n=110
### cognition
- Cognitive_Test_Battery_LEMON__LPS__LPS__LPS_1: cv_corr=+0.443, in_sample_corr=+0.522, n=111
- Cognitive_Test_Battery_LEMON__TAP_Incompatibility__TAP-Incompatibility__TAP_I_13: cv_corr=-0.317, in_sample_corr=+0.222, n=111
- Cognitive_Test_Battery_LEMON__TAP_Incompatibility__TAP-Incompatibility__TAP_I_10: cv_corr=-0.306, in_sample_corr=+0.146, n=110
- Cognitive_Test_Battery_LEMON__TAP_Incompatibility__TAP-Incompatibility__TAP_I_12: cv_corr=-0.298, in_sample_corr=+0.119, n=109
- Cognitive_Test_Battery_LEMON__TMT__TMT__TMT_1: cv_corr=+0.276, in_sample_corr=+0.406, n=111
- Cognitive_Test_Battery_LEMON__TMT__TMT__TMT_2: cv_corr=+0.258, in_sample_corr=+0.394, n=111
- Cognitive_Test_Battery_LEMON__RWT__RWT__RWT_15: cv_corr=+0.234, in_sample_corr=+0.346, n=110
- Cognitive_Test_Battery_LEMON__TAP_Incompatibility__TAP-Incompatibility__TAP_I_24: cv_corr=+0.233, in_sample_corr=+0.371, n=110
- Cognitive_Test_Battery_LEMON__RWT__RWT__RWT_16: cv_corr=-0.227, in_sample_corr=+0.112, n=110
- Cognitive_Test_Battery_LEMON__TAP_Incompatibility__TAP-Incompatibility__TAP_I_9: cv_corr=+0.221, in_sample_corr=+0.413, n=111
- Cognitive_Test_Battery_LEMON__TAP_Working_Memory__TAP-Working_Memory__TAP_WM_5: cv_corr=-0.219, in_sample_corr=+0.142, n=109
- Cognitive_Test_Battery_LEMON__CVLT___CVLT__CVLT_13: cv_corr=-0.214, in_sample_corr=+0.145, n=110
### personality_affect
- Emotion_and_Personality_Test_Battery_LEMON__NYC_Q_lemon__NYC-Q_lemon_28: cv_corr=-0.368, in_sample_corr=+0.086, n=110
- Emotion_and_Personality_Test_Battery_LEMON__CERQ__CERQ_PositiveReappraisal: cv_corr=-0.314, in_sample_corr=+0.103, n=111
- Emotion_and_Personality_Test_Battery_LEMON__UPPS__UPPS_lack_perseverance: cv_corr=+0.303, in_sample_corr=+0.403, n=111
- Emotion_and_Personality_Test_Battery_LEMON__LOT-R__LOT_Optimism: cv_corr=+0.302, in_sample_corr=+0.446, n=111
- Emotion_and_Personality_Test_Battery_LEMON__COPE__COPE_Planning: cv_corr=-0.282, in_sample_corr=+0.136, n=111
- Emotion_and_Personality_Test_Battery_LEMON__FTP__FTP_SUM: cv_corr=+0.273, in_sample_corr=+0.392, n=70
- Emotion_and_Personality_Test_Battery_LEMON__MSPSS__MSPSS_SignificantOthers: cv_corr=-0.267, in_sample_corr=+0.124, n=111
- Emotion_and_Personality_Test_Battery_LEMON__UPPS__UPPS_sens_seek: cv_corr=+0.267, in_sample_corr=+0.409, n=111
- Emotion_and_Personality_Test_Battery_LEMON__STAXI__STAXI_AO: cv_corr=-0.258, in_sample_corr=+0.128, n=111
- Emotion_and_Personality_Test_Battery_LEMON__CERQ__CERQ_Rumination: cv_corr=-0.242, in_sample_corr=+0.166, n=111
- Emotion_and_Personality_Test_Battery_LEMON__COPE__COPE_Alkohol‎_Drogen: cv_corr=-0.241, in_sample_corr=+0.235, n=111
- Emotion_and_Personality_Test_Battery_LEMON__MARS__MARS_Behavioral_Distraction: cv_corr=+0.240, in_sample_corr=+0.382, n=103

## Reading Rule

Use domain-level LOOCV and mean label correlations to compare shared recoverability. Individual label rows explain which labels drive each domain, but the domain summary is the main test.
