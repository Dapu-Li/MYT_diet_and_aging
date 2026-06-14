# MYTH Diet and Aging — Analysis Code

This repository contains the analysis code for the manuscript:

**"A Machine Learning-Derived Dietary Pattern for Aging"**

The study develops and validates the **MYTH (Machine-learning Youthful) Diet**, a 10-component dietary score derived from UK Biobank (N = 191,689) using LightGBM to predict all-cause mortality, and investigates its associations with multi-omics biomarkers and organ-specific aging.

---

## Repository Structure

```
.
├── s1_Diet/                    # Food categorization and food-wide association analysis (R)
│   ├── s0_create_food.R        # Food group construction and categorization
│   ├── s2_food_group.R         # Food group summary statistics
│   ├── s4_food_group_paint.R   # Visualization of food group distributions
│   ├── s5_imp_cox.R            # Cox regression for food-wide association analysis
│   ├── s6_trend_P.R            # P-for-trend analysis across food intake categories
│   ├── s7_RCS.R                # Restricted cubic spline (RCS) dose-response analysis
│   ├── s21_Cox_signature.R     # Cox regression for protein signature analysis
│   ├── s22_med_prosig.R        # Mediation analysis with protein signatures
│   ├── s23_med_Pro_signature.R # Extended protein signature mediation
│   └── s30_table1.R            # Baseline characteristics table
│
├── s2_Machine learning/        # LightGBM model training and evaluation (Python)
│   ├── s1_ImportanceRanking.py # Feature importance ranking via LightGBM
│   ├── s2_SequentalSelection_Loss.py  # Sequential forward selection with custom loss
│   ├── s3_SFS_Plot.py          # Visualization of sequential selection results
│   ├── s4_Prediction.py        # Out-of-sample mortality prediction (10-fold CV)
│   ├── s5_Evaluate.py          # Model evaluation metrics (AUC, Recall, F1)
│   └── Utility/                # Shared utility modules
│       ├── Training_Utilities.py
│       ├── Evaluation_Utilities.py
│       ├── Processing_Utilities.py
│       ├── OneVsRestLightGBMWithCustomizedLoss.py
│       ├── OrdinalClassifier.py
│       ├── FocalLoss.py
│       ├── DelongTest.py
│       └── dl.py
│
├── s3_Disease/                 # Disease association analysis (R)
│   ├── s1_disease_asso.R       # Associations between MYTH score and disease incidence
│   └── s10_disease_aging.R     # MYTH score associations with aging-related diseases
│
├── s4_score/                   # MYTH score calculation and survival analysis (R)
│   ├── s1_MYT.R                # MYTH score construction
│   └── s2_KM.R                 # Kaplan-Meier survival curves
│
├── s5_asso_Organ/              # Organ-specific aging clock associations (R)
│   ├── s1_Organ_cor.R          # Correlation between MYTH score and organ aging
│   └── s20_Organ_asso.R        # Regression of MYTH score on organ aging clocks
│
├── s6_val/                     # External validation (R)
│   └── s11_val.R               # Validation in independent cohort
│
├── s7_med/                     # Multi-omics mediation analysis (R)
│   ├── s11_med_Pro.R           # Proteomic mediation (2,923 proteins)
│   ├── s12_med_Meta.R          # Metabolomic mediation (251 NMR traits)
│   ├── s14_med_Inflammation.R  # Inflammatory marker mediation (9 markers)
│   ├── s15_med_Pro_signature.R # Protein signature mediation
│   └── s40_SEM.R               # Structural equation modeling (diet -> metabolism -> mortality)
│
└── s8_nhanes/                  # NHANES external validation (R)
    ├── s2_nhs_val.R            # NHANES validation of MYTH score
    ├── s4_nhs_MYT.R            # MYTH score calculation in NHANES
    └── s5_nha_val_new.R        # Updated NHANES validation analyses
```

---

## Requirements

### R (>= 4.1.0)
Key packages: `survival`, `rms`, `ggplot2`, `data.table`, `mice`, `lavaan` (SEM)

### Python (>= 3.8)
Key packages: `lightgbm`, `scikit-learn`, `numpy`, `pandas`, `matplotlib`

---

## Data Availability

All analyses were conducted using UK Biobank data and NHANES public-use data. Raw data are not included in this repository due to data access agreements. Researchers may apply for access to UK Biobank data at [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk).

---

## Usage

Each script is self-contained and should be run sequentially within its folder (s1 through s8). Before running, set your working directory and data paths at the top of each file where indicated by `# Set your working directory and data path here`.

---

## Citation

If you use this code, please cite:

> Ma Y, et al. A Machine Learning-Derived Dietary Pattern for Aging. *[Journal]*, 2025.

---

## Contact

Yating Miao — ytmiao@cmu.edu.cn
