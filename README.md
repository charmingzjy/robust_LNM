# Lesion Network Mapping Null-Model Pipeline

The current scripts are configured for the Amnesia dataset and an 8 mm MNI
space.

## Project Structure

Place the scripts and input files in the following structure:

```text
null_model_pipeline/
├── code_pipeline/
│   ├── S1_true_lesion_r_map_generate.m
│   ├── S2_true_lesion_LNM_generate.m
│   ├── S3_generate_null_model_annealing.m
│   ├── S3_generate_null_model_reposition.m
│   ├── S4_null_model_r_map_generate.m
│   ├── S5_null_model_LNM_generate.m
│   └── S6_permutation_test_FDR.m
├── disease_data/
│   └── Amnesia_mask_voxel_indices_v2.xlsx
├── template/
│   ├── valid_indices_standard_wh_MNI8mm_32_32_32.mat
│   ├── rho_z_mean.mat
│   ├── norm.nii.gz
│   └── norm_8mm.nii.gz
└── output/
```

All paths are constructed relative to the script location. The scripts must
therefore remain inside `null_model_pipeline/code_pipeline/`.

## Requirements

- MATLAB with `readtable`, `readcell`, `writecell`, and local functions in
  scripts
- Statistics and Machine Learning Toolbox for `prctile` and `kstest2`
- FreeSurfer MATLAB utilities, including `MRIread`, on the MATLAB path
- Bioinformatics Toolbox for `mafdr`
- Sufficient disk space for 10,000 null datasets and their lesion r-maps

## Input Data

### Lesion index spreadsheet

`disease_data/Amnesia_mask_voxel_indices_v2.xlsx` must contain:

| Column | Content |
|---|---|
| 1 | Subject or lesion ID |
| 2 | Number of lesion voxels (`NumInd`) |
| 3 onward | Whole-volume linear voxel indices |


## Workflow

### S1: Generate true lesion r-maps

Run:

```matlab
S1_true_lesion_r_map_generate
```

For each real lesion, S1 maps its voxel indices to rows of the mean FC matrix
and averages the voxel-wise FC maps.

Output:

```text
output/R0_true_lesion_r_maps_from_meanFC/
└── Amnesia_rmap_meanFC_<subject_id>.mat
```

Each MAT file contains `lesion_r_map`.

### S2: Generate the true LNM

Run:

```matlab
S2_true_lesion_LNM_generate
```

For each lesion r-map, S2 labels the top 30% as `+1` and the bottom 30% as
`-1`. It then calculates the signed majority overlap across real lesions.

Output:

```text
output/R1_true_lesion_rLNM_map_results/
└── Amnesia_LNM_rmap_30.mat
```

The MAT file contains `overlap`.

### S3: Generate null lesion datasets

Choose one of the following methods. Both preserve each lesion's voxel count,
shape, and orientation.

#### Option A: Simulated annealing

```matlab
S3_generate_null_model_annealing
```

This method:

- randomly repositions lesions inside `norm_8mm.nii.gz`;
- optimizes the sorted lesion-frequency distribution;
- retains candidates whose frequency distribution is not significantly
  different from the true distribution by a KS test.

Output:

```text
output/R2_null_model_random_indices_iter10000_sortedfreq/
├── null_00001.mat
├── null_00001.xlsx
└── ...
```

#### Option B: Random repositioning

```matlab
S3_generate_null_model_reposition
```

This method randomly repositions lesions inside the brain mask without
frequency matching or KS filtering.

Output:

```text
output/R2_null_model_random_indices_iter10000_reposition/
├── null_00001.mat
├── null_00001.xlsx
└── ...
```

### S4: Generate lesion r-maps for each null dataset

Set the method at the top of S4:

```matlab
NULL_MODEL_METHOD = 'annealing';
```

or:

```matlab
NULL_MODEL_METHOD = 'reposition';
```

Then run:

```matlab
S4_null_model_r_map_generate
```

Outputs:

```text
output/R3_null10000_generate_r_map_sortedfreq/
```

or:

```text
output/R3_null10000_generate_r_map_reposition/
```

Set `START_NULL` to resume an interrupted run. With
`SKIP_EXISTING_FILES = true`, existing lesion r-map files are not overwritten.

### S5: Generate one LNM for each null dataset

Set `NULL_MODEL_METHOD` to the same method used in S4, then run:

```matlab
S5_null_model_LNM_generate
```

S5 uses exactly the same LNM calculation as S2:

Outputs:

```text
output/R4_LNM_rmaps_null_10000_results_percent30_sortedfreq/
```

or:

```text
output/R4_LNM_rmaps_null_10000_results_percent30_reposition/
```

Each null MAT file contains `overlap`.

### S6: Permutation test and FDR correction

Set `NULL_MODEL_METHOD` to the same method used in S4 and S5, then run:

```matlab
S6_permutation_test_FDR
```

Outputs:

```text
output/R5_Permutation_FDR_results_percent30_sortedfreq/
```

or:

```text
output/R5_Permutation_FDR_results_percent30_reposition/
```

The final MAT file contains:

- `ratio`: proportion of null values less extreme than the true value
- `pval`: one-sided p-values calculated as `1 - ratio`
- `qval`: Benjamini-Hochberg adjusted p-values
- `significant_mask`: `qval < Q_THRESHOLD`
- `true_overlap`
- `less_extreme_count`
- `num_null`

## Recommended Execution Order

For the annealing route:

```text
S1 -> S2 -> S3 annealing -> S4 annealing -> S5 annealing -> S6 annealing
```

For the reposition route:

```text
S1 -> S2 -> S3 reposition -> S4 reposition -> S5 reposition -> S6 reposition
```

S1 and S2 only need to be run once when both null-model methods use the same
real lesion dataset.
