# Antibody Point Mutation Design Pipeline (MMCDP)

## Overview

In this study, we developed:
- A **MicroMutate model** and two graph-based models for antibody-antigen interaction prediction:
  - **DeepGCN_Anti**
  - **Sgraph_AB**
- These models were integrated with evolutionary constraints, statistical potentials, molecular dynamics simulations, and metadynamics methods to construct the **Multi-Method Collaborative Design Pipeline (MMCDP)**.

## Folder Structure

### antibody_design_deep_fragment2_larger_n
Contains examples for evolutionary constraints and statistical potentials screening.

### usage_all_redo_no_norm
Contains examples of the **MicroMutate model** and its usage.

### complex_train_redo_usage
Provides usage examples for **DeepGCN_Anti**.

### complex_train_redo_structure_based_reduce_usage
Contains usage examples for **Sgraph_AB**.

## Conda Environments

### DeepGCN_Anti and Sgraph_AB
- **Environmental File**: `RTMscore_n.yml`

### MicroMutate Model
- **Environmental File**: `keras_n.yml`

### Evolutionary Constraints and Statistical Potentials Screening
- Requires **Python 3.8**, **PDBFixer**, and **NCBI-BLAST 2.16.0+**.
