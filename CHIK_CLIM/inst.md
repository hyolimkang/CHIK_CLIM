# Task: Reporting-Delay Distribution Estimation Pipeline

## Goal
Estimate the reporting-delay distribution for chikungunya cases
using SINAN line-list data, within a Dirichlet-multinomial
Bayesian framework.

## Requirements
- Input: weekly case data at the state level for 11 Brazilian states
- Data directory: readRDS("01_Data/chik_sinan_individual_2015_2024.rds")
- Include a validation method for nowcasting results

## Constraints
- Account for the existing population conservation issue
