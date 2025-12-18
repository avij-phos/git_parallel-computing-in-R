
# Parallel Computing Benchmarking for Microbiome Analysis

This repository contains R code benchmarking 4 parallel computing methods:
1. Sequential (baseline)
2. CAR package
3. Future/Furrr framework
4. Mirai package

Results show Mirai achieves 5× speedup over sequential processing.
To run this Analysis, place  a file named 'feature_table.txt' in this directory.
The file should be a tab-separated feature table with atleast 3 sample column.
see the R script for the exact format requirements.
