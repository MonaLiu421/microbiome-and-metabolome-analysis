#!/bin/bash
source activate fastspar
#step1
fastspar --otu_table $INPUT/sparcc-bacteria-input.tsv --correlation sparcc-bacteria_correlation.tsv --covariance sparcc-bacteria_covariance.tsv
#step2
mkdir bootstrap_counts
fastspar_bootstrap --otu_table $INPUT/sparcc-bacteria-input.tsv --number 1000 --prefix bootstrap_counts/fake_data
#step3
mkdir bootstrap_correlation
fastspar --otu_table $INPUT/${NAME}.tsv --correlation ./bootstrap_correlation/${NAME}_cor.tsv --covariance ./bootstrap_correlation/${NAME}_cov.tsv --iterations 50 --threads 16
#step4
fastspar_pvalues --otu_table $INPUT/sparcc-bacteria-input.tsv  --correlation $INPUT/result/sparcc-bacteria_correlation.tsv --prefix $INPUT/result/bootstrap_correlation/fake_data_*_cor --permuta
tions 1000 -t 16 --outfile pvalues.tsv

