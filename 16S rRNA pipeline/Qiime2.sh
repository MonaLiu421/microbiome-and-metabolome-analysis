#! /bin/bash
source activate miniconda3/envs/qiime2-2021.11

#step1.import PairedEndFastqManifestPhred33 data
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path $INPUT/manifest --output-path ${tmpdir}/${Work_dir}/demux.qza --input-format PairedEndFastqManifestPhre
d33
#Visualization
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

#step2.merge the sequencing
qiime vsearch join-pairs --i-demultiplexed-seqs $INPUT/result/demux.qza --o-joined-sequences ${tmpdir}/${Work_dir}/demux-joined.qza --p-threads 8
#quality control
qiime quality-filter q-score --i-demux $INPUT/result/demux-joined.qza --o-filtered-sequences ${tmpdir}/${Work_dir}/demux-joined-filtered.qza --o-filter-stats ${tmpdir}/${Work_dir}/demux-
joined-filter-stats.qza
#Visualization
qiime demux summarize --i-data demux-joined-filtered.qza --o-visualization demux-joined-filtered.qzv
qiime metadata tabulate --m-input-file demux-joined-filter-stats.qza --o-visualization demux-joined-filter-stats.qzv
#step3.denoise
qiime deblur denoise-16S --i-demultiplexed-seqs $INPUT/result/demux-danduan.qza --p-trim-length 390 --p-sample-stats --o-representative-sequences ${tmpdir}/${Work_dir}/deblur-danduan-rep
-seqs.qza --o-table ${tmpdir}/${Work_dir}/debulr-danduan-table.qza --o-stats ${tmpdir}/${Work_dir}/deblur-danduan-stats.qza --p-jobs-to-start 32

#step4.database-trainning
qiime feature-classifier extract-reads --i-sequences  $INPUT/silva-138-99-seqs.qza --p-f-primer ACTCCTACGGGAGGCAGCAG --p-r-primer GGACTACHVGGGTWTCTAAT --o-reads ${tmpdir}/${Work_dir}/sil
va-138-99-seqs-338-806.qza  --p-n-jobs 32  
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads $INPUT/silva-138-99-seqs-338-806.qza --i-reference-taxonomy $INPUT/silva-138-99-tax.qza --o-classifier ${tmpdir}/${
Work_dir}/silva-138-99-seqs-338-806_classifier.qza

#step5.taxonomy-sliva-annotation
qiime feature-classifier classify-sklearn --i-classifier $INPUT/database/silva-138-99-seqs-338-806_classifier.qza  --i-reads  $INPUT/result/deblur-joined-filter-rep-seqs-390.qza --o-classification ${tmpdir}/${Work_dir}/taxonomy-deblur-sliva.qza --p-n-jobs 32

#step6.filtered
qiime taxa filter-table --i-table debulr-joined-filter-table-390.qza   --i-taxonomy taxonomy-deblur-sliva.qza --p-exclude mitochondria  --o-filtered-table  debulr-joined-filter-table-390-contam.qza
qiime feature-table summarize --i-table debulr-joined-filter-table-390-contam.qza --o-visualization debulr-joined-filter-table-390-contam.qzv
qiime feature-table filter-samples --i-table  debulr-joined-filter-table-390-contam.qza  --p-min-frequency 4000  --o-filtered-table  debulr-joined-filter-table-390-contam-final.qza
qiime feature-table summarize --i-table  debulr-joined-filter-table-390-contam-final.qza --o-visualization  debulr-joined-filter-table-390-contam-final.qzv

#step7.demultiplexing
#stool
qiime feature-table filter-samples --i-table debulr-joined-filter-table-390-contam-final.qza --m-metadata-file ./sample-metadata-new.tsv --p-where "Type= 'stool'" --o-filtered-table ./stool/stool-filter-table.qza
qiime metadata tabulate --m-input-file stool-filter-table.qza --o-visualization stool-filter-table.qzv
qiime tools export --input-path ./stool-filter-table.qza --output-path ./
biom normalize-table -i feature-table.biom -r -o ./feature-table-norm.biom
biom convert -i ./feature-table-norm.biom -o feature-table-norm.tsv --to-tsv --header-key taxonomy
#vagina
qiime feature-table filter-samples --i-table debulr-joined-filter-table-390-contam-final.qza --m-metadata-file ./sample-metadata-new.tsv  --p-where "[Type]='vagina'" --o-filtered-table ./vagina/vagina-filter-table.qza
qiime metadata tabulate --m-input-file vagina-filter-table.qza --o-visualization vagina-filter-table.qzv
qiime tools export --input-path ./vagina-filter-table.qza --output-path ./
biom normalize-table -i feature-table.biom -r -o ./feature-table-norm.biom
biom convert -i ./feature-table-norm.biom -o feature-table-norm.tsv --to-tsv --header-key taxonomy

#step8.phylogeny
#stool
qiime feature-table filter-seqs --i-data ./deblur-joined-filter-rep-seqs-390.qza --i-table stool-filter-table.qza --o-filtered-data stool-filter-rep-seqs.qza && qiime feature-table tabulate-seqs --i-data  stool-filter-rep-seqs.qza  --o-visualization  stool-filter-rep-seqs.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences stool-filter-rep-seqs.qzv --o-alignment stool_rep_seqs_aligned.qza --o-masked-alignment stool_rep_seqs_masked.qza --p-n-threads 4 --o-tree stool-slivaunrooted-tree.qza --o-rooted-tree stool-slivarooted-tree.qza
#vagina
qiime feature-table filter-seqs --i-data ./deblur-joined-filter-rep-seqs-390.qza --i-table vagina-filter-table.qza --o-filtered-data vagina-filter-rep-seqs.qza && qiime feature-table tabulate-seqs --i-data  vagina-filter-rep-seqs.qza  --o-visualization  vagina-filter-rep-seqs.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences vagina-filter-rep-seqs.qza --o-alignment vagina_rep_seqs_aligned.qza --o-masked-alignment vagina_rep_seqs_masked.qza --p-n-threads 4 --o-tree vagina-slivaunrooted-tree.qza --o-rooted-tree vagina-slivarooted-tree.qza





