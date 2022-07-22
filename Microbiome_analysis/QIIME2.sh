#Create a new environment in conda and install qiime2
conda env create -n qiime2 --file qiime2.yml
#Activate the environment
conda activate qiime2
#Set the currecnt directory
cd ~/Microbiome_analysis
#Make a new directory for qza file and visualization file produced by qiime2
mkdir ./Normalized_data/Visualization
mkdir ./Normalized_data/qza



#We will make PCoA plot with the normalized data with 36 samples
#Convert text file (Normalized_new_table.txt) to biom file
biom convert \
-i Normalized_data/Normalized_new_table.txt \
-o Normalized_data/qza/Normalized_new_table.biom \
--to-hdf5 \
--table-type="Taxon table" \
--process-obs-metadata taxonomy

#Import biom into QIIME2
qiime tools import \
--input-path Normalized_data/qza/Normalized_new_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path Normalized_data/qza/Normalized_new_FeatureTable_frequency.qza

#Crease a distant matrix based on Bray Curtis
qiime diversity beta \
--i-table Normalized_data/qza/Normalized_new_FeatureTable_frequency.qza \
--p-metric braycurtis \
--o-distance-matrix Normalized_data/qza/Normalized_new_distance_matrix.qza

#Apply principal coordinate analysis
qiime diversity pcoa \
--i-distance-matrix Normalized_data/qza/Normalized_new_distance_matrix.qza  \
--o-pcoa Normalized_data/qza/Normalized_new_pcoa.qza

#Visualize PCoA plot
qiime emperor plot \
--i-pcoa Normalized_data/qza/Normalized_new_pcoa.qza \
--m-metadata-file metadata.txt \
--o-visualization Normalized_data/Visualization/Normalized_new_pcoa.qzv

#View the PCoA plot
qiime tools view Normalized_data/Visualization/Normalized_new_pcoa.qzv



#We will check the number of observed features and sequencing depth of each sample
#Compute the number of observed features for each sample in a feature table
qiime diversity-lib observed-features \
--i-table OTU_table/qza/nonfiltered_FeatureTable_frequency.qza \
--o-vector OTU_table/qza/nonfiltered_SampleData.qza   

#Export the result as tsv file
qiime tools export --input-path OTU_table/qza/nonfiltered_SampleData.qza  --output-path Normalized_data/data.txt


#Check the number of sequencing 
qiime feature-table summarize \
--i-table OTU_table/qza/nonfiltered_FeatureTable_frequency.qza \
--m-sample-metadata-file metadata.txt \
--o-visualization OTU_table/qza/nonfiltered-summarize.qzv \
--verbose

#Visualise data
qiime tools view OTU_table/qza/nonfiltered-summarize.qzv





#After extraxting two samples, we performed following analysis using QIIME2

#Convert text file (Normalized_new_table.txt) to biom file
biom convert \
-i Normalized_data/0.1filtered_perm_Normalized_table.txt \
-o Normalized_data/qza/0.1filtered_table.biom \
--to-hdf5 \
--table-type="Taxon table" \
--process-obs-metadata taxonomy

#Import biom into QIIME2
qiime tools import \
--input-path Normalized_data/qza/0.1filtered_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path Normalized_data/qza/0.1filtered_FeatureTable_frequency.qza


#Crease a distant matrix based on Bray Curtis
qiime diversity beta \
--i-table Normalized_data/qza/0.1filtered_FeatureTable_frequency.qza \
--p-metric braycurtis \
--o-distance-matrix Normalized_data/qza/0.1filtered_distance_matrix.qza

#Apply principal coordinate analysis
qiime diversity pcoa \
--i-distance-matrix Normalized_data/qza/0.1filtered_distance_matrix.qza  \
--o-pcoa Normalized_data/qza/0.1filtered_pcoa.qza

#Visualize PCoA plot
qiime emperor plot \
--i-pcoa Normalized_data/qza/0.1filtered_pcoa.qza \
--m-metadata-file metadata_34.txt \
--o-visualization Normalized_data/Visualization/0.1filtered_pcoa.qzv

#View the PCoA plot
qiime tools view Normalized_data/Visualization/0.1filtered_pcoa.qzv


#Calculate alpha diversity
#Observed features
qiime diversity alpha \
--i-table Normalized_data/qza/0.1filtered_FeatureTable_frequency.qza \
--p-metric observed_features \
--o-alpha-diversity Normalized_data/qza/0.1filtered_observed_features.qza

#Visualiza observed feature 
qiime diversity alpha-group-significance \
--i-alpha-diversity Normalized_data/qza/0.1filtered_observed_features.qza \
--m-metadata-file metadata_34.txt  \
--o-visualization Normalized_data/Visualization/0.1filtered_observed_features.qzv

#View observed feature
qiime tools view  Normalized_data/Visualization/0.1filtered_observed_features.qzv


#chao1
qiime diversity alpha \
--i-table Normalized_data/qza/0.1filtered_FeatureTable_frequency.qza  \
--p-metric chao1 \
--o-alpha-diversity Normalized_data/qza/0.1filtered_chao1.qza

#Visualize chao1
qiime diversity alpha-group-significance \
--i-alpha-diversity Normalized_data/qza/0.1filtered_chao1.qza \
--m-metadata-file metadata_34.txt \
--o-visualization Normalized_data/Visualization/0.1filtered_chao1.qzv

#View chao1
qiime tools view Normalized_data/Visualization/0.1filtered_chao1.qzv


#shannon index
qiime diversity alpha \
--i-table Normalized_data/qza/0.1filtered_FeatureTable_frequency.qza \
--p-metric shannon \
--o-alpha-diversity Normalized_data/qza/0.1filtered_shannon.qza

#Visualize shannon index
qiime diversity alpha-group-significance \
--i-alpha-diversity Normalized_data/qza/0.1filtered_shannon.qza \
--m-metadata-file metadata_34.txt \
--o-visualization Normalized_data/Visualization/0.1filtered_shannon.qzv

#View shannon index
qiime tools view Normalized_data/Visualization/0.1filtered_shannon.qzv


#ANCOM analysis
#add one count ton every value (pseudocount)
qiime composition add-pseudocount \
--i-table Normalized_data/qza/0.1filtered_FeatureTable_frequency.qza  \
--o-composition-table Normalized_data/qza/0.1filtered_comp-table.qza

#ANCOM generated volcano plots showing differentially abundant features between groups
qiime composition ancom \
--i-table Normalized_data/qza/0.1filtered_comp-table.qza \
--m-metadata-file metadata_34.txt \
--m-metadata-column Pre/Post \
--o-visualization Normalized_data/Visualization/0.1filtered_ancom.qzv

#View ANCOM result
qiime tools view Normalized_data/Visualization/0.1filtered_ancom.qzv

