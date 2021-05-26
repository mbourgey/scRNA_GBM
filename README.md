# scRNA_GBM
Analysis script for GBM single cell RNA seq data

## Data
### Cellranger matrix
The cell_ranger matrices used in the paper analysis could be found [here](https://datahub-262-c54.p.genap.ca/GBM_paper_data/GBM_cellranger_matrix.tar.gz).

### Cell annotation .mat file
The cell annotation matrices used in the paper analysis could be found [here](https://datahub-262-c54.p.genap.ca/GBM_paper_data/annotated_cancer_data.mat).
Cell annotation are provided as a matlab mat file with the corresponding encoding:    

|Cell Type|Encoding|
|---------|--------|
|Unassigned|0|
|Mesenchymal|1|
|Neuronal|2|
|Astro|3|
|Progenitor|4|
|Oligo|5|

### Temporary missing EGA bam
One of the sample bam file (BT400) is temporary missing in the EGA study [EGAS00001004422](https://ega-archive.org/studies/EGAS00001004422) and can be found [here](https://datahub-262-c54.p.genap.ca/GBM_paper_data/BT400.bam) as long with its index file [here](https://datahub-262-c54.p.genap.ca/GBM_paper_data/BT400.bam.bai).
