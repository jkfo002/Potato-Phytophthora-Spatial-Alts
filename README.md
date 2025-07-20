# Potato-Phytophthora-Spatiotemporal-Transcriptome

<img width="1216" height="515" alt="image" src="https://github.com/user-attachments/assets/226bd14a-0d41-4dfb-bd4a-a723c89cb48a" />

We report a comprehensive atlas of spatiotemporal transcriptomes in potato leaves inoculated with Phytophthora infestans at single-cell resolution. By utilizing spatiotemporal atlas, we explore the immune repose at single-cell levels and illustrate the characteristics of cell populations in different immune states and dynamic responses within host microenvironment. This will significantly advance our understanding of the dynamic and heterogeneous nature of plant-pathogen interactions, thereby providing novel strategies for developing disease-resistant crops.

# Environment
```yaml
Singularity sif file based on miniconda3
Avaliable via Google drive: https://drive.google.com/file/d/1WsLoydvRQOpTGL6UIyMCluEcEI4ZLbCc/view?usp=sharing

[ENV name] seuratv43:
  seurat=4.3.0
  monocle3=1.0.0
[ENV name] monocle2:
  monocle=2.22.0
[ENV name] scanpy:
  scipy=1.8.1
  squidpy=1.2.2 
[ENV name] pyscenic:
  pyscenic=0.12.1
  Arboreto=v0.1.6
```

# Content
- Pipeline
  - seurat_pipeline: Basic seurat pipeline for **data processing**, **cell clustering**, **marker gene and DEG analysis**
  - giotto_pipeline: Construct Delaunay network for spatial neighbor identification
- Analysis
  - Spatiotemporal_transcriptomic_atlas: Preprocessing Spatial transcriptome data and supervised Cells clustering
  - Cell-type_expression_patterns: Cell type specifics immune response
  - Pathogen_Targeted_Cells: Characteristics of Pathogen Targeted Cells (PTCs)
  - Host_microenvironment: Dynamic responses within host microenvironment

# Acknowledgements
Thanks to Yuying Li and Zhaonian Dong for developing scripts for these analysis.
