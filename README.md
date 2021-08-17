# High-throughput single-organelle mass spectrometry measurements

<p align="center">
<img src="https://github.com/richardxie1119/DCV_DA/blob/main/coverart.png" width="800",align="middle">
</p>

This is the code repository to perform data analysis in the paper: 
[Image-Guided MALDI Mass Spectrometry for High-Throughput Single-Organelle Characterization](accepted), D.C. Castro, Y.R. Xie, S.S. Rubakhin, E.V. Romanova, J.V. Sweedler, Nature Methods 2020, accepted

### File descriptions
- [**reproduce_figures.ipynb**](reproduce_figures.ipynb): Interactive Jupyter Notebook to reproduce figures for the dense-core vesicle dataset.
- **dcv_analysis.py**: pipeline to perform data analysis on the dense-core vesicle dataset
- **vesicle_classification.py**: complete pipeline for vesicle type classification using machine learning
- **vesicle_classification.py**: code for single-cell classification through interpretable machine learning published in our previous study ([Xie el al, Anal.Chem. 2020](https://pubs.acs.org/doi/10.1021/acs.analchem.0c01660)])

### Data availability
- The MS1 peak data (.xml) that supports the findings of this study are publicly available via the Illinois Data Bank (https://doi.org/10.13012/B2IDB-5949772_V1)
- Data containing multivariate peak intensities processed by MATLAB (.mat) are provided in the [data](data) folder.
