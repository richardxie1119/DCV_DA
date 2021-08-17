# High-throughput single-organelle mass spectrometry measurements

<p align="center">
<img src="https://github.com/richardxie1119/DCV_DA/blob/main/coverart.png" width="800",align="middle">
</p>

This is the code repository to perform data analysis in the paper: 
[Image-Guided MALDI Mass Spectrometry for High-Throughput Single-Organelle Characterization](accepted), D.C. Castro, Y.R. Xie, S.S. Rubakhin, E.V. Romanova, J.V. Sweedler, Nature Methods 2020, accepted

### To reproduce the analysis
- Interactive Jupyter Notebook to reproduce the results for the dense-core vesicle dataset can be found in [this notebook](https://github.com/richardxie1119/DCV_DA/blob/main/reproduce_figures.ipynb)
- **dcv_analysis.py**: pipeline to perform data analysis on the dense-core vesicle dataset
- **vesicle_classification.py**: complete pipeline for vesicle type classification using machine learning

### Data availability
- The MS1 peak data (.xml) that supports the findings of this study are publicly available via the Illinois Data Bank (https://doi.org/10.13012/B2IDB-5949772_V1)
- Data containing multivariate peak intensities processed by MATLAB (.mat) are provided in the [data](data) folder.
