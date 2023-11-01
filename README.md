# essd2023

Plotting code for figures in paper "A Global Gridded Dataset for Cloud Vertical Structure from Combined CloudSat and CALIPSO Observations" in Earth System Science Data, 2023.

Repository contains plotting code and data files. More information is stored in these files than is shown in figues (e.g. timeseries vs all-time means), so there is some opportunity to adjust plots to explore (especially the validation of the dataset) in more detail.

## setup
1. Download repository
```
$ git clone https://github.com/bertrandclim/essd2023
```
2. Install dependencies
```
$ conda install -c conda-forge jupyter-lab matplotlib numpy xarray cartopy cmcrameri scipy
```
3. Launch Jupyter Lab (or Jupyter notebook)
```
$ jupyter-lab
```
4. Open `9-9_ESSD-figs-export_R2R_oct26.ipynb`. Click `kernel -> run all cells`.
