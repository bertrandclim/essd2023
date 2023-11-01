# essd2023

Plotting code for figures in paper "A Global Gridded Dataset for Cloud Vertical Structure from Combined CloudSat and CALIPSO Observations" in Earth System Science Data, 2023.

Repository contains plotting code and data files. More information is stored in these files than is shown in figues (e.g. timeseries vs all-time means), so there is some potential for adjusting plots to explore (especially the validation of the dataset) the data in more detail.

The data product is available for download on [zenodo](https://zenodo.org/records/8057791), and code examples for working with the data product are available at [https://github.com/bertrandclim/3S-GEOPROF-COMB](https://github.com/bertrandclim/3S-GEOPROF-COMB).

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
