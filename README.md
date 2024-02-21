# essd2023
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bertrandclim/essd2023/HEAD?labpath=9-9_ESSD-figs-export_R2R_oct26.ipynb)

Plotting code for figures in paper "A Global Gridded Dataset for Cloud Vertical Structure from Combined CloudSat and CALIPSO Observations" in Earth System Science Data, 2023.

Repository contains interactive Jupyter notebook for plotting and data files. More information is stored in these files than is shown in figues (e.g. timeseries vs all-time means), so there is some potential for adjusting plots to explore the data in more detail, especially with the validation figures.

The data product is available for download on [zenodo](https://zenodo.org/records/8057791), and code examples for working with the data product are available [here](https://github.com/bertrandclim/3S-GEOPROF-COMB).

## setup
### on Binder
1. Click the Binder badge to use the notebook remotely in browser. Note it can take a while to spin up if nobody has launched it recently (~20 minutes from testing).

### on local machine
1. Download repository
```bash
git clone https://github.com/bertrandclim/essd2023
```
2. Create environment
```bash
conda env create --prefix essd2023 --file ./essd2023/environment.yml
```
3. Launch environment
```bash
conda activate essd2023
```
4. Install Jupyter Lab
```bash
conda install -c conda-forge jupyterlab
```
5. Launch Jupyter Lab (or Jupyter notebook)
```bash
jupyter-lab
```
6. Open `9-9_ESSD-figs-export_R2R_oct26.ipynb`. Click `run -> run all cells`.
