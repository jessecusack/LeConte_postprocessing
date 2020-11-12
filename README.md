# LeConte data postprocessing

The code in this repository pulls in various data sources from the LeConte
dropbox folder and performs some postprocessing before compiling it into 
netcdf files. 

## How to run the code

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Install the conda environment `conda env create -f environment.yml`.
3. Activate the conda environment `conda activate lcpp`.
4. Run the postprocessing scripts located in the `code` folder.

## How to develop the code

Install the development environment `conda env create -f environment_dev.yml`. 