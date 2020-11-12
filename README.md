# LeConte data postprocessing

The code in this repository pulls in various data sources from the LeConte
dropbox folder and performs some postprocessing before compiling it into 
netcdf files. 

## How to run the code

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).


[Clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) this repository. 


Using the terminal from within the downloaded repository, install the conda environment: `conda env create -f environment.yml`.


Activate the conda environment: `conda activate lcpp`.


Run the postprocessing scripts located in the `code` folder, for example:

```
cd code
python -W ignore sep2018.py
```

where the `-W ignore` option specifies that warnings should not be displayed. This is not strictly necessary and the code will run without it.

The processed output will be placed, by default, in the `proc` directory. 

<!-- ## How to develop the code

Install the development environment `conda env create -f environment_dev.yml` and look in the `tests` folder. -->