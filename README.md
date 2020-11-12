# LeConte data postprocessing

The code in this repository pulls in various data sources from the LeConte
dropbox folder and performs some postprocessing before compiling it into 
netcdf files. 

## How to run the code

Install the python package manager [Miniconda](https://docs.conda.io/en/latest/miniconda.html).


[Clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) this repository (you need a [github.com](github.com) account).
```
git clone https://github.com/jessecusack/LeConte_postprocessing
```
The above can also be achieved more easily with ssh keys (recommended) but requires some [initial setup](https://jdblischak.github.io/2014-09-18-chicago/novice/git/05-sshkeys.html). 

Using the terminal from within the downloaded repository, install the conda environment.
```
conda env create -f environment.yml
```
This will create a new python environment called `lcpp` with all the packages needed to run the code in this repository. You can see what environments you have with `conda info --env`.

Activate the environment: 
```
conda activate lcpp
```
Note this is not a permanent activation and environments need to be reactivated in every new terminal window.

Run the postprocessing scripts located in the `code` folder, for example:
```
cd code
python -W ignore sep2018.py
```

where the `-W ignore` option specifies that warnings should not be displayed. This is not strictly necessary but the text output to terminal may be cleaner. Note that currently, the code makes some assumptions about where the Dropbox folder is located, but different paths can be specified if this fails. 

The processed output will be placed, by default, in the `proc` directory. 

The processing scripts can accept arguments, for example,
```
python sep2018.py --help
```
will display the optional argument list. You should see something like,
```
usage: sep2018.py [-h] [-ld LECONTE] [-ad ADCP] [-sd SAVE]

optional arguments:
  -h, --help            show this help message and exit
  -ld LECONTE, --leconte LECONTE
                        path to LeConte Dropbox directory
  -ad ADCP, --adcp ADCP
                        path to LeConte ADCP directory
  -sd SAVE, --save SAVE
                        path to save processed data
```
which indicates that you can specify differet paths to the data and save directories. This might be useful if the defaults fail on your computer. 

<!-- ## How to develop the code

Install the development environment `conda env create -f environment_dev.yml` and look in the `tests` folder. -->