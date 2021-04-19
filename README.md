# LeConte data postprocessing

The code in this repository pulls in various data sources from the LeConte
dropbox folder and performs some postprocessing before compiling it into 
netcdf files. 

**This code will not work without access to the LeConte dropbox!**

The instructions below are for _macOS_, but should mostly work for other unix systems. 

## How to run the code

Install the python package manager [Miniconda](https://docs.conda.io/en/latest/miniconda.html).


[Clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) this repository.
```
git clone https://github.com/jessecusack/LeConte_postprocessing.git
```
The above can also be achieved with ssh keys but requires some [initial setup](https://jdblischak.github.io/2014-09-18-chicago/novice/git/05-sshkeys.html) (recommended, but you need a [github.com](github.com) account). Alternatively, download the repository directly from github.

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

where the `-W ignore` option specifies that warnings should not be displayed. This is not strictly necessary but the text output to terminal may be cleaner. The code makes some assumptions about where the Dropbox folder is located, but different paths can be specified if this fails. 

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
which indicates that you can specify differet paths to the data and save directories, e.g.
```
python sep2018.py -ad some/path/to/adcp_final
```
This might be useful if the defaults fail on your computer.

<!-- ## How to develop the code

Install the development environment `conda env create -f environment_dev.yml` and look in the `tests` folder. -->

## How to install the development environment

The development environment contains more modules for testing. 

```
conda activate lcpp-dev
python -m ipykernel install --user --name lcpp-dev --display-name lcpp-dev
conda deactivate
```

## ADCP processing with R and oce

The [oce](https://dankelley.github.io/oce/) package in R has a lot of tools for processing oceanographic data. To make use of it we have to first install R e.g.

```
brew install r
```

Optionally install the GUI editor 'Rstudio', which is a bit like matlab or spyder (`brew install rstudio`). 

Run R and install the required packages:

```
R
install.packages("oce")
install.packages("ncdf4")
install.packages("RSQLite")  # For loading rsk files.
q()
```

Then the processing scripts can be run from the terminal, e.g.:

```
Rscript ABLE_deep_2018.R
```

### Installing the jupyter R kernel

This is a convenient way to run R from within jupyter lab. I followed the installation instructions [here](https://github.com/IRkernel/IRkernel).

```
R
install.packages('IRkernel')
IRkernel::installspec()
q()
```

## ADCP processing with MATLAB

Dylan Winters at OSU has written some great code for parsing ADCP data. This requires matlab with the navigation and robotics toolboxes installed (and perhaps others, I'm not sure)

RBR [provides a MATLAB toolbox](https://rbr-global.com/support/matlab-tools) for working with their instruments (such as the soloT).

[Rich Pawlowicz provides](https://www.eoas.ubc.ca/~rich/) several MATLAB toolboxes for reading ADCP data from teledyne RDI, as well as CTD data from Seabird among others. 