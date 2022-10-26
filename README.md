# Introduction

Python software for simulating celluar/microbial clonal evolution with Gillespie algorithm usage

# Installation and Usage

The simpliest installation:
```
python -m pip install microbialClonalEvolution
```

Usage in python:
```
from microbialClonalEvolution import mainView
mainView.run()
```

Usage from console:
```
python -m clonalEvolution 
#   for simply run graphic verion of software

python -m clonalEvolution -h
#   to see instructions for command line parametrization
``` 

# Usage Instructions

## General Tab: 
![general tab](https://user-images.githubusercontent.com/110567171/185866346-14d7ae9e-6e87-42ed-8afa-baa764d18db6.jpg) \

file name - used in multiple functions: saving files, plots, parameters; specify name of file (according to windows acceptable characters), data files will have cycle number (name_number.csv/.txt) \
file path - used in multiple functions: saving files, plots, parameters; specify absolute path to folder where simulation data should be saved \
file description - 'TODO' add description in params/data file about simulation specyfied by user \
\
clone diagram - button, fire the muller plot creation. It is needed to specify path with muller plot data generated from simulaiton \
file structure (file-name_muller_plot_'binned'_data.txt): 

```
clone, cells, previous clone
(0, 5000, 0)
```

mutations histograms - button, fire the VAF plot creation for all clones in specified file. File is needed to be one of generated simulation data \
file structure (file-name_cycle.csv/.txt):

```
,Clone number,Cells number,Mean fitness,Mean mutation number,Driver mutation list,Passener mutation list,Previous clone number
0,0,5000,1,0,[],[],0
```

## Generated tab:
![generated tab](https://user-images.githubusercontent.com/110567171/185866348-65f4044f-e12d-46e4-99f8-0de815d5e4ba.jpg) \
select which files should be saved automaticaly after one cycle specified in parameters (skip). Not all types works for both algorithm versions 'TODO'

## Parameters tab:

![parameters](https://user-images.githubusercontent.com/110567171/185866353-a148e844-0aeb-4352-b64e-b1f9062dae39.jpg)

## Plot tab:

![plot](https://user-images.githubusercontent.com/110567171/185866356-aa3da5fe-0942-45dc-b6a9-829d5f40c856.jpg)

## Threads tab:

![threads](https://user-images.githubusercontent.com/110567171/185866359-1eceb813-d664-4f96-abd7-66c702ff2a16.jpg)

## Main controls:

![main controls](https://user-images.githubusercontent.com/110567171/185866351-c42f006e-8023-470e-8210-1adf54d5cac3.jpg)

# License
  
GNU GENERAL PUBLIC LICENSE  Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.
 
 ## Disclaimer
 
 Cellular/Microbial Clonal Evolution simulations basing on Gillespie algorithm.
    Copyright (C) 2022 by Jarosław Gil

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Author

Jarosław Gil, Silesian Univeristy of Technology, Department of Computer Graphics, Vision and Digital Systems.

# Changelog

v1.0.0 - Working single cell algorithm, binned algorithm \
v1.1.0 - fixes \
v1.1.5 - fixes \
v1.2.0 - Added matrix algorithm, process solution (instead of thred, one simulation wokring as new process - much faster) \
