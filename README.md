Spot-On (cli)
-------

This repository collects a series of scripts to analyze data from single particle tracking experiments. The code was initially written by Anders Sejr Hansen and translated to Python by Maxime Woringer. A [Matlab version](https://gitlab.com/anders.sejr.hansen/spot-on-matlab) exists and is maintained by Anders Sejr Hansen.

This repository only includes the commandline analysis pipeline. A graphical user interface (GUI) is available for a more user-friendly analysis, but is not included in this repository. This repository contains the command-line version, that can be used independently from the GUI.

Although the functions and methods can be called directly, we provide a walk-through tutorial as a [Jupyter](http://jupyter.org) notebook. 

# Dependencies

- numpy
- scipy
- lmfit

Optional: jupyter

Your package manager may provide precompiled versions of `numpy` and `scipy`. In that case, it might be worth using those libraries, because compilation can take a significant amount of time.

# Installation

## Install dependencies
`pip install -r requirements.txt`

Alternatively, you can install the dependencies manually by typing:
`pip install numpy scipy lmfit`

Optional : `pip install jupyter`

## Install fastSPT

Simply run (as root): `python setup.py install`

Check that it worked: `python -c "import fastspt"`


# Tutorial
A short tutorial, demonstrating the capabilities of the software, is available as a Jupyter notebook. 

See `fastSPT_tutorial.ipynb`. The tutorial is also available online: **address of the page**

**Document here how to open a Jupyter notebook**

# Usage
## Main functions
## Input file format
## Caveats

# References

# License
This program is released under the GNU General Public License version 3 or upper (GPLv3+).


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.



# Authors

# Bugs/suggestions
Send to bugtracker or to email.
