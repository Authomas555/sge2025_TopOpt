# Topology optimization of a transformer

[![GitHub license](https://img.shields.io/github/license/Authomas555/sge2025_TopOpt)](https://github.com/Authomas555/sge2025_TopOpt) [![GitHub release](https://img.shields.io/github/release/Authomas555/sge2025_TopOpt.svg)](https://github.com/Authomas555/sge2025_TopOpt/releases/) [![GitHub stars](https://img.shields.io/github/stars/Authomas555/sge2025_TopOpt)](https://github.com/Authomas555/sge2025_TopOpt/stargazers) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15422223.svg)](https://doi.org/10.5281/zenodo.15422223)


## 1) Contents
This repository contains various implementations of topology optimization method to design a magnetic transformer using [NGSolve](https://ngsolve.org/) and [Python](https://www.python.org/). This repository contains several files including [Jupyter notebooks](https://jupyter.org/). The codes associated to our SGE 2025 article are listed in the folder `SGE2025`. We hope to enrich this repo with new algorithms and method in the future.

## 2) Installation 

To run the codes, please install `NGSolve`, please follow the procedure :
1. Create a new Python environnement (for instance using conda) : `conda create -n myEnvName python=3`
2. Install the following packages `pip install jupyter numpy matplotlib ngsolve`
3. Upgrade webgui widgets : `pip install --upgrade webgui_jupyter_widgets`
3. You should now be able to run the content of the folder `SGE2025`.


## 3) License

Copyright (C) 2025 Thomas GAUTHEY (thomas.gauthey@centralesupelec.fr), Théodore CHERRIERE (theodore.cherriere@centralesupelec.fr), Stéphane Gaydier (stephane.gaydier@g2elab.grenoble-inp.fr)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
