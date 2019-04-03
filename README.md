# Isbi19-tutorial

This repository contains code and slides that are dedicated to reproduce part of the examples presented at ISBI'19 in Venice during the tutorial entitled: "Recent advances in acquisition and reconstruction for Compressed Sensing MRI". The code focuses on basics and recent advances in MR acquisition or design of k-space sampling schemes.  
All aspects related to MRI reconstruction in this tutorial will be teached by [Prof. Jeff Fessler](https://github.com/JeffFessler/MIRT.jl) with code and examples in [Julia](https://julialang.org/) language.

## Slides

Slides area available in pdf format in the [slides](https://github.com/philouc/isbi19-tutorial/tree/master/slides) folder. 

## Code for ISBI'19 Tutorial

### Test data

Some 2D MR synthetic images are available for testing at various image resolution (ie fixed FOV= 256mm and pixel size varying from 1mm down to 125 µm as in plane resolution). Matrix size thus varies from 256x256 up to 2048x2048. Each filename is called **BrainPhantomXXX[X].png** with _XXX_ referring to the image dimension (e.g., 256). A single 2D MRI image collected at 7 Tesla is also available.  

### Matlab

Here, The repository gathers original Matlab code developed by my former PhD students, Nicolas Chauffert and Carole Lazarus and in tight collaboration with [Pierre Weiss](https://www.math.univ-toulouse.fr/~weiss/) and located in the [matlab](https://github.com/philouc/isbi19-tutorial/tree/master/matlab) folder.

### Python

As we're moving to Python langauge, I and my student Nicolas Chartier started to recode the basics for MRI sampling in python. You can find some ipython Notebooks in the [Python](https://github.com/philouc/isbi19-tutorial/tree/master/python) folder.

Note that we also develop our own image reconstruction python package for multiple Fourier imaging modalities, namely [PySAP](https://github.com/CEA-COSMIC/pysap) and its dedicated plugin for MR image: [pysap-mri](https://github.com/CEA-COSMIC/pysap-mri). The main contributors to this plugin are [Loubna El Gueddari](https://www.linkedin.com/in/loubna-el-gueddari-428a64108/) and [Zaccharie Ramzi](https://www.linkedin.com/in/zaccharie-ramzi-043476a5/?originalSubdomain=fr), two PhD candidates under my supervision at [CEA/NeuroSpin](http://joliot.cea.fr/drf/joliot/en/Pages/research_entities/NeuroSpin.aspx).


