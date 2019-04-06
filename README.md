# Isbi19-tutorial

This repository contains code and slides that are dedicated to reproduce part of the examples presented at ISBI'19 in Venice during the tutorial entitled: "Recent advances in acquisition and reconstruction for Compressed Sensing MRI". The code focuses on basics and recent advances in MR acquisition or design of k-space sampling schemes.  
All aspects related to MRI reconstruction in this tutorial will be teached by [Prof. Jeff Fessler](https://github.com/JeffFessler/MIRT.jl) with code and examples in [Julia](https://julialang.org/) language.

## Slides

Slides area available in pdf format in the [slides](https://github.com/philouc/isbi19-tutorial/tree/master/slides) folder. 

## Code for ISBI'19 Tutorial

### Test data

Some 2D MR synthetic images are available for testing at various image resolution (ie fixed FOV= 256mm and pixel size varying from 1mm down to 125 µm as in plane resolution). Matrix size thus varies from 256x256 up to 2048x2048. Each filename is called **BrainPhantomXXX[X].png** with _XXX_ referring to the image dimension (e.g., 256). A single 2D MRI image collected at 7 Tesla is also available. We also provide some SPARKLING (_Spreading Projection Algorithm for Rapid K-space sampLING_) trajectories for illustration purposes for two different image sizes (256x256 and 512x512). 

### Matlab

Here, The repository gathers original Matlab code developed by my former PhD student, [Nicolas Chauffert](http://chauffertn.free.fr/) in collaboration with [Pierre Weiss](https://www.math.univ-toulouse.fr/~weiss/) and located in the [matlab](https://github.com/philouc/isbi19-tutorial/tree/master/matlab) folder. It is useful to reproduce the figures based on TSP-sampling (Travelling Salesman Problem), Markov chain sampling (subfolder: [TSP_Markov](https://github.com/philouc/isbi19-tutorial/tree/master/matlab/TSP_Markov)) and illustrate the projection algorithm onto the hardware constraints imposed on the magnetic field gradients (subfolder [gradient_waveform_design](https://github.com/philouc/isbi19-tutorial/tree/master/matlab/gradient_waveform_design)).

### Python

As we're moving to Python language, I and my student Nicolas Chartier started to recode the basics for MRI sampling in python. You can find some ipython Notebooks in the [Python](https://github.com/philouc/isbi19-tutorial/tree/master/python) folder. Note that we illustrate both Cartesian and non-Cartesian sampling, regular and irregular undersampling. Irregular undersampling can be produced using either pseudo-random generation or incoherent optimization-driven sampling like SPARKLING. The code of the latter approach, designed by [Carole Lazarus](https://www.linkedin.com/in/carole-lazarus-b44907a6/?originalSubdomain=fr), [Nicolas Chauffert](http://chauffertn.free.fr/) and [Pierre Weiss](https://www.math.univ-toulouse.fr/~weiss/), is actually not disclosed.

Importantly, we also develop our own image reconstruction python package for multiple _Fourier imaging_ modalities, namely [PySAP](https://github.com/CEA-COSMIC/pysap). These developments are done in collaboration with the [CosmoStat](https://cosmostat.org) team ([J. L. Starck](http://jstarck.cosmostat.org/) in the context of the [COSMIC](https://cosmic.cosmostat.org) project. The two core developers of [PySAP](https://github.com/CEA-COSMIC/pysap) are A. Grigis (antoine.grigis@cea.fr) and [S. Farrens](http://www.cosmostat.org/people/sfarrens). The new organization of [PySAP](https://github.com/CEA-COSMIC/pysap) relies on on separate plugin for each imaging modality, for instance for MRI: [pysap-mri](https://github.com/CEA-COSMIC/pysap-mri). The main contributors to this plugin are [Loubna El Gueddari](https://github.com/LElgueddari) and [Zaccharie Ramzi](https://github.com/zaccharieramzi), two PhD candidates under my supervision at [CEA/NeuroSpin](http://joliot.cea.fr/drf/joliot/en/Pages/research_entities/NeuroSpin.aspx). 

To reproduce the last block of the 7th notebook, you must first install pyNFFT in the following way (on Linux):

* sudo apt install -y libnfft3-dev 
* export CPLUS_INCLUDE_PATH=/usr/include/python3.5; 

Then, if you use Python3.x (x=5,6,7) ans asuming pip3 is installed on your system, run:
* pip3 install cython numpy
* pip3 install git+https://github.com/pyNFFT/pyNFFT.git

Then go [here](https://github.com/CEA-COSMIC/pysap) to check PySAP's package dependencies (including PyQt4).
Run: 
* sudo apt install cmake to compile the C++ dependencies

Last, install the following branch of PySAP (old system with embedded plug-ins):
* pip3 install git+https://github.com/zaccharieramzi/pysap.git@pogm_addition

Note that the brand new release of PySAP has a new organization with separate plug-ins for MRI, astrophyics, tomography, etc.
You will get it if you run:

* pip3 install --user python-pySAP
