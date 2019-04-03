MRI DEMONSTRATOR

N. Chauffert, 2014
(last update Aug. 5 2015)

Matlab codes relative to the published paper:

Variable Density Sampling with Continuous Trajectories
N. Chauffert, P. Ciuciu, J. Kahn, P. Weiss, SIAM Imaging Sciences, 2014:7(4):1962-1992.
Please cite this reference if you use the code.

Provided dependences, shipped with this package:

% concorde solver to compute approximate TSP solutions
http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm
(section Executable Programs)
-> See the concorde subfolder for compilation instructions.

% toolboxes for wavelets (you can use your own wavelet implementation 
by changing the two lines myWT = .... and myIWT = ... in the MRIdemonstrator*D.m).
We used Gabriel Peyré's toolboxes available there:
- https://github.com/gpeyre/numerical-tours/raw/master/matlab/toolbox_general.zip
- https://github.com/gpeyre/numerical-tours/raw/master/matlab/toolbox_signal.zip

MRIdemonstrator2D.m (/3D.m) allows to reproduce the experiments of the paper in 2D (resp. 3D)
In 3D, the code execution takes some (about 3) hours to run.
