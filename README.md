<p align="center">
	<img src="pic_a.PNG" width="450">
</p>

# BioMOR 

This directory contains a code based on FFT and Model Order Reduction for the parametruc analysis of electromagnetic fields and voxel-based human body interaction

-------------------------------------------------------------------

# Description
 
- "FFT_SingleSolution.m"  solves a single problem with a frequency value and nominal resistivity
- "FFT_MORfreq.m"         generates a Reduced Order Model (ROM) with frequency as parameter
- "FFT_MORres.m"          generates a Reduced Order Model (ROM) with resistivity as parameter

All user-settable quantities, e.g. frequency, are contained in the block identified by the 
```BEGIN USER SETTINGS``` / ```END USER SETTINGS``` comments.

Available test cases
--------------------
A simple test case is contained in separate directories under the folder "data". 
Set the "name_dir" variable in "FFT_SingleSolution.m", "FFT_MORfreq.m", "FFT_MORres.m"
to the appropriate directory.

User-defined test cases
-----------------------
Follow the instuctions given in "README.txt" inside the "data" directory.

Results and Visualization
--------------------
In the "res_para" directory, the results generated by 
"FFT_SingleSolution.m", "FFT_MORfreq.m", "FFT_MORres.m" are exported in 
Paraview format.
The Reduced Order Models can be used as shown in the "USE MOR" sections of  "FFT_MORfreq.m" and "FFT_MORres.m"

Credits
--------------------
If you use BioMOR, please consider citing:
 
  [1] R. Torchio et al., "A Fast Tool for the Parametric Analysis of Human Body Exposed to LF Electromagnetic Fields in Biomedical Applications," in  Computer Methods and Programs in Biomedicine, in publishing
 
 [2] [R. Torchio et al., "FFT-PEEC: A Fast Tool From CAD to Power Electronics Simulations," in IEEE Transactions on Power Electronics, doi: 10.1109/TPEL.2021.3092431](https://ieeexplore.ieee.org/document/9465649)
 
 [3] [P. Bettini et al., "FFT-VI: a smart approach for the electromagnetic design of complex systems in large fusion devices," in 2020 Plasma Phys. Control. Fusion](https://doi.org/10.1088/1361-6587/abce8f)
 
and BioMOR itself

 [3] R. Torchio, "BioMOR toolbox", https://github.com/UniPD-DII-ETCOMP/BioMOR  
 
 BioMOR toolbox contains parts or modifications of the following codes:
 
 [1] [VoxHenry (A. C. Yucel et al., "VoxHenry: FFT-Accelerated Inductance Extraction for Voxelized Geometries," in IEEE Transactions on Microwave Theory and Techniques, doi: 10.1109/TMTT.2017.2785842.)](https://github.com/acyucel/VoxHenry) 
 
 [2] [MATAMG - MATLAB Algebraic Multigrid toolbox](https://github.com/parkmh/MATAMG)
 
 [3] [Mesh voxelisation](https://it.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation)
 
 [4] [lgwt](https://it.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes)
 
 [5] [NIfTI_20140122](https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
 
 [6] [find graph conected components](https://it.mathworks.com/matlabcentral/fileexchange/33877-find-graph-conected-components)
 
 [7] [Generate maximally perceptually-distinct colors](https://it.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)

Contacts & Authors
-----------------------
Riccardo Torchio (riccardo.torchio@unipd.it)

