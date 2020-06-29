# RodiCS

<p align="center">
    <img src="gallery/oscillatory_omeg12alpha065.gif" height="250" alt="Swaying elastic rod"/>
    <img src="gallery/Cy100_wave.gif" height="250" alt="Vortex-induced vibrations with reconfiguration"/>
</p>
<p align="center">
  Elastic rods in action. On the left, a rod sways under an oscillatory flow. Right, a surface wave passes through the rod.
</p>

RodiCS is a Python code that simulates static and dynamic deformation of elastic rods. It solves the Kirchhoff equations with finite elements using the [FEniCS](fenicsproject.org) platform.

This module was implemented during my Masters research project at Polytechnique Montréal:

[Mouad Boudina, Biomechanics of vibrating soft corals, Polytechnique Montréal, 2020]().

Rod simulations are the subject of Chapter 7. One can find the description of Kirchhoff equations, the mechanical and fluid-dynamical load expressions, verification and validation of the code, and finally simulation cases. I also included in it an appendix on the FEniCS implementation in Python describing the syntax and functions used in the code.

## Dependencies

RodiCS works on Python (version 3.8.2), and requires FEniCS (version 2019.1.0), as well as NumPy (version 1.17.4) and Matplotlib (version 3.1.2) modules. Only if data postprocessing or code validation is desired, you'll also need SciPy (version 1.4.1).

## Usage

### Mesh

First you will need a mesh file of the rod. The mesh consists of one-dimensional intervals embedded in a three-dimensional space. It can be saved into a `.xml` file. We give some mesh files in the `xmf_files/` folder. Also, you can generate your own `.xml` file with the desired number of elements with the script `create_xml.py`.

A quick check whether the mesh is well loaded and represented the shape desired, there is the script `read_mesh.py` which reads the mesh and plots it in 3D.

### Execution

The execution file for each simulation type are called:
- `static_exec.py` for static deformation,
- `dynamic_exec.py` for transient simulations,

For each case, the principal functions `run_static` and `run_dynamic` are written respectively in the files `static.py` and `dynamic.py`.

### Arguments

#### Static case
The function `run_static` requires:

- `mesh`: a Mesh object from the FEniCS library,
- `u0`: direction of the load/flow,
- `which_force`: available load case (for now there are the two options `distributed` for a uniform load and `drag` for the fluid-dynamic drag),
- `force_mag`: dimensionless magnitude of the force.

The script `static_test_drag.py` provides a test case of a rod under static drag.

#### Dynamic case
The function `run_dynamic` requires the same argument as in the static case, in addition to:

- `dt`: time step size,
- `Nt`: number of time steps,
- `Ur`: reduced velocity,
- `Gamma`: aspect ratio of the rod,
- `speed_function`: time function of the flow speed magnitude. It takes only the time as an argument. For other parameters, they are accounted as global variables, i.e. defined in the execution script before the function definition.

The script `dynamic_test_oscillatory.py` simulates the motion of a rod under an oscillatory flow.

#### Wake-oscillator model
The function `run_wake_oscillator` in `wake_oscillator.py` is written in the same fashion as `run_dyamic`. Because the wake-oscillator model assumes that the vortex-induced lift is an additional unknown governed by the van der Pol equation (see [Facchinetti et al. (2004)](https://www.sciencedirect.com/science/article/abs/pii/S0889974603001853)), we consider a new mixed space solution with additional subspaces, and the van der Pol equation is included in the final variational form.

The script `wake_oscillator_exec.py` is the execution file where you can input the desired speed profile and material properties of the rod. 

### Verification
For the static case, the script `static_verification.py` simulates the deformation of rods, under a uniformly distributed load, for meshes of different refinement levels. The discretisation error and the observed order of convergence are then calculated and plotted.

For the dynamic case, the script `dynamic_verification.py` also simulates the deformation of rods under a uniformly distributed load in time, with this time fixing the mesh and refining the time step. Here again, the discretisation error and the observed order of convergence are calculated and plotted.

## Contributing
Please feel free to pull requests to improve this module by optimising the solver, proposing a more elaborate rod model, test simulations, visualisation, etc. In the actual code, we only consider uniform load, fluid-dynamic drag, added mass force, and vortex-induced lift. We wish to enrich the code with other forces like gravity/buoyancy, electric/magnetic forces. Another improvement would be to include parameters that vary in space, such as a non-uniform flow speed, heterogeneous rod material, varying rod cross-section, etc. Finally, we want this code to be generalised for branched structures. I also welcome reports of any kind of issues and bugs that are encountered in the code.

Hope you'll enjoy deforming rods!

## License
[MIT](https://choosealicense.com/licenses/mit/)

