---
title: 'RodiCS: a finite element solver of Kirchhoff rods under fluid flow and more'
tags:
  - Python
  - mechanics
  - fluid-structure interaction
  - slender structures
  - Kirchhoff rods
authors:
  - name: Mouad Boudina
    orcid: 0000-0002-4908-4589
    affiliation: "1, 2"
  - name: Stéphane Étienne
    affiliation: 1
  - name: Frédérick P. Gosselin
    orcid: 0000-0003-0639-7419
    affiliation: "1, 2"
affiliations:
 - name: Department of Mechanical Engineering, Polytechnique Montréal, Montréal, Québec, Canada
   index: 1
 - name: Laboratory for Multiscale Mechanics (LM2), Polytechnique Montréal, Montréal, Québec, Canada
   index: 2
date: 42 Month 2020
bibliography: paper.bib

---

# Summary

Elastic rods are long and thin structures, and are the object of many dynamical systems in nature, from bending and swaying plants [@gosselin_mechanics_2019] to swimming flagellated and ciliated bacteria [@brennen_fluid_1977; @lauga_bacterial_2016] and propulsing carangiform fish [@cheng_continuous_1998]. In industry as well, elastic rods serve as models to study drilling risers in offshore sites and tall chimney stacks in constructions that may vibrate in presence of water or air flows [@paidoussis_fluid-structure_2010], and even to investigate the cooking process of spaghetti [@goldberg_mechanics-based_2020]. With the prevalence of rods in mechanical problems, it is important to have a simple to use efficient numerical tool reproducing faithfully the observed dynamics. RodiCS is a Python code that fulfills this task, and simulates the motion of elastic rods under defined forces in general, and their interaction with fluid flows more specifically. While numerous researchers [@violette_computation_2007; @lei_wave_2019; @gosselin_drag_2010; @hassani_large_2016; @leclercq_reconfiguration_2018; @leclercq_vortex-induced_2018; @leclercq_drag_2016; @luhar_wave-induced_2016; @buchak_clapping_2010; @goldberg_mechanics-based_2020] spend valuable time and effort in writing *ad hoc* codes with limited verification and validation, RodiCS comes as a versatile and practical open source library that scientific projects can benefit from and build upon.
 
# Description

With the FEniCS platform [@alnaes_fenics_2015], RodiCS solves the Kirchhoff equations of elastic rods [@landau_theory_1986; @audoly_elasticity_2010] using finite elements. The first step in each simulation is to import the mesh of the rod from an `.xml` file. The mesh consists of one-dimensional first order Lagrange elements embedded in a three-dimensional space. Then, external forces applied on the rod need to be listed. In FEniCS, the Unified Form Language (UFL) [@alnaes_unified_2014] ensures a convenient syntax to write mathematical expressions. Algebraic, vectorial, and differential operators, as well as tensor algebra and index notation, are all included in UFL, so users can enter force expressions exactly as they figure in equations without any difficulty. After applying boundary conditions, Kirchhoff equations are solved with the Newton method. In the case of transient simulations, these equations are discretised in time using the backward Euler scheme. \autoref{fig:examples}(a) shows the static reconfiguration of a 2D beam under a semi-empirical non-linear drag formulation for several load magnitudes, whereas \autoref{fig:examples}(b) records the dynamic deformation of a similar beam under an oscillatory flow over a period. Numerical results match well experiments on thin plates, which are also governed by Kirchhoff equations, conducted by @gosselin_drag_2010 and @leclercq_reconfiguration_2018. A more detailed description of Kirchhoff equations, solving method, and code verification and validation can be found in chapter 7 of @boudina_corals_2020.

![Deformation profiles of experiments on thin plates (black) from @gosselin_drag_2010 and @leclercq_reconfiguration_2018 compared to numerical simulations on elastic rods using RodiCS (purple). (a) Static deformation under drag. (b) Dynamic deformation under an oscillatory flow.\label{fig:examples}](rod_static_dynamic.pdf)

# Statement of Need

The traditional method to simulate the motion of elastic rods in FEniCS is to solve the elasticity equation in three dimensions. This method takes a mesh with a large number of volumetric elements, hence is costly, time-consuming, and particularly inappropriate for slender structures. A more efficient option is to consider a rod theory that describes deformation based only on the curvilinear coordinate of the rod, which transposes numerically into using one-dimensional elements in three-dimensional space. Such a code was proposed by @bleyer_numerical_2018 solving Timoshenko beams. However, it is restricted to static deformations, and is valid only for small displacements due to the linearity of equations. In RodiCS, we simulate both static and dynamic deformations of Kirchhoff rods, which are more representative of thin rods, and we account for large displacements owing to its nonlinear character. Therefore, RodiCS fills a need in open source finite element codes for a way to solve the one-dimensional Kirchhoff rod equation embedded in a three-dimensional space.

Another elastic rod solver worth to mention is the Discrete Elastic Rods (DER) code of @jd2014. Because the DER code is written in C++ and has an intricate file structure, only programmer users are able to modify the source code, excluding thus users new in coding or having basic knowledge in programming. Moreover, it is not cross-platform and has several dependencies, limiting again access to a number of users. Last but not least, @jd2014 acknowledge themselves that they cannot maintain the code or provide support anymore. These drawbacks were incentives for the authors of the present paper to think of writing their own code based on the FEniCS project, an open source finite element solver with an active community, continuously developed, and cross-platform; all significant advantages. Also as an essential modular feature, FEniCS imposes a standard in the mesh definition, making possible either the mesh generation with external meshing software, or the combination of elements in RodiCS with other types of finite elements. Finally, written in Python, RodiCS is inclusive and accessible for users of all programming levels.

RodiCS is intended to solve elastic rod dynamics subjected to any kind of external forces, should they be constant, time-dependent, space-dependent, or a mix of them. In fact, the original idea of designing such a code came from the shortage in commercial software considering customised forces. Even in the handful of those enabling customised entries, the only way is to write subroutines coded in low-level, user-unfriendly languages. Not only does it imply a familiarity with these languages, but it also requires users to be well acquainted with the available user-defined function environment in the software, which is rarely documented in a practical manner. Coded in a high-level language, RodiCS does away with this difficulty and offers users, in academia and industry alike, a user-friendly and easy-to-implement framework to simulate elastic rod deformation pertaining to a wide range of mechanical problems.

# Complex Dynamics

Furthermore, RodiCS can tackle complex mechanical phenomena coupling elastic rod motion with dynamical variables. In fluid-structure interaction, an example of such phenomena is vortex-induced vibrations (VIV), arising due to the periodic vortex shedding downstream from the rod. One of the most successful phenomenological models of VIV is the wake-oscillator model proposed by @facchinetti_coupling_2004, where the vortex-induced lift has no explicit expression, but is rather a self-excited variable governed by the van der Pol equation with a feedback proportional to the rod acceleration. In RodiCS, the lift is an additional unknown of the problem, and the van der Pol equation with the acceleration feedback is included in the variational formulation. RodiCS can even combine vortex-induced lift with other forces and make simulations closer to reality. To reproduce the motion of a vibrating soft coral branch, @boudina_jfm included the added mass force to account for inertial effects, as well as the semi-empirical non-linear fluid-dynamic drag so that the rod bends with large amplitude and undergoes VIV at the same time, as shown in \autoref{fig:example_drag_VIV}.

![An elastic rod undergoing fluid-dynamic drag along with vortex-induced vibrations. Views from left to right: perspective, lateral, and frontal. Arrows indicate the direction of the flow.\label{fig:example_drag_VIV}](rod_VIV.pdf)

Finally, we encourage mechanicians and programmers to contribute in developing RodiCS, whether by incorporating more simulation cases with specific forces, optimising computational methods, or even improving visualisation tools, thence bringing valuable benefit in research and education.

# Acknowledgements

The authors acknowledge the financial support from the Simulation-Based Engineering Science (SBES) program, funded through the CREATE grant of the National Sciences and Engineering Research Council of Canada (NSERC), as well as from Discovery Grants Nos. RGPIN-2019-05335 and RGPIN-2019-07072.

# References
