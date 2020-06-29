---
title: 'RodiCS: deform Kirchhoff rods with FEniCS'
tags:
  - Python
  - mechanics
  - dynamics
  - slender structures
  - Kirchhoff rods
authors:
  - name: Mouad Boudina
    affiliation: "1"
  - name: Frédérick P. Gosselin
    affiliation: 1
  - name: Stéphane Étienne
    affiliation: 1
affiliations:
 - name: Department of Mechanical Engineering, Polytechnique Montréal, Québec, Canada
   index: 1
date: 42 Month 2020
bibliography: paper.bib

---

# Summary

Elastic rods are long and thin structures, and are the object of many dynamical systems in nature, from beding and swaying plants [@hassani_large_2016, @luhar_wave-induced_2016, @gosselin_mechanics_2019] to swimming flagellated and ciliated bacteria [@brennen_fluid_1977, @lauga_bacterial_2016] and propulsing caringiform fish [@cheng_continuous_1998]. In industry as well, elastic rods serve as models to study drilling risers in offshore sites and tall chimney stacks in constructions that may vibrate in presence of water or air flows [@paidoussis_fluid-structure_2010]. With the prevalence of rods in mechanical problems, it is important to have an efficient numerical tool reproducing faithfully the observed dynamics. RodiCS is a Python code that fulfills this task, and simulates the motion of elastic rods under specified forces.

With the FEniCS platform [@alnaes_fenics_2015], RodiCS solves the Kirchhoff equations of elastic rods [@landau_theory_1986; @audoly_elasticity_2010] using finite elements. The first step in each simulation is to import the mesh of the rod from an `.xml` file. The mesh consists of one-dimensional first order Lagrange elements embedded in a three-dimensional space. Then, external forces applied on the rod need to be listed. In FEniCS, the Unified Form Language (UFL) [@alnaes_unified_2014] ensures a convenient way to write mathematical expressions. Algebraic, vectorial, and differential operators, as well as tensor algebra and index notation, are all included in UFL, so users can enter force expressions are they figure in equations without any difficulty. After applying boundary conditions, the Kirchhoff equations are solved with the Newton method. In the case of transient simulations, these equations are discretised in time using the backward Euler scheme. \autoref{fig:examples}(a) shows static deformation profiles under drag for several load magnitudes, whereas \autoref{fig:examples}(b) records dynamic deformation profiles in an oscillatory flow over a period. Numerical results match well experiments on thin plates, which are also governed by Kirchhoff equations, conducted by @gosselin_drag_2010 and @leclercq_reconfiguration_2018.

![Deformation profiles of experiments on thin plates (black) from @gosselin_drag_2010 and @leclercq_reconfiguration_2018 compared to numerical simulation on elastic rod using RodiCS (purple). (a) Static deformation under drag. (b) Dynamic deformation under an oscillatory flow.\label{fig:examples}](rod_static_dynamic.pdf)

The traditional method to simulate the motion of elastic rods in FEniCS is to solve the elasticity equation in three-dimensions. This method needs a mesh with a large enough number of volumetric elements, hence is costly, time-consuming, and particularly inappropriate for slender structures. A more efficient option is to consider a rod theory that describes deformation based only on the curvilinear coordinate of the rod, which transposes numerically into using one-dimensional elements in three-dimensional space. The only code using one-dimensional elements we found is the @bleyer_numerical_2018 code solving the linear Timoshenko equations. Timoshenko theory is useful when cross-sections experience shear, which is usually the case for thick structures. Thin structures, nevertheless, have more tendency to bend and twist, and seldom experience shear. In RodiCS, the deformation is governed by the Kirchhoff theory, which is suitable for thin rods and accounts for large deflection owing to its nonlinear character.

RodiCS is intended to solve elastic rod simulations subjected to any kind of external forces, should they be constant, time-dependent, space-dependent, or a mix of them. In fact, the original idea of designing such a code came from the shortage in commercial software that include customised forces. Even in the handful of those enabling customised entries, the only way is to write subroutines coded in low-level, user-unfriendly languages. Not only it implies a familiarity with these languages, but also requires users to be well acquainted with the available functions in the software, which is very often impractical to handle after getting lost in lengthy documentations. Coded in a high-level language, RodiCS fades away this difficulty and offers users, in academia and industry alike, a user-friendly and easy-to-implement framework to simulate elastic rod deformation pertaining to a wide range of mechanical problems.

Furthermore, RodiCS can tackle complex mechanical phenomena coupling elastic rod motion with dynamical variables. In fluid-structure interaction, an example of such phenomena is vortex-induced vibrations (VIV), arising due to the periodic vortex shedding downstream from the rod. One of the most successful models of VIV is the wake-oscillator model proposed by @facchinetti_coupling_2004, where the vortex-induced lift has no explicit expression, but is rather a self-excited variable governed by the van der Pol equation with a feedback proportional to the rod acceleration. In RodiCS, the lift is an additional unknown of the problem, and the van der Pol equation with the acceleration feedback is included in the variational formulation. To make simulations closer to reality, RodiCS can combine vortex-induced lift with other forces. For instance, @boudina_corals_2020 included the added mass force to account for inertial effects, as well as the fluid-dynamic drag, so that the rod streamlined with the flow and underwent VIV at the same time, as shown in \autoref{fig:example_drag_VIV}.

![Frontal and lateral deformation profiles of an elastic rod undergoing fluid-dynamic drag along with vortex-induced vibrations. The vertical line in the left refers to the initial configuration of the rod.\label{fig:example_drag_VIV}](rod_drag_VIV_face_profile.pdf)

Finally, we encourage mechanicians and programmers to contribute in developing RodiCS, whether by incorporating more simulation cases with specific forces, optimising computational methods, or even improving visualisation tools, thence bringing valuable benefit in research and education.

# Acknowledgements

The authors acknowledge the financial support from the Simulation-Based Engineering Science (SBES) program, funded through the CREATE program of the National Sciences and Engineering Research Council of Canada (NSERC).

# References
