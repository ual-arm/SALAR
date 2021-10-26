# SALAR
SALAR Mechanism Synthesizer

SALAR is an open-source software tool written in Matlab for four-bar mechanism synthesis and distributed under the MIT License. It has been developed by the research groups SAL (Supercomputación - Algoritmos, Spanish for High-performance Computing and Algorithms) and ARM (Automatic Control, Robotics and Mechatronics), from the University of Almería, Spain. The name of SALAR comes from combining those of both groups, and it also happens to be a word in Spanish (salt pan) that refers to some zones in Almería. The software package can be found in: (Not available until acceptance. Attached copy to the reviewers), and it is linked to the paper in [5].

SALAR features:

- A simple but powerful graphical user interface (GUI)
- Problem saving and restoring as plain Matlab Workspaces (.MAT files)
- Two error estimation metrics for optimization:
  - Absolute error: Direct point-to-point comparison between the vertexes of the target path and that achieved by computing the cuadratic error of their coordinates.
  - Normalized Shape-Descriptor Vector (NSDV) error: Comparison between the shape of the target path and that achieved by computing the Euclidean of their NSDV vectors, a new shape encoding methodology included in this tool [5].
- Four state-of-the-art optimization methods avaiable:
  - FMinCon's Interior-Point (FMC): The popular method provided by the Matlab's Optimization Toolbox [1]. Note: This option is not available without a valid license of the referred toolbox.
  - Teaching-Learning-based Optimization (TLBO): A recent and simple-to-configure stochastic and population-based meta-heuristic [2].
  - Differential Evolution (DE): A widely-used stochastic and population-based meta-heuristic with multiple configuration options [3].
  - Málaga University Mechanism Synthesis Algorithm (MUMSA): A recent population-based meta-heuristic specially designed for mechanism synthesis built upon Differential Evolution [4].
- Contest between optimizers: The same problem can be solved by all the optimizers to automatically choose the most promising result.
- Detailed optimization log including average results, standard deviation runtime.
- Cooperation/hybridization between optimizers: The current result can be provided to FMC as its initial point, which might improve the results. The referred current result can come either from a previous optimization contest or from the literature.
- Kinematic analysis: The mechanisms can be simulated to collect the velocities and accelerations of their points of interest as well as those of all their links.

See SALAR in action:

- https://youtu.be/hKKiJfMR18I
- https://youtu.be/ua3B0hjKDhg
- https://youtu.be/uR0xwhLQx40
- https://youtu.be/7wAnM6mSfz4
- https://youtu.be/IdBT-RL8nk4
- https://youtu.be/34YHP0A1EXU

References

* [1] M.A. Branch and A. Grace. MATLAB: Optimization Toolbox: User’s Guide. Math Works, 2020.
* [2] R.V. Rao, V.J. Savsani and D.P. Vakharia. Teaching–learning-based optimization: an optimization method for continuous non-linear large scale problems. Information sciences, 183(1), 1-15, 2012.
* [3] R. Storn and K. Price. Differential evolution–a simple and efficient heuristic for global optimization over continuous spaces. Journal of Global Optimization, 11(4), 341-359, 1997.
* [4] J.A. Cabrera, A. Ortiz, F. Nadal and J.J. Castillo. An evolutionary algorithm for path synthesis of mechanisms. Mechanism and Machine Theory, 46(2), 127-141, 2011.
* [5] J.L. Torres-Moreno, N.C. Cruz, J.D. Álvarez, J.L. Redondo and A. Giménez-Fernandez. An open-source tool for path generation synthesis of four-bar mechanisms. Mechanism and Machine Theroy, Accepted, 2021.
