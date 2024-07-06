# Graph Matching via Convex Relaxations to the simplex
Repository of Matlab code for the paper Graph Matching via Convex Relaxations to the Simplex by Ernesto Araya and Hemant Tyagi

Contains the following folders:

**Algos Graph Matching: 
  - Graph matching algorithms discussed in the paper: proposed Mirror Descent and Projected Gradient descent, Grampa, Umeyama, quadratic programming constrained to the Birkhoff, degree profile.

**Synthetic data:
  - experiments_wigner.m: experiments under the CGW model. Generates fig. 2 (a) and (c) and part of table 2
  - experiments_er.m : experiments under the CER model. Generates fig. 2 (b) and (d) and part of table 2
  - comparison_different_methods.m: comparison between different convex relaxations. Generates fig.1 and table 1.

**Computer vision: 
  - Sim_SHREC16.m: experiments for the computer vision SHREC16 for mathching shapes. Generates fig.4
  - TOPKIDS folder: SHREC16 dataset.

**Facebook: 
  - Sim_facebook.m :  experiments for the facebook dataset. Generates fig.6
  - Standfor3.m: data for facebook experiments.

**Auto sys:
  - auto_sys_dynamics.m : experiments for the autonomous systems dataset. Generates fig.5
  - mat_files: contains the data for the autonomous systems dataset. 
