This computer code is for the implementation of "A joint model for multivariate longitudinal and survival data to discover the conversion to Alzheimer's disease"

Authors: Kai KANG and Xinyuan SONG, The Chinese University of Hong Kong

Files for simulation study:
1. simulation.r
2. generate_data.r
3. mcmc.cpp

Files for real analysis:
1. real_analysis.r
2. mcmc.cpp

Other files:
1. header files, such as f_last.h, f_G.h, f_pre.h and so on. These are specific functions used in MCMC iterations and are incorporated in mcmc.cpp.
2. Txt files, such as DELTA.txt, NT.txt, OT.txt and so on. These are collected ADNI data used in the real analysis.


Simulation study: Use simulation.r
Real analysis: Use real_analysis.r

1. Estimation
1.1 Proposed model: set NU = 6, IF_SEM = 1
1.2 Linear trajectory (M_linear): set NU = 2
1.3 conventional joint model with single indicator (M_ADAS11, M_ADAS13 and M_FAQ): set IF_SEM = 0, NY = 1, and change the dimension of Y
1.4 joint model with independent indicators and trajectories (M_ind): IF_SEM = 0, NY = 3

2. Dynamic prediction
2.1 AUC: set IF_PRE = 1
2.2 Pearson correlation: set IF_LAST = 1
