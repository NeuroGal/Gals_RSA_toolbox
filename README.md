# Gals_RSA_toolbox
Matlab code for various analyses which fall under the umbrella of Representational Similarity Analysis (RSA) \ multivariate pattern analysis, or extensions of it. The toolbox was built with temporally-resolved data in mind (specifically intracranial EEG), but I suppose it can be used with other data types too.

For an introduction to the RSA approach see for example: Kriegeskorte and Kievit (2013), https://doi.org/10.1016/j.tics.2013.06.007. The basic idea is to study neural representations by examining all pairwise distances between responses to different stimuli \ conditions, which are grouped together in a Representational Dissimilarity Matrix (RDM).

This was written as part of the analysis for the following manuscript (so please cite if you use it):
Vishne et al., Cell Reports 2023, 'Distinct Ventral Stream and Prefrontal Cortex Representational Dynamics during Sustained Conscious Visual Perception' https://doi.org/10.1016/j.celrep.2023.112752.


The toolbox includes three parts:
1. Code to compute the RDM, including code for noise normalization of the segments prior to the calculation, which has been shown to improve the RDM reliability. For details about the noise normalization procedure see: Guggenmos, Sterzer and Cichy (2018), “Multivariate pattern analysis for MEG: A comparison of dissimilarity measures”, https://doi.org/10.1016/j.neuroimage.2018.02.044. Functions for plotting the RDM in specific time-points and the mean dissimilarity over time (using the entire RDM or parts of it) are also enclosed.
2. Code for computing Exemplar Discriminability Index (EDI) and Category Discriminability Index (CDI), running statistical quantification (see below) and plotting the output. These metrics measure how distinct is the representation for different exemplars \ condition (how unique and consistent the response is relative to other cases). For more details see Nili, Walther, Alink and Kriegeskorte (2020), “Inferring exemplar discriminability in brain representations”,  https://doi.org/10.1371/journal.pone.0232551.
3. Code for calculating single-exemplar information reliability across stimulus repetitions \ stability across time. These are novel correlation-based metrics conceived as part of the study in Vishne et al., 2023 (https://doi.org/10.1101/2022.08.02.502469). The core function is ‘rdm_rel_calc.m’, which computes both Item Reliability and Geometry Reliability used in the paper (see Figure 5A-B), and also enables calculating reliability after removing some of the representational structure (Figure 5C, gray lines in Figure 5D-E). This function calls ‘rdm_corr_stats.m’ for statistical quantification. The other functions in this folder are helper functions used by one of these.


For example usage of the RDM calculation and Reliability parts (1 and 3 above) see https://github.com/NeuroGal/PersistentViewing_paper (folder: ‘4. Exemplars (RSA)’). This repository also includes additional plotting functions for the Reliability metrics (see ‘Helper Functions’: nice_line_plot.m for cross-repetition reliability and plot_results_mat.m for cross-time-point stability).

I also have another repository, dedicated to dimensionality reduction and plotting multivariate patterns which you may want to check out: **state_space_plot** (https://github.com/NeuroGal/state_space_plot). This implements both plotting of each time-point separately, and plotting all together (state-space trajectories), and a lot more, so I encourage you to check it out yourself.


Statistical quantification for the EDI and Reliability metrics (2-3 above) uses functions implemented in my repository **time_resolved_stats** (https://github.com/NeuroGal/time_resolved_stats), which includes correction for multiple comparisons. Briefly, the options are:
- Point-by-point permutations (followed by FDR correction using Benjamini Hochberg (1995) method).
- Max-statistic correction for multiple comparisons. See Nichols & Holmes, 2002; https://doi.org/10.1002/hbm.1058 for more details.
- Cluster-based permutations correction. See Maris & Oostenveld, 2007; https://doi.org/10.1016/j.jneumeth.2007.03.024 for more details.
- If reliability measures are kept as correlations you can also use t-test (as implemented in Matlab for ‘corr’\’partialcorr’), and follow this by FDR correction (similarly to point-by-point permutations).
See ‘rdm_corr_stats.m’ and the statistics repository for more details.


Gal Vishne, June 2023

Twitter: @neuro_gal
