################################################################################
# fResponse
# Functions for simulating and fitting functional responses
################################################################################
# Developed by Peter Hellström, Department of Zoology, Stockholm University
# 2005-2013
# Last updated: November 22nd 2013
################################################################################
# NOTE: It is difficult to get nls to converge for the TypeIIIb function
# for some data sets. By using mle, I have realized that non-convergence for nls
# occurs when the mle-estimate for the exponent theta is less than one. That functional
# response is not sigmoid, but looks like a Type II response. So you better be careful
# when looking at AICc weights - if the Type IIIb response gets the highest weight
# one should also look at the values of the estimated parameter theta. It could be
# that the functional response is best described by a Type IIIb curve, but is not
# sigmoid in shape!
# It also seems that mle finds sensible estimates when nls fails. However,
# mle2 is slower than nls.

# Functions for numerical response analyses
# Updated 2014-03-31

# POPULATION SIZE ESTIMATION BY USING ACCUMULATION CURVES
# Code for R
# By Peter Hellström, Department of Zoology, Stockholm University
# e-mail: peter.hellstrom@zoologi.su.se
# Last updated: September 1 2011

# Suggested updates:
# Sample from different distributions: assign probs. for Poisson & negative binomial. Check package sampling
# Enter recapture probability, and sample according to that.
# Parametric bootstrap sampling?
# Add output of variation for simulations

# IMPORTANT: Simulation part is not working at the moment, code is under re-development.
# Coding, try to remove excessive for-loops in simulation part, and clear up the code.
# Write a specific function that generates data, one that fits the models,
# and a couple of summary and plot-functions.
# How do for-loops compare with sapply and lapply? Not much difference on huge datasets!

# REFERENCES:
# Kohn, M.H., York, E.C., Kamradt, D.A., Haught, G., Sauvajot, R.M., & Wayne, R.K. (1999) Estimating population size by genotyping faeces. Proceedings of the Royal Society B: Biological Sciences, 266, 657-663.
# Eggert, L.S., Eggert, J.A., & Woodruff, D.S. (2003) Estimating population sizes for elusive animals: the forest elephants of Kakum National Park, Ghana. Molecular Ecology, 12, 1389-1402.
# Frantz, A.C., Schaul, M., Pope, L.C., Fack, F., Schley, L., Muller, C.P., & Roper, T.J. (2004) Estimating population size by genotyping remotely plucked hair: the Eurasian badger. Journal of Applied Ecology, 41, 985-995.
# Bellemain, E., Swenson, J.E., Tallmon, D., Brunberg, S., & Taberlet, P. (2005) Estimating population size of elusive animals with DNA from hunter-collected faeces: four methods of brown bears. Conservation Biology, 19, 150-161.
# Petit, E. & Valiere, N. (2006) Estimating population size with noninvasive capture-mark-recapture data. Conservation Biology, 20, 1062-1073.
# Frantz, A.C. & Roper, T.J. (2006) Simulations to assess the performance of different rarefaction methods in estimating population size using small datasets. Conservation Genetics, 7, 315-318.

