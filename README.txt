This is the README for Compositional Maximum Entropy (CME), a general modeling method for compositional data analysis.  The corresponding repository contains 3 main Julia functions (described below) and several example scripts.

The functions require a total of four inputs:

1. Data - a matrix with N (the number of components) rows and D (the number of samples) columns (with each column summing to 1).

2. lambda - a simulation parameter (see the examples).  

3. gamma - an L2 regularization parameter (all examples use gamma=0)

4. t - the simulation length


The four functions include:

1. CME_Fit(Data,gamma) - Computes the CME model parameters: K,h

2. CME(Data,lambda,gamma,t) - Computes the CME model parameters and performs Monte-Carlo validation.  Useful when either N or the entries in K are small.

3. Compare_to_logit(Data) - Computes the CME model parameters as well as the logit-normal model parameters (an alternative method).