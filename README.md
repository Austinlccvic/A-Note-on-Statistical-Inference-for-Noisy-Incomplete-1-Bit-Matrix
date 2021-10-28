# A-Note-on-Statistical-Inference-for-Noisy-Incomplete-1-Bit-Matrix


## Description
In this project, we consider the statistical inference for noisy incomplete 1-bit matrix. Instead of observing a subset of real-valued entries of a matrix M, we only have one binary (1-bit) measurement for each entry in this subset, where the binary measurement follows a  Bernoulli distribution whose success probability is determined by the value of the entry. Under a unidimensional nonlinear factor analysis model with the logit link, we obtain a point estimator and derive its asymptotic distribution for any linear form of M and latent factor scores. Moreover, our analysis adopts a flexible missing-entry design that does not require a random sampling scheme as required by most of the existing asymptotic results for matrix completion. The proposed estimator is statistically efficient and optimal, in the sense that the Cramer-Rao lower bound is achieved asymptotically for the model parameters. Two applications are considered, including (1) linking  two forms of an educational test and  (2) linking the roll call voting records from multiple years in the United States senate. The first application enables the comparison between examinees who took different test forms, and the second application allows us to compare the liberal-conservativeness of senators who did not serve in the senate at the same time. For more details, please refer to our paper at https://arxiv.org/abs/2105.01769



## Usage
Functions.R file consists of main functions used for model training, data preprocessing, post-processing and inference. 
Use Case -- Simulation Studies.R contains all the codes for the simulation studies in our paper.
Real Data Application--Senate Voting.R contains R codes for real data applications in ranking the senator's conservative scores using our method.
Real Data Application--University Admission Test Data.R contains codes for real data applications in linking educational testing, i.e. estimation and inference for examinees' ability scores and test items' difficulty scores.

