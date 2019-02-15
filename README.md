# MoE
Optimization tools for Mixture of Experts model

In this project we developed a c++ library to ﬁt the Mixture of Experts (ME) model. ME addresses the heterogeneity of functional forms in high dimensional data space, and encourages individual expert to specialize on smaller local problems. To train the ME, we used EM algorithm to ﬁnd the maximum likelihood estimators for the probabilistic models. Since the gating functions are generally nonlinear, the updating equation can not be solved analytically. We developed a library that provides three optimizers to iteratively solve for the MLE in the M-step: gradient descent, a variant of Newton’s method, and a quasi-Newton solver using BroydenFletcher- Goldfarb-Shanno (BFGS) algorithm. We also provided time and convergence analysis comparing these methods and with an existing R package.

