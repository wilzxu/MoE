# MoE
Optimization tools for Mixture of Experts model

In this project we developed a c++ library to ﬁt the Mixture of Experts (ME) model. ME addresses the heterogeneity of functional forms in high dimensional data space, and encourages individual expert to specialize on smaller local problems. To train the ME, we used EM algorithm to ﬁnd the maximum likelihood estimators for the probabilistic models. Since the gating functions are generally nonlinear, the updating equation can not be solved analytically. We developed a library that provides three optimizers to iteratively solve for the MLE in the M-step: gradient descent, a variant of Newton’s method, and a quasi-Newton solver using BroydenFletcher- Goldfarb-Shanno (BFGS) algorithm.

## Example
We constructed toy datasets generating from a mixture of linear experts model with 5 components: The dataset consists of 2500 data points generated from the piecewise linear function with overlapping regions.

<p>
<img src="https://github.com/wilzxu/MoE/blob/master/plot/Data_j5.png" width="700">
</p>

By fitting the MoE model with optimizers implemented in the library, we get estimates of the missing variable indicating the true data genereting model in different regions. Convergence of different optimizers are shown below.

<p>
<img src= "https://github.com/wilzxu/MoE/blob/master/plot/convergence.png" width ="700"> 
</p>

# To Use

To compile:

```
g++ -std=c++11 -O -o main -I /path/to/dependencies main.cpp
```

To use hmeEM, please use one of following commands:
```
./main X.txt y.txt 
```
to inputs your own X, y, see X_NO.txt, Y_NO.txt for example) or

```
./main (simulate toy X, y)
```
to show the result of simulation



 



