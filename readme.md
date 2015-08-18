Adaptive GP

An implementation of a fully nonstationary, heteroscedastic Gaussian process for Matlab. The three key components of a squared exponential kernel -- signal variance, lengthscale and noise variance -- are all modeled as functions with separate GP priors.

MAP and full posterior solutions are supported by gradient descent and by HMC sampling of the posterior.

Currently only supports univariate data.

Simple example (See /demos for more)

```
addpath code
addpath data

load datasets
 
x = Dl.x; y = Dl.y;

gp = nsgp(x, y, 'lso', 'grad');

plotnsgp(gp,true);
```
