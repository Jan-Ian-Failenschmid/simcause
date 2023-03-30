# simcause
This is a collection of functions for conducting causal inferences with the fast causal algorithm (FCI) included in the pcalg package. 
This package is work-in-progress and the functions in it are in a minimally working state. 
I chose to make them public in case anyone finds them helpfull when conducting causal inference with the pcalg package and I hope 
to be able to improve and expand them soon. 

Currently two kinds of functions are featured in this package: 

## fci_precision

This function makes it possible to simulate data from a linear, Gaussian structural causal system, 
varying over different effect and sample sizes. 
Afterwards, it analyses the sensitivity with which FCI identifies the features of the estimated PAG. 
This provides information on the precision with which FCI identifies certain features under the specified 
conditions and can be used similarly to classical power analyses. 

## fci_bootstrap

This function conducts a bootstrap analysis by applying FCI to bootstrap resamples of a given data set. 
FCI only provides the PAG that is most probable given the data. Unfortunately, it does not provide 
any estimates of the uncertainity that is associated with any featrue. This is achieved here, by 
recording in how many bootstraps a feature is identified by FCI. 
 
Multiple plotting functions are included in order to visualize the results from the bootstrap analysis. 
It is possible to plot the proportion in which each feature is identified together with a confidence interval. 
Additionally, a undirected graph can be generated in which each edge is coloured based on the evidence that is 
present for its inclusion. Lastly, a graph can be generated in which features that have been identified in too 
little bootstraps are removed. 

## Installation

The functions can be installed and loaded by running this code 

```r 
devtools::install_github("Jan-Ian-Failenschmid/simcause")
library(simcause)
```
## Licence 

[MIT](https://choosealicense.com/licenses/mit/)
