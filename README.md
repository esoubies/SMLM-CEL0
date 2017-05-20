# SMLM-CEL0

## Description
Single molecule localization microscopy algorithm based on the minimization of the CEL0 relaxation of the L0 penalized least-squares criterion. The minimization is performed with an iteratively reweighted L1 method.

<p align="center">
<img src="https://github.com/esoubies/SMLM-CEL0/blob/master/imgs/recons1.png"/>
<img src="https://github.com/esoubies/SMLM-CEL0/blob/master/imgs/recons2.png"/>
</p>

## Repository content
* main function **SMLMCEL0.m** 
* function **ComputeNorm_ai.m** which computes the norm of the columns of the used operator
* folder **ToyExample** containing a simple example of use in **script.m** 

## SMLM Challenge 2016
The algorithm has been tested on the 2D high density datasets of the [SMLM challenge 2016](http://bigwww.epfl.ch/smlm/challenge2016/index.html). For these tests, algorithm parameters have been set as follows:
* coefEch: 4  (i.e. each pixel of data images is divided in 4 => pixel size 25nm)
* itmaxIRL: 200  (max number of iterations for outer loop IRL1)
* itmaxFista: 200 (max number iterations for inner loop FISTA)
* PSF variance: 4.5e-3
* lambda: 1.1 (for dataset ER2.N3.HD-2D) and 0.21 (for dataset MT4.N2.HD-2D)

## Reference
[High density molecule localization for super-resolution microscopy using CEL0 based sparse approximation.](https://hal.inria.fr/hal-01443565) Proc. ISBI, 2017. 
Simon Gazagnes, Emmanuel Soubies and Laure Blanc-Féraud.
