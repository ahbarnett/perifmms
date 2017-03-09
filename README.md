# perifmms

Alex Barnett   3/8/17

Periodized evaluations of Green's function kernels using a black-box fast multipole method combined with low-rank periodizing scheme for general skew unit cell geometries.

Currently just doubly-periodic Laplace in 2D.

Simplifying assumptions for now:
  * source strengths are compatible w/ periodizing (this is not tested; answers of size 1e16 will result if it does not hold)
  * all sources and targs lie in the unit cell, or close to it.
  * unit cell is centered on origin, and its aspect ratio is not too extreme (up to around 10 is fine).

The answer is only defined up to an overall constant in potential (gradients ans hessians are uniquely defined). The interface is designed to be identical to the CMCL FMM, except with extra arguments describing the unit cell and other options.


### Dependencies

1. MATLAB
1. CMCL FMMLIB2 from http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html

### Installation

1. Download using `git`, `svn`, or as a zip (see green button above).
1. Download and install http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html, and make sure the matlab interfaces work. If you have R2016b you may have trouble with openmp  and mex; if so, poke in the fmmlib2d makefiles to build the single-thread mex executable for now.
1. Make sure fmmlib2d/matlab is in your matlab path.
1. Testing: in matlab run `cd lap2d; lfmm2d2ppart`. Should report errors less than 1e-12 then produce plot of a periodic potential due to some dipoles.

Usage: see test/driver code at bottom of `lfmm2d2ppart`

Note: to run the old 3x3 version, do `addpath old; lap2d2p`

### References

The scheme is a distillation of ideas from the following sequence of papers (in reverse chronological order):

A unified integral equation scheme for doubly-periodic Laplace and Stokes boundary value problems in two dimensions, A. H. Barnett, G. Marple, S. Veerapaneni, and L. Zhao, submitted, Comm. Pure Appl. Math., 29 pages (2016). https://arxiv.org/abs/1611.08038

A fast algorithm for simulating multiphase flows through periodic geometries of arbitrary shape, Gary Marple, Alex Barnett, Adrianna Gillman, and Shravan Veerapaneni, SIAM J. Sci. Comput., 38(5), B740-B772 (2016).

Robust fast direct integral equation solver for quasi-periodic scattering problems with a large number of layers, M. H. Cho and A. H. Barnett, Optics Express, 23(2), 1775-1799 (2015)

A new integral representation for quasi-periodic scattering problems in two dimensions, Alex Barnett and Leslie Greengard, BIT Numer. Math. 51, 67-90 (2011)

A new integral representation for quasi-periodic fields and its application to two-dimensional band structure calculations, Alex Barnett and Leslie Greengard, J. Comput. Phys., 229 (19), 6898-6914 (2010). https://arxiv.org/abs/1001.5464

### Notes

* initial version 1/26/17  
* name change, 2x2 instead of 3x3 neighbor lists, square proxy, 3/8/17  
