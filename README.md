# perifmms

Alex Barnett   1/26/17

Periodized evaluations of Green's function kernels using a black-box fast multipole method (FMM) combined with low-rank periodizing method for essentially arbitrary unit cell geometries.

### Dependencies

1. MATLAB
1. CMCL FMMLIB2 from http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html

### Installation

1. Download using `git`, `svn`, or as a zip (see green button above).
1. Download and install http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html, and make sure the matlab interfaces work.
1. Make sure fmmlib2d/matlab is in your matlab path
1. Testing: in matlab run `cd lap2d; lap2d2p`. Should produce plot of a periodic potential due to some dipoles.

Usage: see driver code at bottom of `lap2d2p`

### References

Some of the ideas are based upon those in the following papers:

A unified integral equation scheme for doubly-periodic Laplace and Stokes boundary value problems in two dimensions, A. H. Barnett, G. Marple, S. Veerapaneni, and L. Zhao, submitted, Comm. Pure Appl. Math., 29 pages (2016). https://arxiv.org/abs/1611.08038

A fast algorithm for simulating multiphase flows through periodic geometries of arbitrary shape, Gary Marple, Alex Barnett, Adrianna Gillman, and Shravan Veerapaneni, SIAM J. Sci. Comput., 38(5), B740-B772 (2016).

Robust fast direct integral equation solver for quasi-periodic scattering problems with a large number of layers, M. H. Cho and A. H. Barnett, Optics Express, 23(2), 1775-1799 (2015)

A new integral representation for quasi-periodic scattering problems in two dimensions, Alex Barnett and Leslie Greengard, BIT Numer. Math. 51, 67-90 (2011)

A new integral representation for quasi-periodic fields and its application to two-dimensional band structure calculations, Alex Barnett and Leslie Greengard, J. Comput. Phys., 229 (19), 6898-6914 (2010). https://arxiv.org/abs/1001.5464

