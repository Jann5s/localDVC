# localDVC
A very simple local Digital Volume Correlation (DVC) Matlab MEX implementation in C++


This is a Matlab MEX implementation of a local DVC algorithm. Meaning, that a function is programmed in C++ which, when compiled using Matlab, can be called from Matlab. The advantage of such a system is that the C++ part can be quite optimized using lower level sematics available in C++ while using the function in prototype code from the Matlab side.

The implementation is very simple. If you are looking for a DIC or DVC tool, this is not for you. Consider this as an accademic implementation that can be compared to other DVC implementations. The idea behind the implementation is that, it will return the Hessian (M) and Jacobian (b) matrices as well as the residual for a given subset size and initial guess (a). Solving for the update in the unknowns (e.g. a = a + da) should be done on the Matlab side as da = M \ b. This implementation understands matching functions of order 0 (constant) to 2 (quadratic). Perhaps unexpectedly, it requires all degrees of freedom (for all matching functions) to be given as initial guess. It uses a basic tri-cubic interpolator to obtain the gray values at non-integer locations.

See the localDVC_test.m file for some examples of how to compile and how to use this function.

# The compiled mex function has the following syntax
[M, b, r, R] = localDVC_kernel(f, g, a, C, L, threads)

inputs:
 f    : reference image, (n x m x d, single, use NaN to mask)
 g    : deformed image, (n x m x d, single, use NaN to mask)
 a    : initial guess, (3N x 1)
        if N == 1,  0th order matching functions
        if N == 4,  1st order matching functions
        if N == 10, 2nd order matching functions

        phi0 = [x^0 y^0 z^0]
        phi1 = [x^1 y^0 z^0]
        phi2 = [x^0 y^1 z^0]
        phi3 = [x^0 y^0 z^1]
        phi4 = [x^2 y^0 z^0]
        phi5 = [x^0 y^2 z^0]
        phi6 = [x^0 y^0 z^2]
        phi7 = [x^0 y^1 z^1]
        phi8 = [x^1 y^0 z^1]
        phi9 = [x^1 y^1 z^0]

 C    : the center of the zoi: [Cx, Cy, Cz], double
 L    : the size of the zoi (must be odd): [Lx, Ly, Lz], double
 threads : the number of threads to use, set to 0 for auto, double

 outputs:
 M : DIC matrix (double 3N x 3N)
 b : right hand member (double 3N x 1)
 r : global residual (double 1 x 1)
 R : residual image (single n x m x d)

