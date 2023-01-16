#pragma warning(disable : 4267 4838 4244 4305) // cast warnings

/*
 Matlab usage:
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
 threads : the number of threads to use, set to 0 for auto

 outputs:
 M : DIC matrix (double 3N x 3N)
 b : right hand member (double 3N x 1)
 r : global residual (double 1 x 1)
 R : residual image (single n x m x d)

 note: if you don't ask for R, it will not be computed, saving some time.

 */

#include <limits>
#include <cmath>
#include <mex.h>
#include "dvc.hpp"

#include <thread>
#include <mutex>
#include <vector>

#include <chrono>


enum PolyOrder
{
    poly0,
    poly1,
    poly2
};

// ------------------------------
// the main class to hold all the threaded data and its related member functions
// this includes:
// - image reading functions like interpolation
// - digital image correlation functions
template <typename FP, typename DP, typename SI> // FloatPrecision (FP) DoublePrecision (DP), SignedInteger (SI)
class ThreadData
{
    // std::mutex mutex;
public:
    // each thread will work on their local M and b, which will be summed
    // outside of the threaded workload.
    // Consequently there is no need for locking semantics

    // list of pixels per thread
    SI *list; // list of pixels to process
    SI Nlist; // number of pixels in the list

    // DIC data
    DP *M; // Hessian
    DP *b; // right hand member
    FP *R; // the residual image
    DP *a; // displacement amplitudes

    // image data
    Image3D<FP, DP, SI> f; // the reference image
    Image3D<FP, DP, SI> g; // the deformed image

    // ZOI settings
    DP *C; // ZOI center
    DP *L; // ZOI size

    // shapefun data
    SI Ndof;
    DP grad [3];
    DP *phi, *Lx, *Ly, *Lz;

    // current pixel info
    DP res;        // residual at the current pixel
    SI x, y, z;    // pixel coordinates
    DP ux, uy, uz; // displacement dof

    // RMS
    DP rr;         // accumulating squared residuals
    SI visited_px; // counter for the RMS

    // constructor
    ThreadData()
    {
        rr = 0;
        visited_px = 0;
        Nlist = 0;
    }

    // destructor
    ~ThreadData()
    {
        delete[] list;
        delete[] phi;
        delete[] Lx;
        delete[] Ly;
        delete[] Lz;
        delete[] M;
        delete[] b;
    }

    void initialize(SI N)
    {
        list = new SI[N];

        phi = new DP[Ndof] {1};

        Lx = new DP[Ndof];
        Ly = new DP[Ndof];
        Lz = new DP[Ndof];

        M = new DP[3 * Ndof * 3 * Ndof] {0};
        b = new DP[3 * Ndof * 1] {0};
    }

    void get_Mb() const
    {
        // b
        for (SI i = 0; i < Ndof; i++)
        {
            // multiply image gradient with shapefun
            Lx[i] = grad[0] * phi[i];
            Ly[i] = grad[1] * phi[i];
            Lz[i] = grad[2] * phi[i];

            // right hand member
            b[0 + 3 * i] += Lx[i] * res;
            b[1 + 3 * i] += Ly[i] * res;
            b[2 + 3 * i] += Lz[i] * res;
        }

        // Hessian (upper diagonal)
        for (SI j = 0; j < Ndof; j++)
        {
            for (SI i = 0; i <= j; i++)
            {
                // JN: this will write a few elements in the lower diagonal when i == j,
                // wasting a few CPU cycles
                M[0 + 3 * i + (0 + 3 * j) * 3 * Ndof] += Lx[i] * Lx[j];
                M[0 + 3 * i + (1 + 3 * j) * 3 * Ndof] += Lx[i] * Ly[j];
                M[0 + 3 * i + (2 + 3 * j) * 3 * Ndof] += Lx[i] * Lz[j];

                M[1 + 3 * i + (0 + 3 * j) * 3 * Ndof] += Ly[i] * Lx[j];
                M[1 + 3 * i + (1 + 3 * j) * 3 * Ndof] += Ly[i] * Ly[j];
                M[1 + 3 * i + (2 + 3 * j) * 3 * Ndof] += Ly[i] * Lz[j];

                M[2 + 3 * i + (0 + 3 * j) * 3 * Ndof] += Lz[i] * Lx[j];
                M[2 + 3 * i + (1 + 3 * j) * 3 * Ndof] += Lz[i] * Ly[j];
                M[2 + 3 * i + (2 + 3 * j) * 3 * Ndof] += Lz[i] * Lz[j];
            }
        }
    }

    void get_disp()
    {
        ux = 0;
        uy = 0;
        uz = 0;

        for (SI i = 0; i < Ndof; i++)
        {
            ux += a[0 + 3 * i] * phi[i];
            uy += a[1 + 3 * i] * phi[i];
            uz += a[2 + 3 * i] * phi[i];
        }
    }

    template <PolyOrder P>
    void get_phi()
    {
        DP xn, yn, zn;
        if (P > poly0)
        {
            DP xn = 2 * DP(x - C[0]) / L[0];
            DP yn = 2 * DP(y - C[1]) / L[1];
            DP zn = 2 * DP(z - C[2]) / L[2];
            // DP xn = DP(x - C[0]);
            // DP yn = DP(y - C[1]);
            // DP zn = DP(z - C[2]);

            phi[1] = xn;
            phi[2] = yn;
            phi[3] = zn;
        }
        if (P > poly1)
        {
            phi[4] = xn * xn;
            phi[5] = yn * yn;
            phi[6] = zn * zn;

            phi[7] = yn * zn;
            phi[8] = xn * zn;
            phi[9] = xn * yn;
        }
    }

    // DIC member functions --------------------------------
    template <int want_R, PolyOrder P>
    void dic(SI o)
    {
        // row col indices
        SI i, j, k;
        f.ind2sub(i, j, k, o);

        // current coordinates (in ML def)
        x = j;
        y = i;
        z = k;

        // prepare the shapefunctions
        get_phi<P>();

        // compute the current displacement
        get_disp();

        // deformed coordinates
        DP x2 = x + ux;
        DP y2 = y + uy;
        DP z2 = z + uz;

        // test if out of bounds
        if (f.is_inside(x2, y2, z2) == 0)
        {
            return;
        }

        // gtilde
        DP gt = g.get_val(x2, y2, z2);

        // residual
        res = f.get_val(o) - gt;

        // test for masked pixels
        if (isnan(res))
        {
            return;
        }

        // update the R image
        if (want_R)
        {
            R[o] = FP(res);
        }

        // update the gradient
        f.get_grad(grad, o);

        // Hessian and right hand member
        get_Mb();
        // mexPrintf("o:%12d, [%3d,%3d,%3d], %6.1f, %6.1f, [%6.1f,%6.1f,%6.1f]\n",o,x,y,z,f.get_val(o),res,grad[0],grad[1],grad[2]);

        // update RMS
        rr += res * res;
        ++visited_px;
    }

    // work fun, process a list of pixels
    template <int want_R, PolyOrder P>
    void exec()
    {
        for (SI i = 0; i < Nlist; i++)
        {
            // get the current pixel index
            SI o = list[i];

            // compute the per pixel DIC contribution
            dic<want_R, P>(o);
        }
    }

    // for multithreading -------------------------------------
    template <int want_R, PolyOrder P>
    std::thread spawn()
    {
        return std::thread([=]
                           { exec<want_R, P>(); });
    }
};

// Main function (called by matlab)
// ---------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    // use 64bit integers to allow very large images (typically required in 3D)
    typedef int64_t int_t;

    // error/warning message identifier
    const char messageid[] = "localDVC:main";

    // Processing inputs
    // ------------------------------------
    if ((nrhs < 5) || (nrhs > 6))
    {
        mexErrMsgIdAndTxt(messageid, "incorrect number of inputs");
    }

    // f
    const int f_dims = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *f_siz = mxGetDimensions(prhs[0]);

    // g
    const int g_dims = mxGetNumberOfDimensions(prhs[1]);
    const mwSize *g_siz = mxGetDimensions(prhs[1]);

    // a
    const int Na = mxGetNumberOfElements(prhs[2]);

    // C
    const int NC = mxGetNumberOfElements(prhs[3]);
    if (NC != 3)
    {
        mexErrMsgIdAndTxt(messageid, "C should be have 3 elements");
    }

    // L
    const int NL = mxGetNumberOfElements(prhs[4]);
    if (NL != 3)
    {
        mexErrMsgIdAndTxt(messageid, "L should be have 3 elements");
    }

    mwSize b_siz[2] = {mwSize(2), mwSize(1)};
    mwSize M_siz[2] = {mwSize(2), mwSize(2)};
    int_t Ndof = 0;

    PolyOrder P;
    if (Na == 1 * 3)
    {
        P = poly0;
        Ndof = 1;
        b_siz[0] = 3;
        M_siz[0] = 3;
        M_siz[1] = 3;
    }
    else if (Na == 4 * 3)
    {
        P = poly1;
        Ndof = 4;
        b_siz[0] = 12;
        M_siz[0] = 12;
        M_siz[1] = 12;
    }
    else if (Na == 10 * 3)
    {
        P = poly2;
        Ndof = 10;
        b_siz[0] = 30;
        M_siz[0] = 30;
        M_siz[1] = 30;
    }
    else
    {
        mexErrMsgIdAndTxt(messageid, "unexpected size of <a>, i.e. the initial guess, must be [3,1], [12,1] or [30,1]");
    }

    // threads
    int Nthreads = 0;
    if (nrhs == 6)
    {
        double Nthreads_arg = (double)mxGetScalar(prhs[5]);
        Nthreads = int(Nthreads_arg);
    }
    if (Nthreads == 0)
    {
        Nthreads = std::thread::hardware_concurrency();
    }
    // mexPrintf("Nthreads %d\n",Nthreads);

    // all should be 3D arrays
    if (f_dims != 3)
    {
        mexErrMsgIdAndTxt(messageid, "f should be a 3D matrix");
    }
    if (g_dims != 3)
    {
        mexErrMsgIdAndTxt(messageid, "g should be a 3D matrix");
    }

    if ((f_siz[0] != g_siz[0]) || (f_siz[1] != g_siz[1]) || (f_siz[2] != g_siz[2]))
    {
        mexErrMsgIdAndTxt(messageid, "f and g must be the same size");
    }

    // image size
    int_t n = f_siz[0];
    int_t m = f_siz[1];
    int_t d = f_siz[2];

    // get the pointer to the data
    float *f_img = ((float *)mxGetData(prhs[0]));
    float *g_img = ((float *)mxGetData(prhs[1]));
    double *a = ((double *)mxGetData(prhs[2]));
    double *C = ((double *)mxGetData(prhs[3]));
    double *L = ((double *)mxGetData(prhs[4]));

    // correct for the ML index
    double Cc [3] {C[0]-1,C[1]-1,C[2]-1};

    Image3D<float, double, int_t> f(f_img, n, m, d);
    Image3D<float, double, int_t> g(g_img, n, m, d);

    // M
    plhs[0] = mxCreateNumericArray(2, M_siz, mxDOUBLE_CLASS, mxREAL);
    double *M = ((double *)mxGetData(plhs[0]));

    // b
    plhs[1] = mxCreateNumericArray(2, b_siz, mxDOUBLE_CLASS, mxREAL);
    double *b = ((double *)mxGetData(plhs[1]));

    // r
    mwSize r_siz[2] = {mwSize(1), mwSize(1)};
    plhs[2] = mxCreateNumericArray(2, r_siz, mxDOUBLE_CLASS, mxREAL);
    double *r = ((double *)mxGetData(plhs[2]));

    // R
    int want_R = 0;
    if (nlhs == 4)
    {
        want_R = 1;
    }

    // R
    float *R = nullptr;
    if (want_R)
    {
        mwSize R_siz[3] = {mwSize(n), mwSize(m), mwSize(d)};
        plhs[3] = mxCreateNumericArray(3, R_siz, mxSINGLE_CLASS, mxREAL);
        R = ((float *)mxGetData(plhs[3]));
    }

    // RMS counter
    int_t visited_px = 0;

    // create the threads
    ThreadData<float, double, int_t> *TD = new ThreadData<float, double, int_t>[Nthreads];

    // maximum number of pixels in each list:
    int_t N = int_t(double(L[0] * L[1] * L[2]) / double(Nthreads)) + L[1];
    // mexPrintf("max thread size %d\n",N);

    // fill thread data
    for (int t = 0; t < Nthreads; t++)
    {
        // store the polynomial order
        TD[t].Ndof = Ndof;

        // image data
        TD[t].f = f;
        TD[t].g = g;
        TD[t].R = R;

        // displacement amplitudes
        TD[t].a = a;

        // ZOI centers
        TD[t].C = Cc;
        TD[t].L = L;

        // initialize the internal storage
        TD[t].initialize(N);

    }

    // set entire R image to NaN
    if (want_R)
    {
        for (int_t o = 0; o < n * m * d; o++)
        {
            R[o] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    // distribute the non-NaN pixels of the ZOI
    int t = 0;
    for (int_t k = 0; k < L[2]; k++)
    {
        int_t z = Cc[2] - ((L[2] - 1) / 2) + k;
        for (int_t j = 0; j < L[0]; j++)
        {
            int_t x = Cc[0] - ((L[0] - 1) / 2) + j;
            for (int_t i = 0; i < L[1]; i++)
            {
                int_t y = Cc[1] - ((L[1] - 1) / 2) + i;
                
                if (f.is_inside(x, y, z) == 0)
                {
                    continue;
                }

                int_t o = f.sub2ind(y, x, z);
                if ( std::isnan(f.get_val(o)) )
                {
                    continue;
                }

                // assign the pixel to one of the threads
                TD[t].list[TD[t].Nlist] = o;
                TD[t].Nlist++;
            }
            // cycle the threads
            t = (t + 1) % Nthreads;
        }
    }

    // create a vector of threads
    std::vector<std::thread> v;

    // the multithreaded work
    for (int t = 0; t < Nthreads; t++)
    {
        if (want_R)
        {
            switch (P)
            {
            case poly0:
                v.push_back(std::thread(TD[t].spawn<1,poly0>()));
                break;
            case poly1:
                v.push_back(std::thread(TD[t].spawn<1,poly1>()));
                break;
            case poly2:
                v.push_back(std::thread(TD[t].spawn<1,poly2>()));
                break;
            }            
        }
        else
        {
            switch (P)
            {
            case poly0:
                v.push_back(std::thread(TD[t].spawn<0,poly0>()));
                break;
            case poly1:
                v.push_back(std::thread(TD[t].spawn<0,poly1>()));
                break;
            case poly2:
                v.push_back(std::thread(TD[t].spawn<0,poly2>()));
                break;
            }            
        }
    }

    // combine the threads
    for (int t = 0; t < Nthreads; ++t)
    {
        v.at(t).join();
    }

    // assemble the thread data
    for (int t = 0; t < Nthreads; t++)
    {
        visited_px += TD[t].visited_px;
        r[0] += TD[t].rr;

        // b
        for (int_t i = 0; i < 3 * Ndof; i++)
        {
            // right hand member
            b[i] += TD[t].b[i];
        }

        // symmetric part of M
        for (int_t j = 0; j < 3 * Ndof; j++)
        {
            for (int_t i = 0; i <= j; i++)
            {
                M[i + j * 3 * Ndof] += TD[t].M[i + j * 3 * Ndof];
            }
        }
    }

    // Hessian complete to full matrix
    for (int_t j = 0; j < 3 * Ndof; j++)
    {
        for (int_t i = 0; i < j; i++)
        {
            M[j + i * 3 * Ndof] = M[i + j * 3 * Ndof];
        }
    }

    // finalize the rms
    r[0] = sqrt(r[0] / visited_px);

    // cleanup
    delete[] TD;
}
