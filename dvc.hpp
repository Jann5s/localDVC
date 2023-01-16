#pragma once

#include <limits>
#include <cmath>


template< typename FP, typename DP, typename SI> 
// FloatPrecision (FP) DoublePrecision SignedInteger (SI)
// The idea is that:
// FP: input image data
// DP: coordinate data and interpolated image data
// SI: indices
class Image3D
	// The member functions inside are designed to be fast, not necessarily safe. 
	// i.e. you can ask for pixels outside of the image
{
public:
	SI n, m, d; // image size: n-rows, m-cols, d-slices: n * m * d
	FP *img;  // pointer to the image data
	SI extent = 3; // edge margin

	// constructor
	Image3D(){
	}

	Image3D(FP *img, SI n, SI m, SI d) : img(img), n(n), m(m), d(d) {
	}
	
	// image member functions ----------------------------------
	void ind2sub(SI &i, SI &j, SI &k, const SI &o) const {
		// convert linear index (pointer) to ijk indices
		i = o % n;
		j = (o / n) % m;
		k = (o / (n * m)) % d;
	}

	SI sub2ind(const SI &i, const SI &j, const SI &k) const {
		// convert ijk indices to linear index (pointer)
		return i + n * j + n * m * k;
	}

	SI is_inside(const DP &x, const DP &y, const DP &z) const {
		// test of a coordinate is inside the image
		SI val = 1;
		if ((x < extent) || (x > m - 1 - extent) || (y < extent) || (y > n - 1 - extent) || (z < extent) || (z > d - 1 - extent)) {
			val = 0;
		}
		return val;
	}

	SI is_inside(const SI &x, const SI &y, const SI &z) const {
		// test of a coordinate is inside the image
		SI val = 1;
		if ((x < extent) || (x > m - 1 - extent) || (y < extent) || (y > n - 1 - extent) || (z < extent) || (z > d - 1 - extent)) {
			val = 0;
		}
		return val;
	}

	FP get_val(const SI &o) const {
		// get grayvalue
		return img[o];
	}

	FP get_val(const SI &x, const SI &y, const SI &z) const {
		// get grayvalue
		SI o = sub2ind(y, x, z);
		return img[o];
	}

	DP get_val(const DP &x, const DP &y, const DP &z) const {
		// get grayvalue
		// return interp_lin(x, y, z);
		return interp_cub(x, y, z);
	}

	DP interp_lin(const DP &x, const DP &y, const DP &z) const {
		SI xi = SI(x), yi = SI(y), zi = SI(z);
		DP dx = x - xi, dy = y - yi, dz = z - zi;
		return
			DP(get_val(xi + 0, yi + 0, zi + 0)) * (1 - dx) * (1 - dy) * (1 - dz) +
			DP(get_val(xi + 1, yi + 0, zi + 0)) * (0 + dx) * (1 - dy) * (1 - dz) +
			DP(get_val(xi + 0, yi + 1, zi + 0)) * (1 - dx) * (0 + dy) * (1 - dz) +
			DP(get_val(xi + 1, yi + 1, zi + 0)) * (0 + dx) * (0 + dy) * (1 - dz) +
			DP(get_val(xi + 0, yi + 0, zi + 1)) * (1 - dx) * (1 - dy) * (0 + dz) +
			DP(get_val(xi + 1, yi + 0, zi + 1)) * (0 + dx) * (1 - dy) * (0 + dz) +
			DP(get_val(xi + 0, yi + 1, zi + 1)) * (1 - dx) * (0 + dy) * (0 + dz) +
			DP(get_val(xi + 1, yi + 1, zi + 1)) * (0 + dx) * (0 + dy) * (0 + dz) ;
	}

    DP interp_cub1D( DP vm, DP v0, DP v1, DP v2, DP x ) const {
        return
            vm * ( 0 - x ) * ( 1 - x ) * ( 2 - x ) / 6 +
            v0 * ( 1 + x ) * ( 1 - x ) * ( 2 - x ) / 2 +
            v1 * ( 1 + x ) * ( 0 + x ) * ( 2 - x ) / 2 +
            v2 * ( 1 + x ) * ( 0 + x ) * ( x - 1 ) / 6;
    }

	DP interp_cub(const DP &x, const DP &y, const DP &z) const {
		SI xi = SI(x), yi = SI(y), zi = SI(z);
		DP dx = x - xi, dy = y - yi, dz = z - zi;
		return interp_cub1D(
            interp_cub1D(
                interp_cub1D( get_val(xi - 1, yi - 1, zi - 1), get_val(xi - 1, yi + 0, zi - 1), get_val(xi - 1, yi + 1, zi - 1), get_val(xi - 1, yi + 2, zi - 1), dy ),
                interp_cub1D( get_val(xi + 0, yi - 1, zi - 1), get_val(xi + 0, yi + 0, zi - 1), get_val(xi + 0, yi + 1, zi - 1), get_val(xi + 0, yi + 2, zi - 1), dy ),
                interp_cub1D( get_val(xi + 1, yi - 1, zi - 1), get_val(xi + 1, yi + 0, zi - 1), get_val(xi + 1, yi + 1, zi - 1), get_val(xi + 1, yi + 2, zi - 1), dy ),
                interp_cub1D( get_val(xi + 2, yi - 1, zi - 1), get_val(xi + 2, yi + 0, zi - 1), get_val(xi + 2, yi + 1, zi - 1), get_val(xi + 2, yi + 2, zi - 1), dy ),
                dx
            ),
            interp_cub1D(
                interp_cub1D( get_val(xi - 1, yi - 1, zi + 0), get_val(xi - 1, yi + 0, zi + 0), get_val(xi - 1, yi + 1, zi + 0), get_val(xi - 1, yi + 2, zi + 0), dy ),
                interp_cub1D( get_val(xi + 0, yi - 1, zi + 0), get_val(xi + 0, yi + 0, zi + 0), get_val(xi + 0, yi + 1, zi + 0), get_val(xi + 0, yi + 2, zi + 0), dy ),
                interp_cub1D( get_val(xi + 1, yi - 1, zi + 0), get_val(xi + 1, yi + 0, zi + 0), get_val(xi + 1, yi + 1, zi + 0), get_val(xi + 1, yi + 2, zi + 0), dy ),
                interp_cub1D( get_val(xi + 2, yi - 1, zi + 0), get_val(xi + 2, yi + 0, zi + 0), get_val(xi + 2, yi + 1, zi + 0), get_val(xi + 2, yi + 2, zi + 0), dy ),
                dx
            ),
            interp_cub1D(
                interp_cub1D( get_val(xi - 1, yi - 1, zi + 1), get_val(xi - 1, yi + 0, zi + 1), get_val(xi - 1, yi + 1, zi + 1), get_val(xi - 1, yi + 2, zi + 1), dy ),
                interp_cub1D( get_val(xi + 0, yi - 1, zi + 1), get_val(xi + 0, yi + 0, zi + 1), get_val(xi + 0, yi + 1, zi + 1), get_val(xi + 0, yi + 2, zi + 1), dy ),
                interp_cub1D( get_val(xi + 1, yi - 1, zi + 1), get_val(xi + 1, yi + 0, zi + 1), get_val(xi + 1, yi + 1, zi + 1), get_val(xi + 1, yi + 2, zi + 1), dy ),
                interp_cub1D( get_val(xi + 2, yi - 1, zi + 1), get_val(xi + 2, yi + 0, zi + 1), get_val(xi + 2, yi + 1, zi + 1), get_val(xi + 2, yi + 2, zi + 1), dy ),
                dx
            ),
            interp_cub1D(
                interp_cub1D( get_val(xi - 1, yi - 1, zi + 2), get_val(xi - 1, yi + 0, zi + 2), get_val(xi - 1, yi + 1, zi + 2), get_val(xi - 1, yi + 2, zi + 2), dy ),
                interp_cub1D( get_val(xi + 0, yi - 1, zi + 2), get_val(xi + 0, yi + 0, zi + 2), get_val(xi + 0, yi + 1, zi + 2), get_val(xi + 0, yi + 2, zi + 2), dy ),
                interp_cub1D( get_val(xi + 1, yi - 1, zi + 2), get_val(xi + 1, yi + 0, zi + 2), get_val(xi + 1, yi + 1, zi + 2), get_val(xi + 1, yi + 2, zi + 2), dy ),
                interp_cub1D( get_val(xi + 2, yi - 1, zi + 2), get_val(xi + 2, yi + 0, zi + 2), get_val(xi + 2, yi + 1, zi + 2), get_val(xi + 2, yi + 2, zi + 2), dy ),
                dx
            ),
            dz );
	}

	void get_grad(DP *gra, const SI &o) const {
		// at integer pixel location gradient
		SI x, y, z;
		ind2sub(y, x, z, o);
		gra[0] = 0.5*(get_val(x + 1, y + 0, z + 0) - get_val(x - 1, y - 0, z - 0));
		gra[1] = 0.5*(get_val(x + 0, y + 1, z + 0) - get_val(x - 0, y - 1, z - 0));
		gra[2] = 0.5*(get_val(x + 0, y + 0, z + 1) - get_val(x - 0, y - 0, z - 1));
		gra[0] = (isnan(gra[0])) ? 0 : gra[0];
		gra[1] = (isnan(gra[1])) ? 0 : gra[1];
		gra[2] = (isnan(gra[2])) ? 0 : gra[2];
	}
};

