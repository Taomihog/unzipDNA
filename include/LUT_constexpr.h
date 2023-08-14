#pragma once

#include <array>

#include "Math_constexpr.h"
#include "Constant_constexpr.h"
#include "DNAmech_constexpr.h"

//==========================================constexpr LUT============================================
//define resolution and size of lut_f and lut_e
//i = arm_length in bp (not used now, otherwise the lut is too big)
//j = unzipped #bp
//k = an index to indicate the total length of the system

//use a pre-calucated constexpr LUT will speed up the calculation by at least 1 order!!
constexpr int j_size = 200;//size of the j index dimension of the LUT. basically the max length of trunk
constexpr int ext_size = 200;//size of the extension dimension of LUT, the larger the number, the more precise the result.
constexpr int ext_resolution = 40;//resolution of the total extension dimention: actual total extension is index * resolution
constexpr int j_resolution = 40;//resolution of the j-index dimention: actual j is index * resolution

//coarse settings for testing purposes
// constexpr int j_size = 20;
// constexpr int ext_size = 20;
// constexpr int ext_resolution = 400;
// constexpr int j_resolution = 400;

using lut_type = std::array<std::array<double,ext_size>,j_size>;


//=================================utility functions to create constexpr lut=================================
constexpr double lz_ss (double force, int j) {
//ssDNA's length per base
	return 2 * j * Const::L0SS * alpha2phi_Smith95_hf(force * Const::LPSS / Const::kT, Const::KSS * Const::LPSS / Const::kT);
}

constexpr double lz_ds (double force) {
//dsDNA's length per base
	return Const::ArmLength * Const::L0DS * alpha2phi_Odijk95(force * Const::LPDS / Const::kT, Const::KDS * Const::LPDS / Const::kT);
}

constexpr double le_ss (double force, int j) {
//function version of ssDNA's energy per bp:
	return 2 * j * Const::kT * Const::L0SS * integ_alphadphi_Odijk95(force * Const::LPSS / Const::kT, Const::KSS * Const::LPSS / Const::kT) / Const::LPSS;
}

constexpr double le_ds (double force) {
//function version of dsDNA's energy per bp:
	return Const::ArmLength * Const::kT * Const::L0DS * integ_alphadphi_Smith95_hf(force * Const::LPDS / Const::kT, Const::KDS * Const::LPDS / Const::kT) / Const::LPDS;
}

//=====================find the force for certain total extension======================
constexpr double delta_ext(double force, double j, double ext) {//func used to find force so the system total extension = ext
	return ext - force/Const::PillarStiffness - lz_ds (force) - lz_ss (force, j);
}

constexpr double find_force(int j, double ext) {//j == length of unzipped trunk, ext total extension
	//simple binary search to get force so calc_z_diff(force) == 0
	//force function must be monotonic
	double f1 = 0.00001;
	double f2 = 10000.0;

double y1 = delta_ext(f1, j, ext);
double y2 = delta_ext(f2, j, ext);

    if (y1 * y2 >= 0) {
        return 0.0;//force is out of the range (f1,f2), return 0.0 (usually this is due to the force less than f1)
    }

	int cnt = 0;//in case the root is not found
	while (++cnt < 10000) {
		
		double fm = (f1 + f2) * 0.5;
		double ym = delta_ext(fm, j, ext);
		
		if (Math_constexpr::Abs(ym) <= tor_binary_search) {
			return f2;
		}

		//(&& has higher precedence)
		if (y1 < y2 && ym > 0 || y1 > y2 && ym < 0) {
			f2 = fm;
			y2 = ym;
		} else if (y1 < y2 && ym < 0 || y1 > y2 && ym > 0) {
			f1 = fm;
			y1 = ym;
		} else {
			return Inf;//means some weird error
		}
	}

	return Inf;//meaning that the root is not found
}

constexpr auto Lut_force = []{
	
	lut_type arr;

	for (int j = 0; j < j_size; ++j) {
		for (int k = 0; k < ext_size; ++k) {
			arr[j][k] = find_force(j * j_resolution, k * ext_resolution);
		}
	}
	return arr;
}();

constexpr auto Lut_energy = []{
	
	lut_type arr;

	double f = 0.0;
	for (int j = 0; j < j_size; ++j) {
		for (int k = 0; k < ext_size; ++k) {
			f = Lut_force[j][k];
			arr[j][k] = 0.5 * f * f / Const::PillarStiffness + le_ds(f) + le_ss(f, j * j_resolution);
		}
	}

	return arr;
}();
