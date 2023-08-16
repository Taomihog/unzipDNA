#pragma once

#include <array>

#include "../include/Math_constexpr.h"
#include "../include/Constant_constexpr.h"
#include "../include/DNAmech_constexpr.h"

//==========================================constexpr LUT============================================
//define resolution and size of the look-up tables
//j = unzipped #bp
//ext = an index to indicate the total length of the system

#ifndef J_SIZE
#define J_SIZE 256
#endif

#ifndef J_RESELUTION
#define J_RESELUTION 8192/J_SIZE
#endif

#ifndef EXT_SIZE
#define EXT_SIZE 256
#endif

#ifndef EXT_RESELUTION
#define EXT_RESELUTION 8192/EXT_SIZE
#endif

//use a pre-calucated constexpr LUT will speed up the calculation by at least 1 order!!
constexpr int j_size = J_SIZE;//size of the j index dimension of the LUT. basically the max length of trunk
constexpr int ext_size = EXT_SIZE;//size of the extension dimension of LUT, the larger the number, the more precise the result.
constexpr int ext_resolution = J_RESELUTION;//resolution of the total extension dimention: actual total extension is index * resolution
constexpr int j_resolution = EXT_RESELUTION;//resolution of the j-index dimention: actual j is index * resolution

//coarse settings for testing purposes
// constexpr int j_size = 20;
// constexpr int ext_size = 20;
// constexpr int ext_resolution = 400;
// constexpr int j_resolution = 400;

using lut_type = std::array<std::array<double,ext_size>,j_size>;

//=================================utility functions to create constexpr lut=================================
constexpr double lz_ss (double force, int j) {
//ssDNA's length per base
    return 2.0 * j * L0SS * alpha2phi_Smith95(force * LPSS / Condition::kT, KSS * LPSS / Condition::kT);
}

constexpr double lz_ds (double force) {
//dsDNA's length per base
    return Condition::ArmLength * L0DS * alpha2phi_Odijk95(force * LPDS / Condition::kT, KDS * LPDS / Condition::kT);
}

constexpr double le_ss (double force, int j) {
//function version of ssDNA's energy per bp:
    return 2.0 * j * Condition::kT * L0SS * integ_alphadphi_Smith95(force * LPSS / Condition::kT, KSS * LPSS / Condition::kT) / LPSS;
}

constexpr double le_ds (double force) {
//function version of dsDNA's energy per bp:
    return Condition::ArmLength * Condition::kT * L0DS * integ_alphadphi_Odijk95(force * LPDS / Condition::kT, KDS * LPDS / Condition::kT) / LPDS;
}

//=====================find the force for certain total extension======================
constexpr double delta_ext(double force, double j, double ext) {//func used to find force so the system total extension = ext
    return ext - force/Condition::PillarStiffness - lz_ds (force) - lz_ss (force, j);
}

constexpr double find_force(int j, double ext) {//j == length of unzipped trunk, ext total extension
    //simple binary search to get force so calc_z_diff(force) == 0
    //force function must be monotonic
    double f1 = 0.001;
    double f2 = 1000.0;

    double y1 = delta_ext(f1, j, ext);
    double y2 = delta_ext(f2, j, ext);

    if (y1 * y2 >= 0) {
        if (y1 > 0){
            return f1;//force is too small
        } else {
            return f2;//force is too large
        }
    }

	int cnt = 0;//in case the root is not found
    while (++cnt < 10000) {
        
        double fm = (f1 + f2) * 0.5;
        double ym = delta_ext(fm, j, ext);
        
        if (Math_constexpr::Abs(ym) <= tor_binary_search) {
            return fm;
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
            arr[j][k] = (0.5 * f * f / Condition::PillarStiffness + le_ds(f) + le_ss(f, j * j_resolution))/Condition::kT;
        }
    }

    return arr;
}();



