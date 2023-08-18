#pragma once

#include <array>

#include "../shared_libs/math_lib.h"
#include "../shared_libs/constants_lib.h"
#include "../shared_libs/models_lib.h"

//==========================================constexpr LUT============================================
//define resolution and size of the look-up tables
//j = unzipped #bp
//ext = an index to indicate the total length of the system

#ifndef J_SIZE
#define J_SIZE 200
#endif

#ifndef J_RESELUTION
#define J_RESELUTION 8000/J_SIZE
#endif

#ifndef EXT_SIZE
#define EXT_SIZE 200
#endif

#ifndef EXT_RESELUTION
#define EXT_RESELUTION 8000/EXT_SIZE
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

namespace {

    using lut_type = std::array<std::array<float,ext_size>,j_size>;

    //=================================utility functions to create constexpr lut=================================
    constexpr double lz_ds (double force) {//dsDNA's length per base
        return Condition::ArmLength * DNAParams::L0DS * 
                alpha2phi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT);
    }
    constexpr double lz_ss (double force, int j) {//ssDNA's length per base
        return 2.0 * j * DNAParams::L0SS * 
                alpha2phi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT);
    }
    constexpr double le_ds (double force) {//function version of dsDNA's energy per bp:
        return Condition::ArmLength * Condition::kT * DNAParams::L0DS * 
                integ_alphadphi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT) / DNAParams::LPDS;
    }
    constexpr double le_ss (double force, int j) {//function version of ssDNA's energy per bp:
        return 2.0 * j * Condition::kT * DNAParams::L0SS * 
                integ_alphadphi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT) / DNAParams::LPSS;
    }
    constexpr double delta_ext(double force, double j, double ext) {//func used to find force so the system total extension = ext
        return ext - force/Condition::PillarStiffness - lz_ds (force) - lz_ss (force, j);//increasing function with force
    }

    constexpr bool False = false;
    constexpr double tor_binary_search = 1.0e-10;
    constexpr double find_force(int j, double ext) {//j == length of unzipped trunk, ext total extension
        //simple binary search to get force so calc_z_diff(force) == 0
        //force function must be monotonic
        double f1 = ValidRange::ValidMinForce;
        double f2 = ValidRange::ValidMaxForce;

        double y1 = delta_ext(f1, j, ext);
        double y2 = delta_ext(f2, j, ext);

        if (y1 * y2 >= 0) {
            if (y1 < 0){
                return ValidRange::ValidMaxForce + 1.0;//force is too large
            } else {
                return ValidRange::ValidMinForce - 1.0;//force is too small
            }
        }

        int cnt = 0;//in case the root is not found
        while (++cnt < 10000) {
            
            double fm = (f1 + f2) * 0.5;
            double ym = delta_ext(fm, j, ext);
            
            if (MyMath::Abs(ym) <= tor_binary_search) {
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
                return MyMath::VeryLargeNumber;//means some weird error
            }
        }
        return MyMath::VeryLargeNumber;//meaning that the root is not found
    }

}


constexpr lut_type Lut_force = []{

    lut_type arr;

    for (int j = 0; j < j_size; ++j) {
        for (int k = 0; k < ext_size; ++k) {
            arr[j][k] = find_force(j * j_resolution, k * ext_resolution);
        }
    }
    return arr;
}();

constexpr lut_type Lut_energy = []{

    lut_type arr;

    double f = 0.0;
    for (int j = 0; j < j_size; ++j) {
        for (int k = 0; k < ext_size; ++k) {
            f = Lut_force[j][k];
            if (f >= ValidRange::ValidMaxForce || f <= ValidRange::ValidMinForce) {
                arr[j][k] = MyMath::VeryLargeNumber;//make this a large number, meaning that do not use the value
            } else {
                arr[j][k] = (0.5 * f * f / Condition::PillarStiffness + le_ds(f) + le_ss(f, j * j_resolution))/Condition::kT;
            }
        }
    }

    return arr;
}();


class lookup_class {
public:
    constexpr lookup_class(lut_type lut_in) : lut(lut_in) {};

    constexpr double operator() (double j0, double extension) const {
        return lookup (j0, extension);
    };

    constexpr lut_type  get_lut() const {
        return lut;
    }

private:

    const lut_type lut;

    constexpr double lookup(double j0, double extension) const {

        if (j0 < 0 || j0 >= j_size * j_resolution || extension < 0 || extension >= ext_size * ext_resolution) {
            //error
            return -1.0;
        }

        double j = j0 / static_cast<double>(j_resolution);
        double k = extension / static_cast<double>(ext_resolution);
        
        //because j and k are positive.
        //static_cast<int> has the same result as std::floor() from cmath lib
        //Legend says that cast is 3 times faster
        double j1 = static_cast<int>(j);
        double k1 = static_cast<int>(k);

        double j2 = j1 + 1;
        double k2 = k1 + 1;

        //https://en.wikipedia.org/wiki/Bilinear_interpolation
        return lut[j1][k1] * (j2 - j) * (k2 - k) + 
                lut[j1][k2] * (j2 - j) * (k - k1) + 
                lut[j2][k1] * (j - j1) * (k2 - k) + 
                lut[j2][k2] * (j - j1) * (k - k1); 
        
    }
};


constexpr lookup_class Force {Lut_force};
constexpr lookup_class Energy {Lut_energy};



