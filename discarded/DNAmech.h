#pragma once

#include <array>
#include <limits>
#include <cmath>

#include "Math_constexpr.h"

//The constexpr sqrt func is copied from
//https://gist.github.com/alexshtf/eb5128b3e3e143187794
namespace {
    //limit access of these functions from outside of this file
    constexpr double sqrtNewtonRaphson(double x, double curr, double prev) {
        return curr == prev
            ? curr
            : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }
    /*
    * Constexpr version of the square root
    * Return value:
    *	- For a finite and non-negative value of "x", returns an approximation for the square root of "x"
    *   - Otherwise, returns NaN
    */
    constexpr double Sqrt(double x){
        return x >= 0 && x < std::numeric_limits<double>::infinity()
            ? sqrtNewtonRaphson(x, x, 0)
            : std::numeric_limits<double>::quiet_NaN();
    }

    constexpr double Square(double x){
        return x * x;
    }

}


//============================DNA parameters===========================
//These are not constants..and should be salt dependent and temperature dependent
//But I didn't consider these as dependent variables (for simplicity)

constexpr double LPDS = 51.97;//dsDNA persistence length
constexpr double KDS = 1318;//dsDNA elastic modulus
constexpr double L0DS = 0.338;//dsDNA contour length per bp
constexpr double LPSS = 0.765;//ssDNA persistence length
constexpr double KSS = 470;//ssDNA elastic modulus  
constexpr double L0SS = 0.554;//ssDNA contour length per nt

//==============================extensible models and ENERGY=======================================
//WLC high force, Odijk 1995 macromolecules
//phi = x/L (L = contour length)
//alpha = fA/kT (A = persistence length)
//k0_eff=k0A/kT (K0 = elastic modulus)
//for simplicity, use Odijk95 for all both ssDNA and dsDNA
constexpr double alpha2phi_Odijk95(double alpha, double k0_eff) { 
    return 1.0 - 0.5 / Sqrt(alpha) + alpha / k0_eff; 
}
constexpr double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
    return alpha - Sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
}
constexpr double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
    return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
}

// constexpr double alpha2phi_Smith95(double alpha, double k0_eff) {
//     return ((1.0 / Tanh(2.0 * alpha) - 0.5 / alpha) * (1.0 + alpha / k0_eff));
// }
// constexpr double integ_alphadphi_Smith95_numeric(double alpha, double k0_eff) {
//     //numerica integration of phi * d(alpha) for Smith95 model

//     double incre = alpha / 100.0;
//     double accum = 0.0;//the first value
//     for (double a = incre; a <= alpha; a += incre) {
//         accum += alpha2phi_Smith95(a, k0_eff);
//     }
//     accum -= 0.5 * alpha2phi_Smith95(alpha, k0_eff);
//     return alpha2phi_Smith95(alpha, k0_eff) * alpha - accum * incre;
// }


//thermal const
constexpr double kT = 4.1;

//define max/min/increment of force
constexpr double force_min = 0.01;
constexpr double force_max = 50.0;
constexpr double force_interval = 0.001;

//calculate the size for all the luts.
constexpr size_t lsize = static_cast<size_t>((force_max-force_min)/force_interval);

//lut of force
constexpr auto lforce = [] {
    std::array<double, lsize> arr = {};
    for (int i = 0; i < lsize; ++i){
        arr[i] = force_min + i * force_interval;
    }
    return arr;
}();

//lut of ssDNA's length per base
constexpr auto lz_ss = [] {
    std::array<double, lsize> arr = {};
    double alpha;
    double keff = KSS * LPSS / kT;
    for (int i = 0; i < lsize; ++i) {
        alpha = (force_min + i * force_interval) * LPSS / kT;
        arr[i] = L0SS * alpha2phi_Odijk95(alpha, keff);
        // arr[i] = L0SS * alpha2phi_Smith95(alpha, keff);
    }
    return arr;
}();

//lut of dsDNA's length per base
constexpr auto lz_ds = [] {
    std::array<double, lsize> arr = {};
    double alpha;
    double keff = KDS * LPDS / kT;
    for (int i = 0; i < lsize; ++i) {
        alpha = (force_min + i * force_interval) * LPDS / kT;
        arr[i] = L0DS * alpha2phi_Odijk95(alpha, keff);
    }
    return arr;
}();

//lut of ssDNA's energy per bp:
constexpr auto le_ss = [] {
    std::array<double, lsize> arr = {};
    double alpha;
    double keff = KSS * LPSS / kT;
    double a = kT * L0SS / LPSS;
    for (int i = 0; i < lsize; ++i) {
        alpha = (force_min + i * force_interval) * LPSS / kT;
        arr[i] = a * integ_alphadphi_Odijk95(alpha, keff);
        // arr[i] = a * integ_alphadphi_Smith95_numeric(alpha, keff);
    }
    return arr;
}();

//lut of dsDNA's energy per bp:
constexpr auto le_ds = [] {
    std::array<double, lsize> arr = {};
    double alpha;
    double keff = KDS * LPDS / kT;
    double a = kT * L0DS / LPDS;
    for (int i = 0; i < lsize; ++i) {
        alpha = (force_min + i * force_interval) * LPDS / kT;
        arr[i] = a * integ_alphadphi_Odijk95(alpha, keff);
    }
    return arr;
}();

static_assert(lforce.size() > 0);
static_assert(lz_ds.size() > 0);
static_assert(lz_ss.size() > 0);
static_assert(le_ds.size() > 0);
static_assert(le_ss.size() > 0);
static_assert(lforce[lsize-1]);
static_assert(lz_ds[lsize-1]);
static_assert(lz_ss[lsize-1]);
static_assert(le_ds[lsize-1]);
static_assert(le_ss[lsize-1]);

//not used, combined to lookup()
inline double f2idx(double force) {
    //returned value is a fractional index;
    return (force - force_min) * (lsize - 1) / (force_max - force_min); 
}

//not used, combined to lookup()
inline double interpolate(double idx, const decltype(lforce) & lut) {
    //use a fractional index to interpolate lut
    int i = std::floor(idx);
    return (idx - i) * lut.at(i + 1) + (i + 1 - idx) * lut.at(i); 
}

//lookup by lin interp
inline double lookup (const decltype(lforce) & lut, double force) {
    double idx = (force - force_min) * (lsize - 1) / (force_max - force_min);
    double i = std::floor(idx);
    if (i == idx) {//corner case
        return lut.at(i);
    }
    return (idx - i) * lut.at(i + 1) + (i + 1 - idx) * lut.at(i); 
}

//look up by quadratic interp
inline double lookup2nd (const decltype(lforce) & lut, double force) {
    double idx = (force - force_min) * (lsize - 1) / (force_max - force_min);
    if (std::floor(idx) == idx) {//corner case
        return lut.at(idx);
    }

    size_t i = static_cast<size_t> (idx + 0.5); //the index the idx closest to
    double delta = idx - i;

    double a = 0.5 * (lut.at(i + 1) + lut.at(i - 1) - 2.0 * lut.at(i));
    double b = 0.5 * (lut.at(i + 1) - lut.at(i - 1)) - 2.0 * a * i;

    return lut.at(i) + 2.0 * a * i * delta + a * delta * delta + b * delta;
}

//this test shows that quadratic interp is not much better than lin interp. the improvement is unimportant
// for (double f = 0.3; f < 5.0; f += 0.1){
//     double x1 = lookup(lz_ds, f);
//     double x2 = lookup2nd(lz_ds, f);
//     double x3 = L0DS * alpha2phi_Odijk95(f * LPDS / kT, KDS * LPDS / kT);
//     std::cout << x1 - x3 << '\t' << x2 - x3 << '\t' << x3 << std::endl;
// }

//====================================test func=================================
//function version of ssDNA's length per base
double lz_ss_f (double force) {
    return L0SS * alpha2phi_Odijk95(force * LPSS / kT, KSS * LPSS / kT);
    // return L0SS * alpha2phi_Smith95(force * LPSS / kT, KSS * LPSS / kT);
}

//function version of dsDNA's length per base
double lz_ds_f (double force) {
    return L0DS * alpha2phi_Odijk95(force * LPDS / kT, KDS * LPDS / kT);
}

//function version of ssDNA's energy per bp:
double le_ss_f (double force) {
    return kT * L0SS * integ_alphadphi_Odijk95(force * LPSS / kT, KSS * LPSS / kT) / LPSS;
    // return kT * L0SS * integ_alphadphi_Smith95_numeric(force * LPSS / kT, KSS * LPSS / kT) / LPSS;
}

//function version of dsDNA's energy per bp:
double le_ds_f (double force) {
    return kT * L0DS * integ_alphadphi_Odijk95(force * LPDS / kT, KDS * LPDS / kT) / LPDS;
}

