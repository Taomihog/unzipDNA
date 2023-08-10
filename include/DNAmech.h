#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <limits>

using namespace std;

//https://gist.github.com/alexshtf/eb5128b3e3e143187794
namespace {
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

    constexpr double reciprocal(double x) {
        return 1.0 / x;
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
}

//thermal const
constexpr double kT = 4.1;

//define max/min/increment of force
constexpr double fmin = 5.0;
constexpr double fmax = 25.0;
constexpr double finterval = 1.0;

//calculate the size for all the luts.
constexpr size_t lsize = static_cast<size_t>((fmax-fmin)/finterval);

//lut of force
constexpr auto lforce = [] {
    std::array<double, lsize> arr = {};
    for (int i = 0; i < lsize; ++i){
        arr[i] = fmin + i * finterval;
    }
    return arr;
}();

//lut of ssDNA's length per base
constexpr auto lz_ss = [] {
    std::array<double, lsize> arr = {};
    double alpha;
    double keff = KSS * LPSS / kT;
    for (int i = 0; i < lsize; ++i) {
        alpha = (fmin + i * finterval) * LPSS / kT;
        arr[i] = L0SS * alpha2phi_Odijk95(alpha, keff);
    }
    return arr;
}();

//lut of dsDNA's length per base
constexpr auto lz_ds = [] {
    std::array<double, lsize> arr = {};
    double alpha;
    double keff = KDS * LPDS / kT;
    for (int i = 0; i < lsize; ++i) {
        alpha = (fmin + i * finterval) * LPDS / kT;
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
        alpha = (fmin + i * finterval) * LPSS / kT;
        arr[i] = a * integ_alphadphi_Odijk95(alpha, keff);
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
        alpha = (fmin + i * finterval) * LPDS / kT;
        arr[i] = a * integ_alphadphi_Odijk95(alpha, keff);
    }
    return arr;
}();