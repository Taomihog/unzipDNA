#pragma once

#include "Math_constexpr.h"
using namespace Math_constexpr;

//====================================DNA mechanical parameters======================================
//These are pseudo constants..they are salt dependent and temperature dependent

constexpr double LPDS = 51.97;//dsDNA persistence length
constexpr double KDS = 1318;//dsDNA elastic modulus
constexpr double L0DS = 0.338;//dsDNA contour length per bp
constexpr double LPSS = 0.765;//ssDNA persistence length
constexpr double KSS = 470;//ssDNA elastic modulus  
constexpr double L0SS = 0.554;//ssDNA contour length per nt

//================================WLC/FJC model: parameter definitions==================================
//phi = x/L (L = contour length)
//alpha = fA/kT (A = persistence length)
//k0_eff=k0A/kT (K0 = elastic modulus)

//====================================Marko-Siggia 1995 WLC========================================
constexpr double phi2alpha_MS(double phi){ return phi + 0.25 / Square(1.0 - phi) - 0.25; }//need to explicitly declare "inline" functions
//For now I haven't implement constexpr cos() and acos(), so I cannot make this function constexpr.
// double alpha2phi_MS(double alpha) {//My ppt: xxxxxxxx
    
//     alpha = alpha - 0.75;
//     double p = -Square(alpha) / 3.0;
//     double q = -2.0 * Cubic(alpha) / 27.0 + 0.25;
    
//     //Cadano's formula is only correct when fA/kT-0.75 < 3.0/np.cbrt(16.0)
//     if (alpha < 3.0/Cbrt(16.0)) {
//         double m = Sqrt(Square(q)/4.0 + Cubic(p)/27.0);
//         double n= -q/2.0;
//         return 1.0 + alpha/3.0 + Cbrt(n + m) + Cbrt(n - m);
//     } else { //use trigonometric solution
//         return 1.0 + alpha / 3.0 + 2.0 * Sqrt(-p / 3.0) * cos(acos(1.5 * q * Sqrt(-3.0 / p) / p) / 3.0 - 2.0 * Pi * 2.0 / 3.0);
//     }
// }

//==========================extensible MS, or Modified-MS model, Wang et al. 1997===============================
//NO additional function to calculate x from f, because if we can substitute the x/L-f/K by:
//phi = x/L-f/K and alpha = fA/kT
//Then we can calculate:
//constexpr double alpha2phi_MMS(){}
//need to finish alpha2phi_MS first
//todo: phi2alpha_MMS() will definitely need extra work.

//==============================WLC high force, Odijk 1995 macromolecules===========================
constexpr double alpha2phi_Odijk95(double alpha, double k0_eff) { 
    
    return 1.0 - 0.5 / Sqrt(alpha) + alpha / k0_eff; 
}
constexpr double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
    return alpha - Sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
}
constexpr double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
    return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
}

//================================FJC, Smith 1995 macromolecules================================
constexpr double alpha2phi_Smith95(double alpha, double k0_eff) {
    return (Coth(2.0 * alpha) - 0.5 / alpha) * (1.0 + alpha / k0_eff);
}
constexpr double integ_phidalpha_Smith95(double alpha, double k0_eff) { 
    const int n = 100;//integration from alpha/n to alpha!!!, not from 0 so there is a minor error
    const double delta = alpha/n;
    double sum = 0.0;
    double a = 0.0;
    for (int i = 1; i < n + 1; ++i) {//means i can = n.
        a = delta * i;
        sum += alpha2phi_Smith95(a, k0_eff);
    }
    sum -= 0.5 * alpha2phi_Smith95(delta, k0_eff);
    sum -= 0.5 * alpha2phi_Smith95(alpha, k0_eff);
    return sum * delta;
}
constexpr double integ_alphadphi_Smith95(double alpha, double k0_eff) {
    //integ actually starts from 1, but it's OK since it is for partition function calculation
    return alpha * alpha2phi_Smith95(alpha, k0_eff) - integ_phidalpha_Smith95(alpha, k0_eff);
}

// constexpr double alpha2phi_Smith95_hf(double alpha, double k0_eff) {
//     //high force version, the tan(..)  term equals to 1
//     return (1.0 - 0.5 / alpha) * (1.0 + alpha / k0_eff);
// }
// constexpr double integ_phidalpha_Smith95_hf(double alpha, double k0_eff) { 
//     //high force version, the tan(..)  term equals to 1
//     return (0.5 * alpha + k0_eff - 0.5)  * alpha / k0_eff - 0.5 * Ln(alpha);
// }
// constexpr double integ_alphadphi_Smith95_hf(double alpha, double k0_eff) {
//     //integ actually starts from 1, but it's OK since it is for partition function calculation
//     return alpha * alpha2phi_Smith95_hf(alpha, k0_eff) - integ_phidalpha_Smith95_hf(alpha, k0_eff);
// }
