#pragma once

#include "math_lib.h"

//================================WLC/FJC model: parameter definitions==================================
//phi = x/L (L = contour length)
//alpha = fA/kT (A = persistence length)
//k0_eff=k0A/kT (K0 = elastic modulus)

//======================================valid force range for the model==========================================
//to increase speed and to make the models constexpr feasible, I simplified the models
//the modified models are only precious in a certain range
//energy/force calculation above the range may not be accurate
namespace ValidRange{
    constexpr double ValidMaxForce = 10000000.0;
    constexpr double ValidMinForce = 0.2;
}


//====================================Marko-Siggia 1995 WLC========================================
constexpr double phi2alpha_MS(double phi){ return phi + 0.25 / MyMath::Square(1.0 - phi) - 0.25; }//need to explicitly declare "inline" functions
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
    // if (alpha < 0.25) {
    //     return 0.0;// undefined at alpha == 0, just give it a small value
    //     //I can do this because this is force, energy must be calculated correctly!!
    // }
    return 1.0 - 0.5 / MyMath::Sqrt(alpha) + alpha / k0_eff; 
}
constexpr double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
    return alpha - MyMath::Sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
}
constexpr double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
    return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
}
// ================================MODIFIED VERSION OF FJC, Smith 1995 macromolecules================================
//Modified version specific for ssDNA force region, and keeps accuracy
//For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
//My homemade Langevin_integ function should be accurate enough in this region.
constexpr double alpha2phi_Smith95_m(double alpha, double k0_eff) {//"m" means modified
    return MyMath::Langevin(2.0 * alpha) + alpha / k0_eff;
}
constexpr double integ_phidalpha_Smith95_m(double alpha, double k0_eff) { 
    return 0.5 * MyMath::Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
}
constexpr double integ_alphadphi_Smith95_m(double alpha, double k0_eff) {
    //integ actually starts from 1, but it's OK since it is for partition function calculation
    return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
}

//================================FJC, Smith 1995 macromolecules================================
//the resolution is not very good, and the calculation is too slow
// constexpr double alpha2phi_Smith95(double alpha, double k0_eff) {
//     return (MyMath::Coth(2.0 * alpha) - 0.5 / alpha) * (1.0 + alpha / k0_eff);
// }
// constexpr double integ_phidalpha_Smith95(double alpha, double k0_eff) { 
//     const int n = 100;//integration from alpha/n to alpha!!!, not from 0 so there is a minor error
//     const double delta = alpha/n;
//     double sum = 0.0;
//     double a = 0.0;
//     for (int i = 1; i < n + 1; ++i) {//means i can = n.
//         a = delta * i;
//         sum += alpha2phi_Smith95(a, k0_eff);
//     }
//     sum -= 0.5 * alpha2phi_Smith95(delta, k0_eff);
//     sum -= 0.5 * alpha2phi_Smith95(alpha, k0_eff);
//     return sum * delta;
// }
// constexpr double integ_alphadphi_Smith95(double alpha, double k0_eff) {
//     //integ actually starts from 1, but it's OK since it is for partition function calculation
//     return alpha * alpha2phi_Smith95(alpha, k0_eff) - integ_phidalpha_Smith95(alpha, k0_eff);
// }

// ================================High force OF FJC, Smith 1995 macromolecules================================
// The langevian function is ill formed (at least for the computer) when alpha -> 0 so avoid it=========
// These are high force version of Smith's extensible FJC, the coth(..)  term equals to 1 when force is large
// For ssDNA, since we also calculate dE/dJ, the low force region energy contribution thus becomes important, 
// we cannot use this high force approximation

// constexpr double alpha2phi_Smith95_hf(double alpha, double k0_eff) {
//     return (1.0 - 0.5 / alpha) * (1.0 + alpha / k0_eff);
// }
// constexpr double integ_phidalpha_Smith95_hf(double alpha, double k0_eff) { 
//     //remove the error of integration from zero to ValidMinAlpha and replace it with Langevin_0_to_MinAlpha
//     return (0.5 * alpha + k0_eff - 0.5)  * alpha / k0_eff - 0.5 * MyMath::Ln(alpha);
// }
// constexpr double integ_alphadphi_Smith95_hf(double alpha, double k0_eff) {
//     return alpha * alpha2phi_Smith95_hf(alpha, k0_eff) - integ_phidalpha_Smith95_hf(alpha, k0_eff);
// }
