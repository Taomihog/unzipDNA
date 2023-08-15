#pragma once

#include "Math_constexpr.h"
using namespace Math_constexpr;

//the DNA structure is something like this:

//             ↑Force
//             ||
//             ||
//               =================
//             ||
//             ||
//             ↓Force

//This is called a "Y structure". A Y structure has 2 dsDNA arms and a dsDNA trunk.
//when the two arms are stretched as shown, the trunk will be unzipped.
//the force during this "unzip" precess is measured against the end-to-end distance of the arms. 
//The sequence of the trunk determines the "force vs extension" profile.
//And can be calculated theoretically.

//for more information, see:

//[1] B. Essevaz-Roulet, U. Bockelmann, and F. Heslot (1997) PNAS  
//[2] Bockelmann, Ulrich, et al. (2002)Biophysical journal 
//[3] Huguet, Josep M., et al. (2010) PNAS


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
    return (1.0 / Tanh(2.0 * alpha) - 0.5 / alpha) * (1.0 + alpha / k0_eff);
}
constexpr double alpha2phi_Smith95_hf(double alpha, double k0_eff) {
    //high force version, the tan(..)  term equals to 1
    return (1.0 - 0.5 / alpha) * (1.0 + alpha / k0_eff);
}
constexpr double integ_phidalpha_Smith95_hf(double alpha, double k0_eff) { 
    //high force version, the tan(..)  term equals to 1
    return (0.5 * alpha + k0_eff - 0.5)  * alpha / k0_eff - 0.5 * Ln(alpha);
}
constexpr double integ_alphadphi_Smith95_hf(double alpha, double k0_eff) {
    return alpha * alpha2phi_Smith95_hf(alpha, k0_eff) - integ_phidalpha_Smith95_hf(alpha, k0_eff);
}
