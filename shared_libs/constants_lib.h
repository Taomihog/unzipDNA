#pragma once

//Constants
namespace Const {
    constexpr double Joul = 4184;//Joul-Calorie conversion
    constexpr double Avogadro = 6.022E+23;
    constexpr double Boltzmann = 0.0138065;
    constexpr double pNnm = 1.0e21;
}

//Experiment conditions
namespace Condition {
    constexpr double Temperature = 298.15;
    constexpr double kT = Temperature * Const::Boltzmann;
    constexpr int ArmLength = 2200;//total length of the 2 dsDNA arms, unit is base-pair.
    constexpr double PillarStiffness = 0.07406;//spring constant of the pillar/optical trap/micro-needle/etc that is used for stretching.
    constexpr double SaltConc = 100; //salt concentration in mM, this the the salt concentraion or the unzip experiment;
}

//====================================DNA mechanical parameters======================================
//These are pseudo constants..they are salt dependent and temperature dependent
namespace DNAParams{
    constexpr double LPDS = 51.97;//dsDNA persistence length
    constexpr double KDS = 1318;//dsDNA elastic modulus
    constexpr double L0DS = 0.338;//dsDNA contour length per bp
    constexpr double LPSS = 0.765;//ssDNA persistence length
    constexpr double KSS = 470;//ssDNA elastic modulus  
    constexpr double L0SS = 0.554;//ssDNA contour length per nt
}

//====================================basepair energy measured by Huguet et al===============================
//ref: Huguet, Josep M., et al. (2010) PNAS
namespace BPEnergy {
    

    constexpr int LUTsize = 4;

    constexpr double LUTdH[LUTsize][LUTsize] = {
        {-7.28, -4.63, -5.21, -5.80},//aa, at, ag, ac
        {-8.31, -7.28, -8.96, -8.16},//ta, tt, tg, tc
        {-8.16, -5.80, -8.57, -10.1},//ga, gt, gg, gc
        {-8.96, -5.21, -9.66, -8.57} //ca, ct, cg, cc
    };

    constexpr double LUTdS[LUTsize][LUTsize] = {
        {-20.28, -11.62, -12.89, -14.46},//aa, at, ag, ac
        {-25.06, -20.28, -24.48, -22.46},//ta, tt, tg, tc
        {-22.46, -14.46, -22.30, -25.96},//ga, gt, gg, gc
        {-24.48, -12.89, -24.43, -22.30} //ca, ct, cg, cc
    };

    constexpr double LUTm[LUTsize][LUTsize] = {
        {0.145, 0.117, 0.070, 0.099},//aa, at, ag, ac
        {0.091, 0.145, 0.091, 0.155},//ta, tt, tg, tc
        {0.155, 0.099, 0.063, 0.079},//ga, gt, gg, gc
        {0.091, 0.070, 0.132, 0.063} //ca, ct, cg, cc
    };

    //Constexpr class
    //https://stackoverflow.com/questions/61281843/creating-compile-time-key-value-map-in-c
    class constexpr_map_class {

        public:
        constexpr int operator[] (char key) const {
                return bp2idx (key);
        }

        private:


        constexpr int bp2idx(char base) const {
            switch (base) {
            case 'a':
            case 'A':
                return 0;

            case 't':
            case 'T':
                return 1;

            case 'g':
            case 'G':
                return 2;

            case 'c':
            case 'C':
                return 3;

            default:
                return -1;
            }
        };
    };

    constexpr constexpr_map_class bp2idx_map;

    static_assert (bp2idx_map['a'] == 0, "Error.");
    static_assert (bp2idx_map['A'] == 0, "Error.");
    static_assert (bp2idx_map['t'] == 1, "Error.");
    static_assert (bp2idx_map['T'] == 1, "Error.");
    static_assert (bp2idx_map['g'] == 2, "Error.");
    static_assert (bp2idx_map['G'] == 2, "Error.");
    static_assert (bp2idx_map['c'] == 3, "Error.");
    static_assert (bp2idx_map['C'] == 3, "Error.");
    //static_assert (bp2idx_map['C'] == 4, "Error.");//will cause error
}