#pragma once


namespace Const {//physical conctants
    constexpr double Joul = 4184;//Joul-Calorie conversion
    constexpr double Avogadro = 6.022E+23;
    constexpr double Boltzmann = 0.0138065;
    constexpr double pNnm = 1.0e21;
}

namespace Condition {//Experiment conditions
    constexpr double Temperature = 298.15;
    constexpr double kT = Temperature * Const::Boltzmann;
    constexpr int ArmLength = 2200;//total length of the 2 dsDNA arms, unit is base-pair.
    constexpr double PillarStiffness = 0.07406;//spring constant of the pillar/optical trap/micro-needle/etc that is used for stretching.
}