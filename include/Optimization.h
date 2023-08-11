#pragma once

#include <limits>
#include <functional>

//why do not use lambda:
//=============Passing capturing lambda as function pointer====================
//https://stackoverflow.com/questions/28746744/passing-capturing-lambda-as-function-pointer

namespace {
	constexpr double inf = std::numeric_limits<double>::infinity(); 
	constexpr double tor0 = 1.0e-12;
}

namespace Optimization {
	// a simple bisection method which only works on monotonic functions
	double Bisection(double, double, std::function<double(double)>, double tor = tor0);
}