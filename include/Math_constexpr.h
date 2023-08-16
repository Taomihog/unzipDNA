#pragma once

#include <limits>
#include <array>


namespace Math_constexpr{
    
    constexpr double Inf = std::numeric_limits<double>::infinity(); 
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double Pi = 3.1415926535897932384626433832795;
    constexpr size_t precision_ln = 30; 
    constexpr double tor_binary_search = 1.0e-10;

    //Newton-Raphson method for root calculation
    constexpr double Sqrt(double x) {
        //https://gist.github.com/alexshtf/eb5128b3e3e143187794

        if (x < 0) {
            return NaN;
        }

        if (x == 0) {
            return 0.0;
        }

        double prev = 1.0;
        double curr = x;
        while (curr != prev) {
            prev = curr;
            curr = 0.5 * (prev + x / prev);
        }

        return curr;
    }
    
    constexpr double Cbrt(double x) {
        //https://math.stackexchange.com/questions/3586214/newtons-method-for-cube-root
        
        if (x == 0) {
            return 0.0;
        }

        double prev = 1.0;
        double curr = x;

        while (curr != prev) {
            prev = curr;
            curr = 0.3333333333333333333333333333333333333 * (2.0 * prev + x / (prev * prev));
        }

        return curr;
    }

    constexpr double Pow(double x, int n) {
        double res = 1;
        for (int i = 0; i < n; ++i){
            res *= x;
        }
        return res;
    }

    constexpr double Square(double x){
        return x * x;
    }
    
    constexpr double Cubic(double x) {
        return x * x * x;
    }

    constexpr double Abs(double x) {
        return x > 0.0 ? x : -x;
    }

    constexpr double Tanh(double x) {
        //the polyn formula is copied from 
        //https://math.stackexchange.com/questions/107292/rapid-approximation-of-tanhx
        return 
        (-.67436811832e-5+(.2468149110712040+(.583691066395175e-1+.3357335044280075e-1*x)*x)*x)/
        (.2464845986383725+(.609347197060491e-1+(.1086202599228572+.2874707922475963e-1*x)*x)*x);
    }

    constexpr double Ln(double x) {
        //I wrote this! Haha:)
        
        if (x < 0.0) {
            return 1.0 / 0.0;//err;
        }
        
        if (x < 0.25) {
            return -Ln(1.0 / x);
        }
        if (x > 4.0 ) {;
            return Ln(x * 0.5) + 0.69314718055994530941723212145818;
        }
        //calculate x if 1.0 <= x <= 4.0
        double y = (x - 1) / (x + 1);
        const double ysq = y * y;
        double ypow = y; //y^+1
        double sum = ypow;
        for (size_t i = 1; i <= precision_ln; ++i) {
            ypow *= ysq;
            sum += ypow/(2*i + 1);
        }
        return 2.0 * sum;
    }

    
    constexpr double Coth(double x) {
        if ( x > 7.9608220) {
            // I compared Coth(x) and the np.cosh(x)/np.sinh(x), and 1, and figured out this value.
            //above this value, 1.0 is more close to Coth(x).
            //this value depends on how many terms are used of course.
            return 1.0;
        }
        double product = 1.0;//x**0/0!
        double sum1 = product;
        double sum2 = 0.0;
        int i = 1;
        for (; i <= 15; ++i){
            product *= x;
            product /= 2*i - 1;//x**1/1!, ..
            sum2 += product;
            product *= x;
            product /= 2*i; //x**2/2!, ...
            sum1 += product;
        }
        product *= x;
        product /= 2*i + 1;//x**1/1!, ..
        sum2 += product;
        
        return sum1/sum2;
    }


    //some tests

    constexpr double Coth_test_res = Coth(0.9);//1.3960672530300118350929600819912
    static_assert(Coth_test_res > 1.3960672530);//1.3960672525087865
    static_assert(Coth_test_res < 1.3960672531);

    constexpr double Sqrt_test_res = Sqrt(4.0);
    static_assert(Sqrt_test_res > 1.9999999999);
    static_assert(Sqrt_test_res < 2.0000000001);

    constexpr double Cbrt_test_res = Cbrt(8.0);
    static_assert(Cbrt_test_res > 1.9999999999);
    static_assert(Cbrt_test_res < 2.0000000001);

    constexpr double Tanh_test_res = Tanh(0.4);
    static_assert(Tanh_test_res > 0.379922);//This is not what tanh(0.4) really is, this is just a test to make sure it return a valuel
    static_assert(Tanh_test_res < 0.379923);

    constexpr double Ln_test_res = Ln(45);
    static_assert(Ln_test_res > 3.80666);
    static_assert(Ln_test_res < 3.80667);
}
