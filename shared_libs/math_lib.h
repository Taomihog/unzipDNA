#pragma once

#include <limits>

namespace MyMath{
    
    constexpr double Inf = std::numeric_limits<double>::infinity(); 
    constexpr double VeryLargeNumber = 1.0e20;
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double Pi = 3.1415926535897932384626433832795;
    constexpr size_t precision_ln = 30; 

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

    constexpr double Ln(double x) {
        //I wrote this! Haha:)
        
        if (x < 0.0) {
            return 1.0 / 0.0;//err;
        }
        
        if (x < 1.0) {
            return -Ln(1.0 / x);
        }
        if (x > 3.0 ) {;
            return 1.0 + Ln(x * 0.36787944117144232159552377016146);
        }
        //calculate x if 1.0 <= x <= 4.0
        double y = (x - 1) / (x + 1);
        const double ysq = y * y;
        double ypow = y; //y^+1
        double sum = ypow;
        for (size_t i = 1; i <= precision_ln; ++i) {
            ypow *= ysq;
            sum += ypow/(2.0 * i + 1);
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
        double sum2 = 0.0f;
        int i = 1;
        for (; i <= 10; ++i){
            product /= 2*i - 1;//x**1/1!, ..
            product *= x;
            sum2 += product;

            product /= 2*i; //x**2/2!, ...
            product *= x;
            sum1 += product;
        }
        product *= x;
        product /= 2*i + 1;//x**1/1!, ..
        sum2 += product;
        
        return sum1/sum2;
    }

    // constexpr double Tanh_coarse(double x) {
    //     //the polyn formula is copied from 
    //     //https://math.stackexchange.com/questions/107292/rapid-approximation-of-tanhx
    //     return 
    //     (-.67436811832e-5+(.2468149110712040+(.583691066395175e-1+.3357335044280075e-1*x)*x)*x)/
    //     (.2464845986383725+(.609347197060491e-1+(.1086202599228572+.2874707922475963e-1*x)*x)*x);
    // }

    constexpr double Tanh(double x) {
        return 1.0 / Coth(x);
    }

    constexpr double Langevin(double x) {
        if ( x > 7.9608220) {//see my comment in Coth()
            return 1.0 - 1.0 / x;
        }
        double product = x*x;//x**0/0!
        double factorial = 6;//2!
        double sum1 = 1.0/2.0 - 1.0 / 6.0;
        double sum2 = 1.0;
        for (int i = 2; i <= 10; ++i){
            sum2 += product /factorial;
            factorial *= 2.0 * i;
            sum1 += product * (1.0/factorial - 1.0 /factorial/(2.0 * i + 1.0));
            factorial *= (2.0 * i + 1.0);
            product *= x*x;
        }
        return x*sum1/sum2;
    }

    constexpr double Langevin_integ(double x) {
        // = ln(sinh(x)/x)
        
        if ( x > 7.9608220) {//see my comment in Coth()
            return (x - 7.9608220) - Ln(x / 7.9608220) + 5.1931424368760504294482139411658;
        }
        
        double sum = 0.0;
        double factor = 1.0;
        double product = 1.0;
        for (double i = 1.0; i < 20.0; ++i) {
            sum += factor;
            factor *= (x * x /(i * 2.0) / (i * 2.0 + 1.0));
        }
        //return sum;
        return Ln(sum);
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

    constexpr double Ln_test_res = Ln(45);
    static_assert(Ln_test_res > 3.80666);
    static_assert(Ln_test_res < 3.80667);
}
