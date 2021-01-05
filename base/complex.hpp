//
// Created by rokeabbey on 2020/12/21.
//

#ifndef FFT_IMPL_COMPLEX_HPP
#define FFT_IMPL_COMPLEX_HPP

#include <string>


#define COMPLEX_ADD_ASSIGN_MULTI(o, a, b) { \
    (o).re += (a).re * (b).re - (a).im * (b).im; \
    (o).im += (a).im * (b).re + (a).re * (b).im; \
}

/**
 * 请注意 o与a不能相同
 */
#define COMPLEX_MULTI(o, a, b) { \
    (o).re = (a).re * (b).re - (a).im * (b).im; \
    (o).im = (a).im * (b).re + (a).re * (b).im; \
}

#define COMPLEX_MULTI_ASSIGN(o, a, tmp) { \
    tmp = o.re;                           \
    o.re = o.re * a.re - o.im * a.im;     \
    o.im = tmp * a.im + o.im * a.re;      \
}
#define COMPLEX_ADD(o, a, b) { \
    (o).re = (a).re + (b).re;  \
    (o).im = (a).re + (b).re;  \
}
namespace rk {
    template<typename Tp>
    class Complex {
    public:
        Tp re, im;

        Complex() : re(0), im(0) {}

        Complex(Tp re, Tp im) : re(re), im(im) {}

        inline
        Complex<Tp> operator*(Complex<Tp> &a) {
            Complex<Tp> result;
            result.re = re * a.re - im * a.im;
            result.im = re * a.im + im * a.re;
            return result;
        }

        inline
        Complex<Tp> operator*(Tp &n) {
            Complex<Tp> c;
            c.re = n * re;
            c.im = n * im;
            return c;
        }

        inline
        void operator*=(const Complex<Tp> &a) {
            Tp _re = re;
            re = re * a.re - im * a.im;
            im = _re * a.im + im * a.re;
        }

        inline
        void operator*=(const Tp &n) {
            re *= n;
            im *= n;
        }

        inline
        void operator+=(const Complex<Tp> &a) {
            re += a.re;
            im += a.im;
        }

        inline
        Complex<Tp> operator+(Complex<Tp> &a) {
            return Complex<Tp>(re + a.re, im + a.im);
        }

        inline
        Complex<Tp> conj(){
            return Complex<Tp>(re, -im);
        }


        std::string toString() {
            return std::to_string(re) + " + " + std::to_string(im) + "i";
        }


    };
}

#endif //FFT_IMPL_COMPLEX_HPP
