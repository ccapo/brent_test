#include <iostream>
#include <iomanip>
#include <cmath>

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

template<class T, class F> std::pair<T, T> brent(F f, T xmin, T xmax, const T tol) {
	T a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	T e = 0.0;
	std::pair<T, T> res;

	a = (xmin < xmax ? xmin : xmax);
	b = (xmin > xmax ? xmax : xmin);
	x = w = v = xmax;
	fw = fv = fx = (f)(x);
	
	for(int iter = 1; iter <= ITMAX; iter++) {
		xm = 0.5*(a + b);
		tol1 = tol*fabs(x) + ZEPS;
		tol2 = 2.0*tol1;
		if(fabs(x - xm) <= (tol2 - 0.5*(b - a))) {
			res.first = x;
			res.second = fx;
			return res;
		}

		if (fabs(e) > tol1) {
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0) {
        p = -p;
      }
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a - x) || p >= q*(b - x)) {
        if (x >= xm) {
          e = a - x;
        } else {
          e = b - x;
        }
				d = CGOLD*e;
			} else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2) {
					d = SIGN(tol1, xm - x);
        }
			}
		} else {
      if (x >= xm) {
        e = a - x;
      } else {
        e = b - x;
      }
			d = CGOLD*e;
		}

    if (fabs(d) >= tol1) {
      u = x + d;
    } else {
      u = x + SIGN(tol1, d);
    }

		fu = (f)(u);

		if (fu <= fx) {
			if (u >= x) {
        a = x;
      } else {
        b = x;
      }
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) {
        a = u;
      } else {
        b = u;
      }
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}

  std::cout << "Exceeded maximum number of iterations\n";
	res.first = x;
	res.second = fx;
	return res;
}