// complex Lambert.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <complex>

using namespace std;

//z * exp(z)
complex<double> zexpz(complex<double> z)
{
	return z * exp(z);
}

//The derivative of z * exp(z) = exp(z) + z * exp(z)
complex<double> zexpz_d(complex<double> z)
{
	return exp(z) + z * exp(z);
}

//The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
complex<double> zexpz_dd(complex<double> z)
{
	return 2. * exp(z) + z * exp(z);
}

//Determine the initial point for the root finding
complex<double> InitPoint(complex<double> z, int k)
{
	const double pi{ 3.14159265358979323846 };
	const double e{ 2.71828182845904523536 };
	complex<double> I{0, 1};
	complex<double> two_pi_k_I{ 0., 2. * pi * k };
	complex<double> ip{ log(z) + two_pi_k_I - log(log(z) + two_pi_k_I) };// initial point coming from the general asymptotic approximation
	complex<double> p{ sqrt( 2. * (e * z + 1.) ) };// used when we are close to the branch cut around zero and when k=0,-1

	if (abs(z - (-exp(-1.))) <= 1.) //we are close to the branch cut, the initial point must be chosen carefully
	{
		if (k == 0) ip = -1. + p - 1./3. * pow(p, 2) + 11./72. * pow(p, 3);
		if (k == 1 && z.imag() < 0.) ip = -1. - p - 1./3. * pow(p, 2) - 11. / 72. * pow(p, 3);
		if (k == -1 && z.imag() > 0.) ip = -1. - p - 1./3. * pow(p, 2) - 11./72. * pow(p, 3);
	}

	if (k == 0 && abs(z - .5) <= .5) ip = (0.35173371 * (0.1237166 + 7.061302897 * z)) / (2. + 0.827184 * (1. + 2. * z));// (1,1) Pade approximant for W(0,a)

	if (k == -1 && abs(z - .5) <= .5) ip = -(((2.2591588985 +
		4.22096*I) * ((-14.073271 - 33.767687754*I) * z - (12.7127 -
			19.071643*I) * (1. + 2.*z))) / (2. - (17.23103 - 10.629721*I) * (1. + 2.*z)));// (1,1) Pade approximant for W(-1,a)

	return ip;
}

complex<double> LambertW(complex<double> z, int k = 0)
{
	//For some particular z and k W(z,k) has simple value:
	if (z == 0.) return (k == 0) ? 0. : -INFINITY;
	if (z == -exp(-1.) && (k == 0 || k == -1)) return -1.;
	if (z == exp(1.) && k == 0) return 1.;

	//Halley method begins
	complex<double> w{ InitPoint(z, k) }, wprev{ InitPoint(z, k) }; // intermediate values in the Halley method
	const unsigned int maxiter = 30; // max number of iterations. This eliminates improbable infinite loops
	unsigned int iter = 0; // iteration counter
	double prec = 1.E-30; // difference threshold between the last two iteration results (or the iter number of iterations is taken)

	do
	{
		wprev = w;
		w -= 2.*((zexpz(w) - z) * zexpz_d(w)) /
			(2.*pow(zexpz_d(w),2) - (zexpz(w) - z)*zexpz_dd(w));
		iter++;
	} while ((abs(w - wprev) > prec) && iter < maxiter);
	return w;
}

int main()
{
	complex<double> z{ 1.73, 2.118 };
	
	cout << LambertW(z, 2) << endl;

	//for (int k = -10; k <= 10; k++) cout << "W = " << LambertW(z, k) << endl;

    return 0;
}

