#include <core/common.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>

/*
 * Sets the radial Coulomb wave F_l(k*r, eta), with eta = Z/k into data for all radial
 * grid points specified by r
 */
void SetRadialCoulombWave(int Z, int l, double k, blitz::Array<double, 1> r, blitz::Array<double, 1> data)
{
	double eta = Z / k;

	for (int i=0; i<r.size(); i++)
	{	
		double x = k * r(i);
		gsl_sf_result F, Fp, G, Gp;
		double exp_F, exp_G;
		int error = gsl_sf_coulomb_wave_FG_e(eta, x, (double)l, 0., &F, &Fp, &G, &Gp, &exp_F, &exp_G);
		if (error == GSL_EOVRFLW)
		{
			cout << "WARNING: Overflow in SetCoulombWave(" << Z << ", " << l << ", " << k << ", r=" << r(i) << ");" << endl;
			cout << "         exp_F = " << exp_F << ", exp_G = " << exp_G << endl;
		}

		data(i) = F.val;
	}
}

double GetCoulombNormalization(double Z, int l, double k)
{
	double eta = Z / k;
	gsl_sf_result C;
	int error = gsl_sf_coulomb_CL_e(l, eta, &C);
	return C.val;
}

/* 
 * Gets the Coulomb phase sigma_l = arg(gamma(l + 1 + i*eta))
 */
double GetCoulombPhase(int l, double eta)
{
	gsl_sf_result absval, argval;
	if (gsl_sf_lngamma_complex_e(1.0+l, eta, &absval, &argval) == GSL_ELOSS)
	{
		cout << "Overflow error in gsl_sf_lngamma_complex_e, l=" << l << ", eta=" << eta << endl;
	}
	return argval.val;
}

