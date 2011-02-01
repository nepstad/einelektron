#ifndef SPHERICALVELOCITYXBODY_H
#define SPHERICALVELOCITYXBODY_H

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include "sphericalbase.h"




class velocityHelperX
{
public:
	static double sphericalvelocityBodyX(int l, int m, int lp, int mp)	
	{
	

		//SphericalBasis::ClebschGordan cg;

		double coupling = 0;
		double I1_1 = 0;
		double I2_1 = 0; 
		double I3_1 = 0;
		double dlta1 = 0;
		double dlta2 = 0;
		double norms = 0;
		double temp = 0;
		double E_lm;
		double F_lm;
		double G_lm;
		double H_lm;
		double I_lm;
		double J_lm;
		double K_lm;
		double dlta_m;
		double K2_term1;
		double K2_term2;
		double K1_term1;
		double K1_term2;
		double K1_term3;

		double J1 = 0;
		double J2 = 0;






		I1_1 += I1integral(lp,mp,l,m);

		dlta1 = Mykronecker(mp - 1, m);
		dlta2 = Mykronecker(mp + 1, m);	

		/*
		* Integral I2 in r1.
		*/
		norms = LegendreNorm(l,m) * LegendreNorm(lp,mp);

		temp = M_PI * m * norms * (dlta1 - dlta2); 
		temp *= K1(lp,std::abs(mp),l,std::abs(m));
		I2_1 += temp;


		/*
		 * Integral J1 in r1.
		 */
		F_lm = F(l,m);
		G_lm = G(l,m);
		H_lm = H(l,m);
		I_lm = I(l,m);

		K1_term1 = K1(lp,std::abs(mp),l-2,std::abs(m));
		K1_term2 = K1(lp,std::abs(mp),l,std::abs(m));
		K1_term3 = K1(lp,std::abs(mp),l+2,std::abs(m));

		temp = LegendreNorm(l-2,m) * F_lm * K1_term1;
		temp += LegendreNorm(l,m) * (G_lm + H_lm) * K1_term2;
		temp += LegendreNorm(l+2,m) * I_lm * K1_term3;
		temp *= M_PI * (dlta1 + dlta2) * LegendreNorm(lp,mp);


		I3_1 += std::abs(m) * temp;
		J1 += std::abs(m) * temp;

		/*
		 * Integral J2 in r1.
		 */
		J_lm = J(l,m);
		K_lm = K(l,m);
		E_lm = E(l,m);
	
		dlta_m = delta_m(m);


		K2_term1 = K2(l-1,std::abs(m+dlta_m),lp,std::abs(mp));
		K2_term2 = K2(l+1,std::abs(m+dlta_m),lp,std::abs(mp));

		temp = LegendreNorm(l-1,m+dlta_m) * J_lm * K2_term1;
		temp += LegendreNorm(l+1,m+dlta_m) * K_lm * K2_term2;
		temp *= M_PI * LegendreNorm(lp,mp) * (dlta1 + dlta2);



		I3_1 -= dlta_m * E_lm * temp;  
		J2 -= dlta_m * E_lm * temp;

		
		coupling += (-I1_1 - I2_1 + I3_1);
		
		//cout << l << "\t"  << m << "\t" << lp << "\t"  << mp << "\t" << I1_1 << "\t" << I2_1 << "\t" << I3_1 << "\t" << coupling << "\t Js:" << J1 << "\t" << J2 << endl;

		return coupling;
	}
	
	static double LegendreNorm(double l, double m)
	{
		/*
		 * Y_{l,m} = LegendreNorm(l,m) * P_l^m
		 */

		if (std::abs(m)<=l && l>=0)
		{
			double norm;
			norm = std::pow(-1.,.5 * (m - std::abs(m))) * sqrt((2. * l + 1.) / (4. * M_PI));
			norm *= std::sqrt(exp(gsl_sf_lnfact(l - std::abs(m)) - gsl_sf_lnfact(l + std::abs(m))));
			return norm;
		}
		else
		{
			return 0.;
		}
	}
	

	
	static int Mykronecker(int a, int b)
	{
		return a == b ? 1 : 0;
	}

	static double E(double l, double m)
	{
		double j = (l - std::abs(m)) * (l + std::abs(m) + 1.); 
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( (l - std::abs(m)) * (l + std::abs(m) + 1.) );
		}
	}

	static double F(double l, double m)
	{
		double j = ((l+m)*(l-m)*(l+m-1.)*(l-m-1.)) / ((2.*l+1.)*std::pow(2.*l-1.,2)*(2.*l-3.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m)*(l-m)*(l+m-1.)*(l-m-1.)) / ((2.*l+1.)*std::pow(2.*l-1.,2)*(2.*l-3.)) );
		}
	}

	static double G(double l, double m)
	{
		{
			return ((l+m)*(l-m)) / ((2.*l+1.)*(2.*l-1.));
		}
	}

	static double H(double l, double m)
	{
		{
			return ((l+m+1.)*(l-m+1.)) / ((2.*l+1.)*(2.*l+3.)) ;
		}
	}

	static double I(double l, double m)
	{
		double j = ((l+m+1.)*(l-m+1.)*(l+m+2.)*(l-m+2.)) / ((2.*l+1.)*std::pow(2.*l+3.,2)*(2.*l+5.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+1.)*(l-m+1.)*(l+m+2.)*(l-m+2.)) / ((2.*l+1.)*std::pow(2.*l+3.,2)*(2.*l+5.)) );
		}
	}

	static double delta_m(double m)
	{
		return m >= 0. ? 1. : -1.;
	}


	static double J(double l, double m)
	{
		double j = ((l+m+delta_m(m))*(l-m-delta_m(m))) / ((2.*l+1.)*(2.*l-1.));
		if (j < 0.)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+delta_m(m))*(l-m-delta_m(m))) / ((2.*l+1.)*(2.*l-1.)) );
		}
	}



	static double K(double l, double m)
	{
		double j = ((l+m+delta_m(m)+1.)*(l-m-delta_m(m)+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < 0)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( ((l+m+delta_m(m)+1.)*(l-m-delta_m(m)+1.)) / ((2.*l+1.)*(2.*l+3.)) );
		}
	}



	static double K1(int l, int m, int p, int q)
	{
		int c1 = l-p;
		int c2 = m-q;
		int nu = std::min(l,p);
		int mu = std::min(m,q);
		
		if ((l >= std::abs(m)) && (p >= std::abs(q)))
		{
			if ((l+p) % 2==1)
			{
				if (((c1 < 0) && (c2 < 0)) || ((c1 > 0) && (c2 > 0)))
				{
					return -2.0 * exp(gsl_sf_lnfact(nu+mu) - gsl_sf_lnfact(nu-mu));
				}
				else
				{
					return 0.0;
				}
			}
			else
			{
				return 0.;
			}
		}
		else
		{
			return 0.0;
		}
	}


	static double K2(int l, int m, int p, int q)
	{
		if ((l >= std::abs(m)) && (p >= std::abs(q)))
		{
			if ((l+p-m-q) % 2 == 0)
			{
				double outerSum = 0.0;
				int outerMax = std::floor((l-m)/2.);
				int innerMax = std::floor((p-q)/2.);

				for (int i=0; i<=outerMax; i++)
				{
					double innerSum = .0;
					for (int j=0; j<=innerMax; j++)
					{
						innerSum += Cconstant(p,q,j) * exp(gsl_sf_lngamma(.5 * (l+p-m-q -2.*(i+j)+1.)) + gsl_sf_lngamma(.5 * (m+q+2.*(i+j+1.))) - gsl_sf_lngamma(.5*(l+p+3.)));
					}
					outerSum += innerSum * Cconstant(l,m,i); 
				}
				return outerSum;
			
			}
			else
			{
		 		return 0.0;	
			}
		}
		else
		{	
			return 0.0;
		}
	}

	static double Cconstant(double alpha, double beta, double gamma)
	{
		double C;
		C = std::pow(-1.,gamma) * std::pow(2.,-(beta + 2 * gamma)); 
		C *= exp(gsl_sf_lnfact(alpha + beta) - gsl_sf_lnfact(beta + gamma)-gsl_sf_lnfact(gamma) - gsl_sf_lnfact(alpha - beta - 2. * gamma));
		return C;
	}

	static double MyCoefficient(double lp, double l)
	{
		return std::sqrt(.5 * (2. * l + 1.) / (2. * lp + 1.) );
	}

			/*
			 * Integral I1 in r1.
			 */

	static double I1integral(double lp, double mp, double l, double m)
	{
		SphericalBasis::ClebschGordan cg;
		return MyCoefficient(lp,l) * cg(l,1,0,0,lp,0) * (cg(l,1,m,-1,lp,mp) - cg(l,1,m,1,lp,mp));
	}		

};



#endif
