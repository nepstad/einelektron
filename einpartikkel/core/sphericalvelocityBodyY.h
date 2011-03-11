#ifndef SPHERICALVELOCITYYBODY_H
#define SPHERICALVELOCITYYBODY_H

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include "sphericalbase.h"




class velocityHelperY
{
public:
	static double sphericalvelocityBodyY(int l, int m, int lp, int mp)	
	{

		double eps = std::pow(10.,-15);

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
		//norms = LegendreNormDouble(l,m,lp,mp);

		//temp = M_PI * m * norms * (dlta1 + dlta2); 
		//temp *= K1(lp,std::abs(mp),l,std::abs(m));
		//I2_1 += temp;



		/*
		 * I2 stable
		 */

		bool lnNorm_ok = checklnLegendreNorm(l,m,lp,mp);
		bool lnK1_ok = checklnK1(lp,std::abs(mp),l,std::abs(m));


		if ((lnNorm_ok == true) && (lnK1_ok == true))
		{
			double lnNorm = lnLegendreNorm(l,m,lp,mp);
			double lnK1Int = lnK1(lp,std::abs(mp),l,std::abs(m));

			I2_1 += exp(lnNorm + lnK1Int);
			I2_1 *= (dlta1 + dlta2); 
			I2_1 *= 0.5 * m * pow(-1,0.5 * (m + mp + std::abs(m) + std::abs(mp)));
		}

		//cout << "samme " << I2_1 << " " << temp  << " \t " << l << " " << m << " " << lp << " " << mp << endl;


		/*
		 * Integral J1 in r1.
		 */
	/*
		F_lm = F(l,m,eps);
		G_lm = G(l,m,eps);
		H_lm = H(l,m,eps);
		I_lm = I(l,m,eps);

		K1_term1 = K1(lp,std::abs(mp),l-2,std::abs(m));
		K1_term2 = K1(lp,std::abs(mp),l,std::abs(m));
		K1_term3 = K1(lp,std::abs(mp),l+2,std::abs(m));

		temp = LegendreNormDouble(l-2,m,lp,mp) * F_lm * K1_term1;
		temp += LegendreNormDouble(l,m,lp,mp) * (G_lm + H_lm) * K1_term2;
		temp += LegendreNormDouble(l+2,m,lp,mp) * I_lm * K1_term3;
		temp *= M_PI * (-dlta1 + dlta2); 


		I3_1 += std::abs(m) * temp;

	*/
		/*
		 * J1 stable
		 */
		bool lnNorm1_ok = checklnLegendreNorm(l-2,m,lp,mp);
		bool lnNorm2_ok = checklnLegendreNorm(l,m,lp,mp);
		bool lnNorm3_ok = checklnLegendreNorm(l+2,m,lp,mp);

		bool lnF_ok = checklnF(l,m,eps);
		bool lnG_ok = checklnG(l,m,eps);
		bool lnH_ok = checklnH(l,m,eps);
		bool lnI_ok = checklnI(l,m,eps);

		bool lnK1Int_1_ok = checklnK1(lp,std::abs(mp),l-2,std::abs(m));
		bool lnK1Int_2_ok = checklnK1(lp,std::abs(mp),l,std::abs(m));
		bool lnK1Int_3_ok = checklnK1(lp,std::abs(mp),l+2,std::abs(m));


		if ((lnNorm1_ok == true) && (lnK1Int_1_ok))
		{
			if (lnF_ok == true)
			{
				double lnNorm1 = lnLegendreNorm(l-2,m,lp,mp);
				double lnK1Int_1 = lnK1(lp,std::abs(mp),l-2,std::abs(m));
				double lnFconst = lnF(l,m);

				J1 += exp(lnNorm1 + lnK1Int_1 + lnFconst);
			}
		}

		if ((lnNorm2_ok == true) && (lnK1Int_2_ok == true))
		{
			double lnNorm2 = lnLegendreNorm(l,m,lp,mp);
			double lnK1Int_2 = lnK1(lp,std::abs(mp),l,std::abs(m));

			if (lnG_ok == true)
			{
				double lnGconst = lnG(l,m);

				J1 += exp(lnNorm2 + lnK1Int_2 + lnGconst);
			}

			if (lnH_ok == true)
			{
				double lnHconst = lnH(l,m);

				J1 += exp(lnNorm2 + lnK1Int_2 + lnHconst);
			}
		}

		if ((lnNorm3_ok == true) && (lnK1Int_3_ok == true))
		{
			if (lnI_ok == true)
			{
				double lnNorm3 = lnLegendreNorm(l+2,m,lp,mp);
				double lnK1Int_3 = lnK1(lp,std::abs(mp),l+2,std::abs(m));
				double lnIconst = lnI(l,m);

				J1 += exp(lnNorm3 + lnK1Int_3 + lnIconst);
			}
		}

		
		J1 *= (-dlta1 + dlta2);
		J1 *= std::abs(m) * 0.5 * pow(-1,0.5 * (m + mp + std::abs(m) + std::abs(mp)));


		//cout << "samme " << J1 << " "  << std::abs(m) * temp << " \t " << l << " " << m << " " << lp << " " << mp << endl;



		/*
		 * Integral J2 in r1.
		 */
	 
		J_lm = J(l,m,eps);
		K_lm = K(l,m,eps);
		E_lm = E(l,m);
	
	
		dlta_m = delta_m(m);


		K2_term1 = K2(l-1,std::abs(m+dlta_m),lp,std::abs(mp));
		K2_term2 = K2(l+1,std::abs(m+dlta_m),lp,std::abs(mp));

		temp = LegendreNormDouble(l-1,m+dlta_m,lp,mp) * J_lm  * K2_term1;
		temp += LegendreNormDouble(l+1,m+dlta_m,lp,mp) * K_lm * K2_term2;
		temp *= M_PI * (-dlta1 + dlta2);



		double JJ2 = dlta_m * E_lm * temp;  
		
	
		/*
		 * J2 stable
		 */

	
		dlta_m = delta_m(m);

		bool lnNorm4_ok = checklnLegendreNorm(l-1,m+dlta_m,lp,mp);
		bool lnNorm5_ok = checklnLegendreNorm(l+1,m+dlta_m,lp,mp);

		bool lnJ_ok = checklnJ(l,m,eps);
		bool lnK_ok = checklnK(l,m,eps);
		bool lnE_ok = checklnE(l,m,eps);


		if (lnE_ok == true)
		{
			//double lnEconst = lnE(l,m);
			vector<double> lnEconst = lnEvect(l,m);

			if ((lnNorm4_ok == true) && (lnJ_ok == true))
			{
				//double lnNorm4 = lnLegendreNorm(l-1,m+dlta_m,lp,mp);
				//cout << "tt" <<  lnJ(l,m,eps)+lnLegendreNorm(l-1,m+dlta_m,lp,mp)<< endl;

				vector<double> lnNorm4 = lnLegendreNormVect(l-1,m+dlta_m,lp,mp);
				vector<double> lnJconst = lnJvect(l,m,eps);

				//double p = 0;
				//for (int idx = 0; idx < lnJconst.size(); idx++)
				//{
				//	p += lnJconst[idx];
				//}	

				//cout << "ss" << p << endl;

				vector<double> t;
				t.reserve(lnNorm4.size() + lnJconst.size());
				t.insert(t.end(), lnNorm4.begin(), lnNorm4.end()); 
				t.insert(t.end(), lnJconst.begin(), lnJconst.end());

				vector<double> s;
				s.reserve(t.size() + lnEconst.size());
				s.insert(s.end(), t.begin(), t.end()); 
				s.insert(s.end(), lnEconst.begin(), lnEconst.end());

				//J2 += lnSumK2(l-1,std::abs(m+dlta_m),lp,std::abs(mp),lnNorm4 + lnJconst + lnEconst);
				J2 += lnSumK2(l-1,std::abs(m+dlta_m),lp,std::abs(mp),s);
			}
		 
			if ((lnNorm5_ok == true) && (lnK_ok == true))
			{
				//cout << "tt" << lnK(l,m) +  lnLegendreNorm(l+1,m+dlta_m,lp,mp) + lnE(l,m) << endl;

				vector<double> lnNorm5 = lnLegendreNormVect(l+1,m+dlta_m,lp,mp);
				vector<double> lnKconst = lnKvect(l,m);

				//double s = 0;
				//for (int idx = 0; idx < lnKconst.size(); idx++)
				//{
				//	s += lnKconst[idx];
				//}	

				//cout << "ss" << s << endl;

				vector<double> u;
				u.reserve(lnNorm5.size() + lnKconst.size());
				u.insert(u.end(), lnNorm5.begin(), lnNorm5.end()); 
				u.insert(u.end(), lnKconst.begin(), lnKconst.end());

				vector<double> v;
				v.reserve(u.size() + lnEconst.size());
				v.insert(v.end(), u.begin(), u.end()); 
				v.insert(v.end(), lnEconst.begin(), lnEconst.end());

				double s = 0;
				for (int idx = 0; idx < v.size(); idx++)
				{
					s += v[idx];
				}	

				//cout << "--" << s << endl;

				//J2 += lnSumK2(l+1,std::abs(m+dlta_m),lp,std::abs(mp),lnNorm5 + lnKconst + lnEconst);
				J2 += lnSumK2(l+1,std::abs(m+dlta_m),lp,std::abs(mp),v);
			}
		}

		J2 *= dlta_m * 0.25 * (-dlta1 + dlta2) * std::pow(-1.,0.5 * (mp + std::abs(mp) + m + dlta_m + std::abs(m + dlta_m)));

	

		
		//coupling += (-I1_1 + I2_1 + I3_1);
		//cout << "coupling " << coupling << " " <<-I1_1 + I2_1 + J1 + J2 << endl; 

		coupling += (-I1_1 + I2_1 + J1 + J2);
		//cout << "fff" << " " << J2 << " " << JJ2 << " " << " \t " << l << " " << m << " " << lp << " " << mp << endl;

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
			norm = std::pow(-1.,.5 * (m + std::abs(m))) * sqrt((2. * l + 1.) / (4. * M_PI));
			norm *= std::sqrt(exp(gsl_sf_lnfact(l - std::abs(m)) - gsl_sf_lnfact(l + std::abs(m))));
			return norm;
		}
		else
		{
			return 0.;
		}
	}
	
	static double LegendreNormDouble(double l, double m, double lp, double mp)
	{

		if ( (std::abs(m)<=l && l>=0) && (std::abs(mp)<=lp && lp>=0))
		{
			double tmp = gsl_sf_lnfact(l - std::abs(m)) + gsl_sf_lnfact(lp - std::abs(mp)); 
			tmp -= gsl_sf_lnfact(lp + std::abs(mp)) + gsl_sf_lnfact(l + std::abs(m));
			tmp += log(2*l + 1) + log(2*lp+1);

			double norm = std::pow(-1.,.5 * (m + mp + std::abs(m) + std::abs(mp))) / (4 * M_PI) * std::sqrt(exp(tmp));

			return norm;
		}
		else
		{
			return 0.;
		}
	}

	static double lnLegendreNorm(double l, double m, double lp, double mp)
	{
		double tmp = log(2 * l + 1.) + log(2 * lp + 1.);
		tmp += gsl_sf_lnfact(l - std::abs(m)) + gsl_sf_lnfact(lp - std::abs(mp));
		tmp -= gsl_sf_lnfact(l + std::abs(m)) + gsl_sf_lnfact(lp + std::abs(mp));
		tmp *= 0.5;
		return tmp;
	}

	static vector<double> lnLegendreNormVect(double l, double m, double lp, double mp)
	{
		vector<double> v;

		v.push_back(0.5 * log(2 * l + 1.));
		v.push_back(0.5 * log(2 * lp + 1.));
		v.push_back(0.5 * gsl_sf_lnfact(l - std::abs(m)));
		v.push_back(0.5 * gsl_sf_lnfact(lp - std::abs(mp)));
		v.push_back((-0.5) * gsl_sf_lnfact(l + std::abs(m)));
		v.push_back((-0.5) * gsl_sf_lnfact(lp + std::abs(mp)));

		return v;
	}

	static bool checklnLegendreNorm(double l, double m, double lp, double mp)
	{
		if ( (std::abs(m)<=l && l>=0) && (std::abs(mp)<=lp && lp>=0))
		{
			return true;
		}
		else
		{
			return false;
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

	static double lnE(double l, double m)
	{
		return 0.5 * (log(l - std::abs(m)) + log(l + std::abs(m) + 1.));
	}



	static vector<double> lnEvect(double l, double m)
	{
		vector<double> v;

		v.push_back(0.5 * log(l - std::abs(m)));
		v.push_back(0.5 * log(l + std::abs(m) + 1.));
		return v;
	}


	static bool checklnE(double l, double m, double eps)
	{
		double j = (l - std::abs(m)) * (l + std::abs(m) + 1.); 
		if (j < 0.)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	static double F(double l, double m, double eps)
	{
		double j = ((l+m)*(l-m)*(l+m-1.)*(l-m-1.)) / ((2.*l+1.)*std::pow(2.*l-1.,2)*(2.*l-3.));
		if (j < eps)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( exp(log(l+m) + log(l-m) + log(l+m-1) + log(l-m-1) - log(2.*l+1.) - 2 * log(2.*l-1.) - log(2.*l-3.)) );
		}
	}


	static double lnF(double l, double m)
	{
		
		return 0.5 * (log(l+m) + log(l-m) + log(l+m-1) + log(l-m-1) - log(2.*l+1.) - 2 * log(2.*l-1.) - log(2.*l-3.));
	}


	static bool checklnF(double l, double m, double eps)
	{
		double j = ((l+m)*(l-m)*(l+m-1.)*(l-m-1.)) / ((2.*l+1.)*std::pow(2.*l-1.,2)*(2.*l-3.));
		if (j < eps)
		{
			return false;
		}
		else
		{
			return true;
		}
	}


	static double G(double l, double m, double eps)
	{
		double j = ((l+m)*(l-m)) / ((2.*l+1.)*(2.*l-1.));
		if (j < eps)
		{
			return 0.;
		}
		else
		{
			return exp(log(l+m) + log(l-m) - log(2*l+1) - log(2*l-1));
		}
	}


	static double lnG(double l, double m)
	{
		return log(l+m) + log(l-m) - log(2*l+1) - log(2*l-1);
	}


	static bool checklnG(double l, double m, double eps)
	{
		double j = ((l+m)*(l-m)) / ((2.*l+1.)*(2.*l-1.));
		if (j < eps)
		{
			return false;
		}
		else
		{
			return true;
		}
	}


	static double H(double l, double m, double eps)
	{
		double j = ((l+m+1.)*(l-m+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < eps)
		{
			return 0.;
		}
		else
		{
			return exp( log(l+m+1) + log(l-m+1) - log(2*l+1) - log(2*l+3)) ;
		}
	}


	static double lnH(double l, double m)
	{
		return log(l+m+1) + log(l-m+1) - log(2*l+1) - log(2*l+3);
	}


	static bool checklnH(double l, double m, double eps)
	{
		double j = ((l+m+1.)*(l-m+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < eps)
		{
			return false;
		}
		else
		{
			return true;
		}
	}


	static double I(double l, double m, double eps)
	{
		double j = ((l+m+1.)*(l-m+1.)*(l+m+2.)*(l-m+2.)) / ((2.*l+1.)*std::pow(2.*l+3.,2)*(2.*l+5.));
		if (j < eps)
		{
			return 0.;
		}
		else
		{
			return std::sqrt( exp( log(l+m+1) + log(l-m+1) + log(l+m+2) + log(l-m+2) - log(2.*l+1.) - 2 * log(2.*l+3) - log(2.*l+5)) );
		}
	}


	static double lnI(double l, double m)
	{

		return 0.5 * (  log(l+m+1) + log(l-m+1) + log(l+m+2) + log(l-m+2) - log(2.*l+1.) - 2 * log(2.*l+3) - log(2.*l+5));
	}


	static bool checklnI(double l, double m, double eps)
	{
		double j = ((l+m+1.)*(l-m+1.)*(l+m+2.)*(l-m+2.)) / ((2.*l+1.)*std::pow(2.*l+3.,2)*(2.*l+5.));
		if (j < eps)
		{
			return false;
		}
		else
		{
			return true; 
		}
	}


	static double delta_m(double m)
	{
		return m >= 0. ? 1. : -1.;
	}


	static double J(double l, double m, double eps)
	{
		double j = ((l+m+delta_m(m))*(l-m-delta_m(m))) / ((2.*l+1.)*(2.*l-1.));
		if (std::abs(l) < eps)
		{
			return 1.;
		}
		else
		{
			if (j < eps)
			{
				return 0.;
			}
			else
			{
				return std::sqrt( exp(log(l+m+delta_m(m)) + log(l-m-delta_m(m)) - log(2.*l+1.) - log(2.*l-1.)) );
			}
		}
	}

	static double lnJ(double l, double m, double eps)
	{
		if (std::abs(l) < eps)
		{
			return 0.;
		}
		else
		{
			return 0.5 * (log(l+m+delta_m(m)) + log(l-m-delta_m(m)) - log(2.*l+1.) - log(2.*l-1.) );
		}
	}

	static vector<double> lnJvect(double l, double m, double eps)
	{
		vector<double> v;

		if (std::abs(l) < eps)
		{
			v.push_back(0.);
			return v;
		}
		else
		{
			v.push_back(0.5 * log(l+m+delta_m(m)));
			v.push_back(0.5 * log(l-m-delta_m(m)));
			v.push_back((-0.5) * log(2.*l+1.));
			v.push_back((-0.5) * log(2.*l-1.));

			return v;
		}
	}



	static bool checklnJ(double l, double m, double eps)
	{
		double j = ((l+m+delta_m(m))*(l-m-delta_m(m))) / ((2.*l+1.)*(2.*l-1.));
		if (j < eps)
		{	
			return false;
		}
		else
		{
			return true;
		}
	}

	static double K(double l, double m, double eps)
	{
		double j = ((l+m+delta_m(m)+1.)*(l-m-delta_m(m)+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < eps)
		{
			return 0.;
		}
		else
		{
			return std::sqrt(exp(log(l+m+delta_m(m)+1.) + log(l-m-delta_m(m)+1.)  - log(2.*l+1.) - log(2.*l+3.)));
		}
	}


	static double lnK(double l, double m)
	{
		return 0.5 * (log(l+m+delta_m(m)+1.) + log(l-m-delta_m(m)+1.)  - log(2.*l+1.) - log(2.*l+3.));
	}

	static vector<double> lnKvect(double l, double m)
	{
		vector<double> v;

		v.push_back(0.5 * log(l+m+delta_m(m)+1.));
		v.push_back(0.5 * log(l-m-delta_m(m)+1.));
		v.push_back((-0.5) * log(2.*l+1.));
		v.push_back((-0.5) * log(2.*l+3.));

		return v;
	}

	static bool checklnK(double l, double m, double eps)
	{
		double j = ((l+m+delta_m(m)+1.)*(l-m-delta_m(m)+1.)) / ((2.*l+1.)*(2.*l+3.));
		if (j < eps)
		{
			return false;
		}
		else
		{
			return true;
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
					return 2.0 * exp(gsl_sf_lnfact(nu+mu) - gsl_sf_lnfact(nu-mu));
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

	static double lnK1(int l, int m, int p, int q)
	{
		int nu = std::min(l,p);
		int mu = std::min(m,q);
	
		return gsl_sf_lnfact(nu+mu) - gsl_sf_lnfact(nu-mu);
	}


	static bool checklnK1(int l, int m, int p, int q)
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
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}



	static double K2(int l, int m, int p, int q)
	{
		if ((l >= std::abs(m)) && (p >= std::abs(q)))
		{
			if ((l+p-m-q) % 2 == 0)
			{
				double outerSum = 0.0;
				double outerSumOld = 0.0;
				int outerMax = std::floor((l-m)/2.);
				int innerMax = std::floor((p-q)/2.);

				vector<double> v;

				for (int i=0; i<=outerMax; i++)
				{

					for (int j=0; j<=innerMax; j++)
					{

					
						double gmArg = gsl_sf_lngamma(.5 * (l+p-m-q -2.*(i+j)+1.)) + gsl_sf_lngamma(.5 * (m+q+2.*(i+j+1.))) - gsl_sf_lngamma(.5*(l+p+3.));

						gmArg += gsl_sf_lngamma(p+q+1)-gsl_sf_lngamma(q+j+1)-gsl_sf_lngamma(j+1)-gsl_sf_lngamma(p-q-2*j+1);
						gmArg += gsl_sf_lngamma(l+m+1)-gsl_sf_lngamma(m+i+1)-gsl_sf_lngamma(i+1)-gsl_sf_lngamma(l-m-2*i+1);
						gmArg -= log(2) *( m + q + 2 * (i + j));

						double tmp = std::pow(-1,i+j) * exp(gmArg);
						
						v.push_back(tmp);
	

					}
					


				}


				sort(v.begin(),v.end(),magnitude());
				
				for (int idx = 0; idx < v.size(); idx++)
				{
					outerSum += v[idx];
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







	//static double lnSumK2(int l, int m, int p, int q,double lnsum)
	static double lnSumK2(int l, int m, int p, int q,vector<double> lnsum)
	{
		if ((l >= std::abs(m)) && (p >= std::abs(q)))
		{
			if ((l+p-m-q) % 2 == 0)
			{
				double outerSum = 0.0;
				double outerSumOld = 0.0;
				int outerMax = std::floor((l-m)/2.);
				int innerMax = std::floor((p-q)/2.);

				vector<double> v;


				for (int i=0; i<=outerMax; i++)
				{

					for (int j=0; j<=innerMax; j++)
					{

						vector<double> w;
						w.push_back(gsl_sf_lngamma(.5 * (l+p-m-q -2.*(i+j)+1.)));
						w.push_back(gsl_sf_lngamma(.5 * (m+q+2.*(i+j+1.))));
						w.push_back((-1.) * gsl_sf_lngamma(.5*(l+p+3.)));
						w.push_back(gsl_sf_lngamma(p+q+1));
						w.push_back((-1.) * gsl_sf_lngamma(q+j+1));
						w.push_back((-1.) * gsl_sf_lngamma(j+1));
						w.push_back((-1.) * gsl_sf_lngamma(p-q-2*j+1));
						w.push_back(gsl_sf_lngamma(l+m+1));
						w.push_back((-1.) * gsl_sf_lngamma(m+i+1));
						w.push_back((-1.) * gsl_sf_lngamma(i+1));
						w.push_back((-1.) * gsl_sf_lngamma(l-m-2*i+1));
						w.push_back((-1.) * log(2) *( m + q + 2 * (i + j)));


						vector<double> u;
						u.reserve(w.size() + lnsum.size());
						u.insert(u.end(), w.begin(), w.end()); 
						u.insert(u.end(), lnsum.begin(), lnsum.end());

						sort(u.begin(),u.end(),magnitude());

						double tmp = 0.;
	
						for (int idx = 0; idx < u.size(); idx++)
						{
							tmp += u[idx];
						}		


						double gmArg = gsl_sf_lngamma(.5 * (l+p-m-q -2.*(i+j)+1.)) + gsl_sf_lngamma(.5 * (m+q+2.*(i+j+1.))) - gsl_sf_lngamma(.5*(l+p+3.));


						gmArg += gsl_sf_lngamma(p+q+1)-gsl_sf_lngamma(q+j+1)-gsl_sf_lngamma(j+1)-gsl_sf_lngamma(p-q-2*j+1);
						gmArg += gsl_sf_lngamma(l+m+1)-gsl_sf_lngamma(m+i+1)-gsl_sf_lngamma(i+1)-gsl_sf_lngamma(l-m-2*i+1);
						gmArg -= log(2) *( m + q + 2 * (i + j));

						//cout << "SS " << gmArg << endl;

						double tamp = 0;
						for (int idx = 0; idx < lnsum.size(); idx++)
						{
							tamp += lnsum[idx];
						}

						double tame = 0;
						for (int idx = 0; idx < w.size(); idx++)
						{
							tame += w[idx];
						}

						//cout << "AA" << tamp << endl; 
						//cout << "BA" << tame << endl; 

						//double tmp = std::pow(-1,i+j) * exp(gmArg + lnsum);

						double tmp1 = exp(tmp); 
						tmp1 *= std::pow(-1,i+j);


						v.push_back(tmp1);

					}

				}


				sort(v.begin(),v.end(),magnitude());
				
				for (int idx = 0; idx < v.size(); idx++)
				{
					outerSum += v[idx];
				}		
				
				//cout << "outer" << outerSum << endl;
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

	static double I1integral(double lp,double mp,double l,double m)
	{
		SphericalBasis::ClebschGordan cg;		
		return MyCoefficient(lp,l) * cg(l,1,0,0,lp,0) * (cg(l,1,m,-1,lp,mp) + cg(l,1,m,1,lp,mp));
	}

	//testing testing
	struct magnitude {
  		bool operator()(const double& x, const double& y)
    	const
    	{ return abs(x) < abs(y); }
	};




};



#endif
