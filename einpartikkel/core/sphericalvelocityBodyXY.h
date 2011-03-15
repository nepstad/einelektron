#ifndef SPHERICALVELOCITYYBODYXY_H
#define SPHERICALVELOCITYYBODYXY_H
#define numDigitsPrecision1 50


#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include "sphericalbase.h"


#include <arprec/mp_real.h>

class velocityHelperXY
{
public:
	static double sphericalvelocityBodyXY(int l, int m, int lp, int mp,vector<mp_real> vv, bool ImX)	
	{

		// Arbitrary precision library.
		// Initialization should be set to desired precision plus two
		mp::mp_init(numDigitsPrecision1 + 2); 
		mp::mpsetprec(numDigitsPrecision1 ); 
		mp::mpsetoutputprec(numDigitsPrecision1 ); 
		cout.precision(numDigitsPrecision1 ) ; 


		// Define the "eps"; used in zero check.
		double eps = std::pow(10.,-15);

		double coupling = 0;
		double I2 = 0; 
		double dlta1 = 0;
		double dlta2 = 0;
		double dlta_m;

		double J1 = 0;
		double J2 = 0;

		dlta1 = Mykronecker(mp - 1, m);
		dlta2 = Mykronecker(mp + 1, m);	


		
		/*
		 * <Y_lm|\sin\phi\sin\theta|Y_lm> 
		 */

		double I1 = 0;
		if (ImX == true)
		{
			I1 += I1integralX(lp,mp,l,m);
		}
		else
		{
			I1 += I1integralY(lp,mp,l,m);
		}




		/*
		 * <Y_lm|\frac{\cos\phi}{\sin\theta}\frac{\partial}{\partial\phi}|Y_lm>
		 */

		//Check if call is necessary.	
		bool lnNorm_ok = checklnLegendreNorm(l,m,lp,mp);
		bool lnK1_ok = checklnK1(lp,std::abs(mp),l,std::abs(m));


		if ((lnNorm_ok == true) && (lnK1_ok == true))
		{
			vector<double> lnNorm = lnLegendreNormVect(l,m,lp,mp);
			vector<double> lnK1Int = lnK1vect(lp,std::abs(mp),l,std::abs(m));
		
			//Join the vectors and sort
			vector<double> v = joinTwoVectors(lnNorm,lnK1Int);
			sort(v.begin(),v.end(),magnitude());

			I2 += exp(sumVector(v));

			if (ImX == true)
			{
				I2 *= (dlta1 - dlta2);
			}
			else
			{
				I2 *= (dlta1 + dlta2); 
			}

			I2 *= 0.5 * m * pow(-1,0.5 * (m + mp + std::abs(m) + std::abs(mp)));
		}

	
		/*
		 * <Y_lm|\sin\phi\cos\theta\crac{\partial}{\partial\theta}|Y_lm>
		 * This integral is evaluated in two partis: J1 and J2.
		 */

		//Check if call is necessary.
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
				vector<double> lnNorm1 = lnLegendreNormVect(l-2,m,lp,mp);
				vector<double> lnK1Int_1 = lnK1vect(lp,std::abs(mp),l-2,std::abs(m));
				vector<double> lnFconst = lnFvect(l,m);
			
				//Join the vectors and sort
				vector<double> b = joinThreeVectors(lnNorm1,lnK1Int_1,lnFconst);
				sort(b.begin(),b.end(),magnitude());
				
				J1 += exp(sumVector(b));

			}
		}

		if ((lnNorm2_ok == true) && (lnK1Int_2_ok == true))
		{
			vector<double> lnNorm2 = lnLegendreNormVect(l,m,lp,mp);
			vector<double> lnK1Int_2 = lnK1vect(lp,std::abs(mp),l,std::abs(m));

			vector<double> c = joinTwoVectors(lnNorm2,lnK1Int_2);


			if (lnG_ok == true)
			{
				vector<double> lnGconst = lnGvect(l,m);

				//Join the vectors and sort
				vector<double> d = joinTwoVectors(c,lnGconst);
				sort(d.begin(),d.end(),magnitude());

				J1 += exp(sumVector(d));
			}

			if (lnH_ok == true)
			{
				vector<double> lnHconst = lnHvect(l,m);

				//Join the vectors and sort
				vector<double> f = joinTwoVectors(c,lnHconst);
				sort(f.begin(),f.end(),magnitude());
				
				J1 += exp(sumVector(f));

			}
		}

		if ((lnNorm3_ok == true) && (lnK1Int_3_ok == true))
		{
			if (lnI_ok == true)
			{
				vector<double> lnNorm3 = lnLegendreNormVect(l+2,m,lp,mp);
				vector<double> lnK1Int_3 = lnK1vect(lp,std::abs(mp),l+2,std::abs(m));
				vector<double> lnIconst = lnIvect(l,m);

				//Join the vectors and sort
				vector<double> h = joinThreeVectors(lnNorm3,lnK1Int_3,lnIconst);
				sort(h.begin(),h.end(),magnitude());
				
				J1 += exp(sumVector(h));

			}
		}

		if (ImX == true)
		{
			J1 *= (dlta1 + dlta2);
		}
		else
		{
			J1 *= (-dlta1 + dlta2);
		}
		J1 *= std::abs(m) * 0.5 * pow(-1,0.5 * (m + mp + std::abs(m) + std::abs(mp)));


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
			vector<double> lnEconst = lnEvect(l,m);

			if ((lnNorm4_ok == true) && (lnJ_ok == true))
			{
				vector<double> lnNorm4 = lnLegendreNormVect(l-1,m+dlta_m,lp,mp);
				vector<double> lnJconst = lnJvect(l,m,eps);

				vector<double> s = joinThreeVectors(lnNorm4,lnJconst,lnEconst);

				J2 += lnSumK2(l-1,std::abs(m+dlta_m),lp,std::abs(mp),s,vv);
			}
		 
			if ((lnNorm5_ok == true) && (lnK_ok == true))
			{
				vector<double> lnNorm5 = lnLegendreNormVect(l+1,m+dlta_m,lp,mp);
				vector<double> lnKconst = lnKvect(l,m);

				vector<double> v = joinThreeVectors(lnNorm5,lnKconst,lnEconst);

				J2 += lnSumK2(l+1,std::abs(m+dlta_m),lp,std::abs(mp),v,vv);
			}
		}

		J2 *= dlta_m * 0.25;

		if (ImX == true)
		{
			J2 *= (dlta1 + dlta2);
		}
		else
		{
			J2 *= (-dlta1 + dlta2); 
		}
		J2 *= std::pow(-1.,0.5 * (mp + std::abs(mp) + m + dlta_m + std::abs(m + dlta_m)));

		if (ImX == true)
		{
			coupling += (-I1 - I2 + J1 + J2);
		}
		else
		{
			coupling += (-I1 + I2 + J1 + J2);
		}

		return coupling;
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


	static vector<double> lnFvect(double l, double m)
	{
		vector<double> v;
		v.push_back(0.5 * log(l+m));
		v.push_back(0.5 * log(l-m));
		v.push_back(0.5 * log(l+m-1));
		v.push_back(0.5 * log(l-m-1));
		v.push_back((-0.5) * log(2.*l+1.));
		v.push_back((-1.) * log(2.*l-1.));
		v.push_back((-0.5) * log(2.*l-3.));
		
		return v;
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


	static vector<double> lnGvect(double l, double m)
	{
		vector<double> v;
		v.push_back(log(l+m));
		v.push_back(log(l-m));
		v.push_back((-1.) * log(2*l+1));
		v.push_back((-1.) * log(2*l-1));

		return v;
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


	static vector<double> lnHvect(double l, double m)
	{
		vector<double> v;
		v.push_back(log(l+m+1));
		v.push_back(log(l-m+1));
		v.push_back((-1.) * log(2*l+1));
		v.push_back((-1.) * log(2*l+3));

		return v;
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


	static vector<double> lnIvect(double l, double m)
	{
		vector<double> v;
		v.push_back(0.5 * log(l+m+1));
		v.push_back(0.5 * log(l-m+1));
		v.push_back(0.5 * log(l+m+2));
		v.push_back(0.5 * log(l-m+2));
		v.push_back((-0.5) * log(2.*l+1.));
		v.push_back((-1.) * log(2.*l+3));
		v.push_back((-0.5) * log(2.*l+5));
		return v;
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


	static vector<double> lnK1vect(int l, int m, int p, int q)
	{
		int nu = std::min(l,p);
		int mu = std::min(m,q);
	
		vector<double> v;
		v.push_back(gsl_sf_lnfact(nu+mu));
		v.push_back((-1.)*gsl_sf_lnfact(nu-mu));
		return v;
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


	static double lnSumK2(int l, int m, int p, int q,vector<double> lnsum, vector<mp_real> vv)
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

				mp_real sumArbPrec;
				sumArbPrec = mp_real(0.0);

				for (int i=0; i<=outerMax; i++)
				{

					for (int j=0; j<=innerMax; j++)
					{
						//This part must be done in arbitrary precision.
						mp_real gammaArg;
						gammaArg = vv[fixIndex(.5 * (l+p-m-q -2.*(i+j)+1.))] + vv[fixIndex(.5 * (m+q+2.*(i+j+1.)))] - vv[fixIndex(.5*(l+p+3.))];
						gammaArg += vv[fixIndex(p+q+1)]-vv[fixIndex(q+j+1)]-vv[fixIndex(j+1)]-vv[fixIndex(p-q-2*j+1)];
						gammaArg += vv[fixIndex(l+m+1)]-vv[fixIndex(m+i+1)]-vv[fixIndex(i+1)]-vv[fixIndex(l-m-2*i+1)];
						gammaArg -= log(2.0) * (mp_real(m) + mp_real(q) + mp_real(2.0) * (mp_real(i) + mp_real(j)));

						//cout << q+j+1  << " g" << endl;
			
						sort(lnsum.begin(),lnsum.end(),magnitude());

						mp_real sumAll = mp_real(0.0);
						for (int idx = 0; idx < lnsum.size(); idx++)
						{
							sumAll += mp_real(lnsum[idx]);
						}

						mp_real out;
						out = exp(sumAll + gammaArg);
						out *= pow(-1,i+j);

						sumArbPrec += out;
					}

				}

				outerSum = arbPrecToDouble(sumArbPrec);
			
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


	static double sumVector(vector<double> a)
	{
		double sum = 0;
		for (int idx = 0; idx < a.size(); idx++)
		{
					sum += a[idx];
		}		
		return sum;
	}


	static vector<double> joinThreeVectors(vector<double> a, vector<double> b, vector<double> c)
	{
		vector<double> tmpVec;
		tmpVec.reserve(a.size() + b.size());
		tmpVec.insert(tmpVec.end(), a.begin(), a.end()); 
		tmpVec.insert(tmpVec.end(), b.begin(), b.end());

		vector<double> rsltVec;
		rsltVec.reserve(tmpVec.size() + c.size());
		rsltVec.insert(rsltVec.end(), tmpVec.begin(), tmpVec.end()); 
		rsltVec.insert(rsltVec.end(), c.begin(), c.end());

		return rsltVec;
	}


	static vector<double> joinTwoVectors(vector<double> a, vector<double> b)
	{
		vector<double> rsltVec;
		rsltVec.reserve(a.size() + b.size());
		rsltVec.insert(rsltVec.end(), a.begin(), a.end()); 
		rsltVec.insert(rsltVec.end(), b.begin(), b.end());

		return rsltVec;
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


	static double I1integralX(double lp,double mp,double l,double m)
	{
		SphericalBasis::ClebschGordan cg;		
		return MyCoefficient(lp,l) * cg(l,1,0,0,lp,0) * (cg(l,1,m,-1,lp,mp) - cg(l,1,m,1,lp,mp));
	}

	static double I1integralY(double lp,double mp,double l,double m)
	{
		SphericalBasis::ClebschGordan cg;		
		return MyCoefficient(lp,l) * cg(l,1,0,0,lp,0) * (cg(l,1,m,-1,lp,mp) + cg(l,1,m,1,lp,mp));
	}


	struct magnitude {
		//Auxiliary function when sorting
  		bool operator()(const double& x, const double& y)
    	const
    	{ return abs(x) < abs(y); }
	};


    static double arbPrecToDouble(mp_real a)
    {
        //Convert a real arb. prec (mp_real) to double.

		//Represent mp_real as string
        string b = a.to_string(20);

		//Identify exponent and significand
        int f = b.find('^');
        string sub = b.substr(f+1,b.size()+f+1);
        int ff = sub.find('x');
        string subExponent = sub.substr(0,ff);
        string subSignificand = sub.substr(ff+1,sub.size()-ff);

        int exponent;
        double significand;

		//Convert from string to int and double
        stringstream ss(subExponent);
        ss >> exponent;
        stringstream dd(subSignificand);
        dd >> significand;

		//Recombine
        double s = significand * std::pow(10,exponent);

        return s;
    }


	static int fixIndex(double a)
	{
		//Maps input to correct index in 
		//log-gamma vector.
		int b = int(2.0 * a);
		return b;
	}


	static vector<mp_real> genLogGamma(int lmax)
	{
		//Set up log-gamma vector.
		vector<mp_real> v;
		v.push_back(mp_real(0.0));

		double a = 2.0 * lmax + 1; 	

		int maxIdx = 2 * a;

		for (int idx = 1; idx <= maxIdx; idx++)
		{
			double b = idx / 2.0;
			v.push_back(log(gamma(mp_real(b))));
		}
		return v;
	}

};



#endif
