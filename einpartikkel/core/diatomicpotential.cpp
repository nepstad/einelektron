//Include interface
#include "sphericalbase.h"

//Including legendre polynomials
#include<gsl/gsl_sf_legendre.h>

//Include laserhelper for kroenecker
#include "laserhelper.h"

//
// Diatomic Coulomb potential: -1 / |r + R/2| - 1/|r - R/2|
//
template<int Rank>
class DiatomicCoulombPotential : public CustomPotentialSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

	DiatomicCoulombPotential() {}
	virtual ~DiatomicCoulombPotential() {}
		
	double R;
	double ThetaR;
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialSphericalBase<Rank>::ApplyConfigSection
			(config);
		config.Get("inter_nuclear_r", R);
		config.Get("theta_inter_nucl", ThetaR);
	}
	
	virtual void UpdatePotentialData(typename blitz::Array<cplx,Rank> data,
	   typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		
		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef SphericalHarmonicBasisRepresentation SphHarmRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >
			(psi->GetRepresentation());

		SphHarmRepr::Ptr angRepr = boost::static_pointer_cast< SphHarmRepr >
			(repr->GetRepresentation(this->AngularRank));

		int rCount = data.extent(this->RadialRank);
		int angCount= data.extent(this->AngularRank);

		blitz::Array<double, 1> localr = psi->GetRepresentation()->
			GetLocalGrid(this->RadialRank);
	
		ClebschGordan cg;
		BasisPairList angBasisPairs = GetBasisPairList
			(this->AngularRank);
		if(data.extent(this->RadialRank) != rCount)
			throw std::runtime_error("Invalid r size");
		if(data.extent(this->AngularRank) != angBasisPairs.extent(0))
			throw std::runtime_error("Invalid ang size");
		
		blitz::TinyVector<int, Rank> index;
		data = 0;
		double 	R_half =  R/2.0;
		//Looping over angular indices
		for(int angIndex = 0; angIndex < angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;
			
			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex= angBasisPairs(angIndex, 1);

			LmIndex left = angRepr->Range.GetLmIndex
				(leftIndex);
			LmIndex right =angRepr->Range.GetLmIndex
				(rightIndex);

			double cosTheta = std::cos(ThetaR);
			
			//"Left" quantum numbers
			int mp = left.m;
			int lp = left.l;
			
			//"Right" quantum numbers
			int m = right.m;
			int l = right.l;

			int minL3 = std::abs(l - lp);
			int maxL3 = l + lp;
			for(int l3 = minL3; l3<=maxL3; l3++)
			{
				//Selection rules
				if(l3 % 2 == 1) continue;

				double l3Coeff = 1.0;
				l3Coeff *= Coefficient(l,lp);
				l3Coeff *= cg(l, l3, 0, 0, lp, 0);
			
				double l3Sum = 0;
				
				for(int m3 = -l3; m3 <= l3; m3++)
				{
					
					//Angular matrix element
					double cur = 1; 
					cur *= 
						gsl_sf_legendre_sphPlm(l3,std::abs(m3),cosTheta);
					cur *= 2.0;
					cur *= CondonShortleyPhase(-m3);
					cur *= MultipoleCoeff(l3);
					cur *= cg(l,l3,m,m3,lp,mp);
					
					l3Sum += cur;
				}
				
				l3Sum *= l3Coeff;
			
				//Radial matrix elements
				double r, rmin, rmax, rfrac;
				for(int ri=0; ri < rCount; ri++)
				{
					r = localr(ri);
					index(this->RadialRank) = ri;
				
					rmin=std::min(r,R_half);
					rmax=std::max(r,R_half);
					rfrac = rmin / rmax;
						
					data(index) += -1. * l3Sum * std::pow(rfrac, l3) / rmax;
				}

			}//End l3 loop
		
		}//End angular index loop

	}//End UpdatePotentialData


	static double Coefficient(int a, int b)
	{
		return std::sqrt((2. *a + 1.) / (2. *b +1.));
	}

	static double MultipoleCoeff(int c)
	{
		return std::sqrt((4. * M_PI) / (2. * c + 1.));
	}
	
	static double CondonShortleyPhase(int m)
	{
		if(m < 0) return 1.0;
		return std::pow(-1.0, m);
	}

};//End class DiatomicCoulomb
