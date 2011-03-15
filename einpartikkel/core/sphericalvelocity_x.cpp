#define numDigitsPrecision1 50
#include "sphericalbase.h"
#include <arprec/mp_real.h>

#include "sphericalvelocityBodyXY.h"

/* First part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics
 *
 * <Ylm | - \frac{1}{r} \frac{\sin \phi}{\sin \theta} \partialdiff{}{\phi}
 *  	  + \frac{1}{r} \cos \phi \cos \theta \partialdiff{}{\theta}
 *        - \frac{1}{r} \cos \phi \sin \theta | Yl'm'> 
 */
template<int Rank>
class CustomPotential_LaserVelocity_X
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank;

public:
	CustomPotential_LaserVelocity_X() {}
	virtual ~CustomPotential_LaserVelocity_X() {}

	cplx Charge;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("angular_rank", AngularRank);
		//charge with sign
		config.Get("charge", Charge);		
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef SphericalHarmonicBasisRepresentation SphHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		SphHarmRepr::Ptr angRepr = boost::static_pointer_cast< SphHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int rCount = data.extent(this->RadialRank);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr = psi->GetRepresentation()->GetLocalGrid(this->RadialRank);

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);
		blitz::TinyVector<int, Rank> index;
		data = 0;

		//Arbitrary precision library.
		//Initialization should be set to desired precision plus two
		mp::mp_init(numDigitsPrecision1 + 2); 
		mp::mpsetprec(numDigitsPrecision1 ); 
		mp::mpsetoutputprec(numDigitsPrecision1 ); 
		cout.precision(numDigitsPrecision1 ) ; 


		//Find lmax
		int lmax = 0;
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);

			int l = left.l;
			int lp = right.l;

			int tmpMax = std::max(l,lp);
			if (tmpMax > lmax )
			{
				lmax = tmpMax;
			}
		}

		//Setup Log-Gamma 
		vector<mp_real> v;
		v = velocityHelperXY::genLogGamma(lmax+1);


		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);

			//"Left" quantum numbers
			int l = left.l;
			int m = left.m;
			
			//"Right" quantum numbers
			int lp = right.l;
			int mp = right.m;

			//Selection rules 
			if (std::abs(m - mp) != 1) continue;
			if (std::abs(l - lp) != 1) continue;

			double coupling = velocityHelperXY::sphericalvelocityBodyXY(lp,mp,l,m,v,true);

			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				double r = localr(ri);
				
				data(index) = - IM * coupling/r;

				//Charge scaling from config
				data(index) *= (-1.0) * Charge;
			}
		}
	}
};


/* 
 * Second part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics.
 *
 * <Ylm | \cos \phi \sin \theta | Yl'm'>
 *
 * Should be used with first order differentiation in r
 *
 */
template<int Rank>
class CustomPotential_LaserVelocityDerivativeR_X
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank;

public:
	CustomPotential_LaserVelocityDerivativeR_X() {}
	virtual ~CustomPotential_LaserVelocityDerivativeR_X() {}

	cplx Charge;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("angular_rank", AngularRank);
		//charge with sign
		config.Get("charge", Charge);	
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef SphericalHarmonicBasisRepresentation SphHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		SphHarmRepr::Ptr angRepr = boost::static_pointer_cast< SphHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int rCount = data.extent(this->RadialRank);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr = psi->GetRepresentation()->GetLocalGrid(this->RadialRank);

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);
		blitz::TinyVector<int, Rank> index;
		data = 0;

		SphericalBasis::ClebschGordan cg;

		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);

			//"Left" quantum numbers
			int l = left.l;
			int m = left.m;
			
			//"Right" quantum numbers
			int lp = right.l;
			int mp = right.m;

			//Selection rules 
			if (std::abs(m - mp) != 1) continue;
			if (std::abs(l - lp) != 1) continue;


			if (std::abs(m) > l) continue;
			if (std::abs(mp) > lp) continue;

			double coupling = 0;

			/*
			 * Integral I1.
			 */
			coupling = velocityHelperXY::I1integralX(l,m,lp,mp);


			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				data(index) = - IM * coupling;

				//Charge scaling from config
				data(index) *= (-1.0) * Charge;
			}
		}
	}
};

