#include "sphericalbase.h"
#include "laserhelper.h"

/* First part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics
 *
 * <Ylm | - \frac{1}{r} \sin \theta \partialdiff{}{\theta} 
 *	      - \frac{\cos \theta}{r} | Yl'm'>
 */
template<int Rank>
class CustomPotential_LaserVelocity
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank;

public:
	CustomPotential_LaserVelocity() {}
	virtual ~CustomPotential_LaserVelocity() {}

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
			if (m != mp) continue;
			if (std::abs(l - lp) != 1) continue;

			double C = LaserHelper::C(lp, m) * LaserHelper::kronecker(l, lp+1);
			double D = LaserHelper::D(lp, m) * LaserHelper::kronecker(l, lp-1);
			double E = LaserHelper::E(lp, m) * LaserHelper::kronecker(l, lp+1);
			double F = LaserHelper::F(lp, m) * LaserHelper::kronecker(l, lp-1);

			double coupling = -(C + D) - (E + F);


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
 * Should be used with first order differentiation in r
 *
 */
template<int Rank>
class CustomPotential_LaserVelocityDerivativeR
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank;

public:
	CustomPotential_LaserVelocityDerivativeR() {}
	virtual ~CustomPotential_LaserVelocityDerivativeR() {}

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
			if (m != mp) continue;
			if (std::abs(l - lp) != 1) continue;

			double E = LaserHelper::E(lp, m) * LaserHelper::kronecker(l, lp+1);
			double F = LaserHelper::F(lp, m) * LaserHelper::kronecker(l, lp-1);

			double coupling = (E + F);


			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				data(index) = - IM* coupling;

				//Charge scaling from config
				data(index) *= (-1.0) * Charge;
			}
		}
	}
};




