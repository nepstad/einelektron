#include "sphericalbase.h"


#include "sphericalvelocityBodyX.h"

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

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("angular_rank", AngularRank);
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
			if (std::abs(m - mp) != 1) continue;
			if (std::abs(l - lp) != 1) continue;

			double coupling = velocityHelperX::sphericalvelocityBodyX(lp,mp,l,m);


			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				double r = localr(ri);
				
				data(index) = - IM * coupling/r;
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

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("angular_rank", AngularRank);
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
			coupling = velocityHelperX::I1integral(l,m,lp,mp);

			//cout << l << "\t"  << m << "\t" << lp << "\t"  << mp << "\t" << coupling << "\t deriv" << endl;

			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				data(index) = - IM * coupling;
			}
		}
	}
};

