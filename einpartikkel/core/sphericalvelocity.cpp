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
			if (m != mp) continue;
			if (std::abs(l - lp) != 1) continue;

			double E = LaserHelper::E(lp, m) * LaserHelper::kronecker(l, lp+1);
			double F = LaserHelper::F(lp, m) * LaserHelper::kronecker(l, lp-1);

			double coupling = (E + F);


			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				data(index) = - IM * coupling;
			}
		}
	}
};



//template<int Rank>
//class CustomPotential_LaserVelocityDerivativeR1
//{
//public:
//	typedef blitz::Array<int, 2> BasisPairList;
//
//private:
//	BasisPairList AngularBasisPairs;
//	int AngularRank;
//	int RadialRank1;
//	int RadialRank2;
//
//public:
//	CustomPotential_LaserVelocityDerivativeR1() {}
//	virtual ~CustomPotential_LaserVelocityDerivativeR1() {}
//
//	void ApplyConfigSection(const ConfigSection &config)
//	{
//		config.Get("radial_rank1", RadialRank1);
//		config.Get("radial_rank2", RadialRank2);
//		config.Get("angular_rank", AngularRank);
//	}
//
//	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
//	{
//		if (rank != AngularRank)
//		{
//			throw std::runtime_error("Only angular rank supports basis pairs");
//		}
//		AngularBasisPairs.reference(basisPairs.copy());
//	}
//
//	BasisPairList GetBasisPairList(int rank)
//	{
//		if (rank == AngularRank)
//			return AngularBasisPairs;
//		else
//			return BasisPairList();
//	}
//
//	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
//	{
//		using namespace CoupledSpherical;
//
//		typedef CombinedRepresentation<Rank> CmbRepr;
//		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
//		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
//	
//		int r1Count = data.extent(RadialRank1);
//		int r2Count = data.extent(RadialRank2);
//		int angCount = data.extent(AngularRank);
//
//		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
//		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
//		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);
//
//		ClebschGordan cg;
//
//		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
//		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
//		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");
//
//		cplx IM(0,1.0);
//
//		data = 0;
//		blitz::TinyVector<int, Rank> index;
//	
//		for (int angIndex=0; angIndex<angCount; angIndex++)
//		{
//			index(AngularRank) = angIndex;
//
//			int leftIndex = angBasisPairs(angIndex, 0);
//			int rightIndex = angBasisPairs(angIndex, 1);
//	
//			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
//			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);
//
//			//"Left" quantum numbers
//			int l1 = left.l1;
//			int l2 = left.l2;
//			int L = left.L;
//			int M = left.M;
//			
//			//"Right" quantum numbers (Mp = M)
//			int l1p = right.l1;
//			int l2p = right.l2;
//			int Lp = right.L;
//			int Mp = right.M;
//
//			//Rough selection rule
//			if (M != Mp) continue;
//			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));
//
//			double coupling1 = 0;
//			//double coupling2 = 0;
//
//			for (int m=-lStop; m<=lStop; m++)
//			{
//				double curCg = cg(l1, l2, m, M-m, L, M) * cg(l1p, l2p, m, Mp-m, Lp, Mp);
//
//				if (LaserHelper::kronecker(l2, l2p) != 0)
//				{
//					if (std::abs(m) <= l1 && std::abs(m) <= l1p)
//					{
//						double E = LaserHelper::E(l1p, m) * LaserHelper::kronecker(l1, l1p+1);
//						double F = LaserHelper::F(l1p, m) * LaserHelper::kronecker(l1, l1p-1);
//
//						coupling1 += (E + F) * curCg;
//					}
//				}
//			}
//
//			for (int ri1=0; ri1<r1Count; ri1++)
//			{
//				index(RadialRank1) = ri1;
//				//double r1 = localr1(ri1);
//
//				for (int ri2=0; ri2<r2Count; ri2++)
//				{
//					index(RadialRank2) = ri2;
//					//double r2 = localr2(ri2);
//				
//					data(index) = - IM * (coupling1);
//				}
//			}
//		}
//	}
//};
