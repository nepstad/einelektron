#include "sphericalbase.h"

/*
 * Potential evaluator for linearly polarized length gauge electric field (z-direction)
 */
template<int Rank>
class CustomPotential_LaserLength_Z : public CustomPotentialSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotential_LaserLength_Z() {}
	virtual ~CustomPotential_LaserLength_Z() {}

	cplx Charge;

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialSphericalBase<Rank>::ApplyConfigSection(config);
		//charge with sign
		config.Get("charge", Charge);
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

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);
	
			if (std::abs(left.l - right.l) != 1) continue;
			if (std::abs(left.m - right.m) != 0) continue;

			// "Left" quantum numbers
			int l = left.l;
			int m = left.m;
			
			// "Right" quantum numbers 
			int lp = right.l;
			int mp = right.m;
			
			//Angular matrix element
			double I = cg(lp,1,0,0,l,0) * cg(lp,1,mp,0,l,m);
			I *= Coefficient(lp, l);

			//Radial matrix element 
			for (int ri=0; ri<rCount; ri++)
			{
				double r = localr(ri);
				index(this->RadialRank) = ri;

				data(index) += I * r;
			}
		}

		//Charge scaling from config
		data *= (-1.) * Charge;
	}

	static double Coefficient(int a, int b)
	{
		return std::sqrt((2 * a + 1.0 ) / (2 * b + 1.0));
	}
};


/*
 * Potential evaluator for linearly polarized length gauge electric field (x-direction)
 */
template<int Rank>
class CustomPotential_LaserLength_X : public CustomPotentialSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotential_LaserLength_X() {}
	virtual ~CustomPotential_LaserLength_X() {}

	cplx Charge;

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialSphericalBase<Rank>::ApplyConfigSection(config);
		//charge with sign
		config.Get("charge", Charge);
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

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);
	
			if (std::abs(left.l - right.l) != 1) continue;
			if (std::abs(left.m - right.m) != 1) continue;

			// "Left" quantum numbers
			int l = left.l;
			int m = left.m;
			
			// "Right" quantum numbers 
			int lp = right.l;
			int mp = right.m;
			
			//Angular matrix element
			double I = cg(lp,1,0,0,l,0) * (cg(lp,1,mp,-1,l,m) - cg(lp,1,mp,1,l,m));
			I *= Coefficient(lp, l);

			//Radial matrix element 
			for (int ri=0; ri<rCount; ri++)
			{
				double r = localr(ri);
				index(this->RadialRank) = ri;

				data(index) += I * r;
			}
		}

		//Charge scaling form config
		data *= (-1.) * Charge;
	}

	static double Coefficient(int a, int b)
	{
		return std::sqrt((2 * a + 1.0 ) / (2*(2 * b + 1.0)));
	}
};


/*
 * Potential evaluator for linearly polarized length gauge electric field (y-direction)
 */
template<int Rank>
class CustomPotential_LaserLength_Y : public CustomPotentialSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotential_LaserLength_Y() {}
	virtual ~CustomPotential_LaserLength_Y() {}

	cplx Charge;

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialSphericalBase<Rank>::ApplyConfigSection(config);
		//charge with sign
		config.Get("charge", Charge);
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

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
		cplx IM(0,1.0);

	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);
	
			if (std::abs(left.l - right.l) != 1) continue;
			if (std::abs(left.m - right.m) != 1) continue;

			// "Left" quantum numbers
			int l = left.l;
			int m = left.m;
			
			// "Right" quantum numbers 
			int lp = right.l;
			int mp = right.m;

			
			//Angular matrix element
			double I = cg(lp,1,0,0,l,0) * (cg(lp,1,mp,-1,l,m) + cg(lp,1,mp,1,l,m));
			I *= Coefficient(lp, l);

			//Radial matrix element 
			for (int ri=0; ri<rCount; ri++)
			{
				double r = localr(ri);
				index(this->RadialRank) = ri;

				data(index) += I * r;
			}
		}
	
		data *= IM;
		
		//Charge scaling from config
		data *= (-1.) * Charge;
	}

	static double Coefficient(int a, int b)
	{
		return std::sqrt((2 * a + 1.0 ) / (2*(2 * b + 1.0)));
	}
};

