#include "sphericalbase.h"

/*
 * Angular Kinetic Energy
 */
template<int Rank>
class CustomPotential_AngularKineticEnergy_Spherical : public CustomPotentialSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;
	double Mass;

public:
	CustomPotential_AngularKineticEnergy_Spherical() {}
	virtual ~CustomPotential_AngularKineticEnergy_Spherical() {}

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialSphericalBase<Rank>::ApplyConfigSection(config);
		config.Get("mass", Mass);
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		SphericalHarmonicBasisRepresentation::Ptr angRepr = boost::static_pointer_cast< SphericalHarmonicBasisRepresentation >(repr->GetRepresentation(this->AngularRank));
	
		int rCount = data.extent(this->RadialRank);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr = psi->GetRepresentation()->GetLocalGrid(this->RadialRank);
		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			//Quantum numbers
			LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
			LmIndex right = angRepr->Range.GetLmIndex(rightIndex);	

			double centrifugalTerm = 0;

			if (!((left.l == right.l) && (left.m == right.m)))
			{
				continue;
			}
			centrifugalTerm = left.l * (left.l + 1);

			for (int ri=0; ri<rCount; ri++)
			{
				index(this->RadialRank) = ri;
				double r = localr(ri);
				data(index) = centrifugalTerm / (2.0 * Mass * r * r);
			}
		}
	}
};

/*
 * Angular Kinetic Energy
 */
template<int Rank>
class SphericalKineticEnergyEvaluator : public CustomPotentialSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	SphericalKineticEnergyEvaluator() {}
	virtual ~SphericalKineticEnergyEvaluator() {}

	double Mass;

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialSphericalBase<Rank>::ApplyConfigSection(config);
		config.Get("mass", Mass);
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace SphericalBasis;

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
	
		for (int ri=0; ri<rCount; ri++)
		{
			index(this->RadialRank) = ri;
			double r = localr(ri);
	
			for (int angIndex=0; angIndex<angCount; angIndex++)
			{
				index(this->AngularRank) = angIndex;

				int leftIndex = angBasisPairs(angIndex, 0);
				int rightIndex = angBasisPairs(angIndex, 1);
	
				LmIndex left = angRepr->Range.GetLmIndex(leftIndex);
				LmIndex right = angRepr->Range.GetLmIndex(rightIndex);
	
				if (left.l != right.l) continue;
				if (left.m != right.m) continue;

				double V = left.l * (left.l + 1.) / (2. * Mass * r * r);
	
				data(index) = V;
			}
		}
	}
};

