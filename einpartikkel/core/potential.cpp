#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

using namespace blitz;

template<int Rank>
class KineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double mass;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return - 1. / (2. * mass);
	}

};


template<int Rank>
class RadialHarmonicPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		return 0.5*r*r;
	}
};

template<int Rank>
class CoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;
	double Charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
		config.Get("charge", Charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		return Charge / r;
	}
};


template<int Rank>
class SingleActiveElectronPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int angularRank;
	int radialRank;

	//Potential parameters
	double Z;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double a6;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("a1", a1);
		config.Get("a2", a2);
		config.Get("a3", a3);
		config.Get("a4", a4);
		config.Get("a5", a5);
		config.Get("a6", a6);

		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 *
	 * Some general tips for max efficiency:
	 * - If possible, move static computations to ApplyConfigSection.
	 * - Minimize the number of branches ("if"-statements are bad)
	 * - Minimize the number of function calls (sin, cos, exp, are bad)
	 * - Long statements can confuse the compiler, consider making more 
	 *   simpler statements
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = std::fabs(pos(radialRank));
		
		double V = -(Z + a1 * exp(-a2 * r) + a3 * r * exp(-a4 * r)
			+ a5 * exp(-a6 * r)) / r;

		return V;
	}
};


template<int Rank>
class OverlapPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int angularRank;
	int radialRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return 1.0;
	}
};

template<int Rank>
class ComplexAbsorbingPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int radialRank;
	double scalingReal;
	double scalingImag;
	double factorReal;
	double factorImag;
	double absorberStart;
	double absorberLength;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", radialRank);
		config.Get("absorber_start", absorberStart);
		config.Get("absorber_length", absorberLength);
		config.Get("scaling_real", scalingReal);
		config.Get("scaling_imag", scalingImag);
		config.Get("factor_real", factorReal);
		config.Get("factor_imag", factorImag);
	}

	/*
	 * Called for every grid point 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		cplx V = 0;
		if (r > absorberStart)
		{
			double curLength = (r - absorberStart) / absorberLength;
			double Vr = factorReal * std::pow(curLength, scalingReal);
			double Vi = factorImag * std::pow(curLength, scalingImag);
			V += cplx(Vr , Vi);
		}
		return V;
	}
};


/*
 * Manolopoulos transmission free absorbing potential
 *
 *   See D. E. Manolopoulos, J. Chem. Phys 117, 9552, 2002.
 *
 *   energy_cutoff: wavepacket components with energy greater than this is
 *                  absorbed
 *   grid_max:      last point on the grid absorberwise
 *   delta:         accuracy parameter, determines width of absorber
 *
 */
template<int Rank>
class ManolopoulosAbsorber : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Rank info
	int RadialRank;

	//Potential parameters
	double EnergyCutoff;
	double GridMax;
	double Delta;
	double Start;

	//CAP constants
	double A, B, C;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("energy_cutoff", EnergyCutoff);
		config.Get("grid_max", GridMax);
		config.Get("delta", Delta);
		config.Get("radial_rank", RadialRank);

		A = 0.112449;
		B = 0.0082735;
		C = 2.62206;

		//Calculate absorber start
		double kmin = std::sqrt(2*EnergyCutoff);
		Start = GridMax - C / (2 * Delta * kmin);
		
		//cout << "Absorber starts at " << Start << std::endl;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(RadialRank);
		cplx cap = 0;

		// Calculate absorber
		if (r >= Start)
		{
			double u = 2 * Delta * std::sqrt(2*EnergyCutoff) * (r - Start);
			double y = A*u - B*u*u*u + 4.0/((C-u)*(C-u)) - 4.0/((C+u)*(C+u));
			cap += -I * EnergyCutoff * y;
		}

		return cap;
	}
};

