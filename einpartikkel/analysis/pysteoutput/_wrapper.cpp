
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <analysis.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_modules_einpartikkel_analysis_wrapper()
{
    def("SetRadialCoulombWave",  SetRadialCoulombWave);
}

