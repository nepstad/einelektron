
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <analysis.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libeinpartikkelanalysis)
{
    def("SetRadialCoulombWave",  SetRadialCoulombWave);
}

