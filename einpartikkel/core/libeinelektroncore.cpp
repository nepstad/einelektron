// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libeinelektroncore)
{
    Export_wrapper();
}
