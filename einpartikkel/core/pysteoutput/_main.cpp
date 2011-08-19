// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_modules_einpartikkel_core_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libeinelektroncore)
{
    Export_modules_einpartikkel_core_wrapper();
}
