// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_modules_einpartikkel_analysis_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libeinpartikkelanalysis)
{
    Export_modules_einpartikkel_analysis_wrapper();
}
