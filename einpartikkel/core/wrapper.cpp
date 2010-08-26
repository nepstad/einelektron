
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential.cpp>
#include <spherical.cpp>
#include <sphericallength.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct SphericalKineticEnergyEvaluator_2_Wrapper: SphericalKineticEnergyEvaluator<2>
{
    SphericalKineticEnergyEvaluator_2_Wrapper(PyObject* py_self_, const SphericalKineticEnergyEvaluator<2>& p0):
        SphericalKineticEnergyEvaluator<2>(p0), py_self(py_self_) {}

    SphericalKineticEnergyEvaluator_2_Wrapper(PyObject* py_self_):
        SphericalKineticEnergyEvaluator<2>(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        SphericalKineticEnergyEvaluator<2>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,2> p0, boost::shared_ptr<Wavefunction<2> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,2> p0, boost::shared_ptr<Wavefunction<2> > p1, std::complex<double> p2, double p3) {
        SphericalKineticEnergyEvaluator<2>::UpdatePotentialData(p0, p1, p2, p3);
    }

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotentialSphericalBase<2>::SetBasisPairs(p0, p1);
    }

    blitz::Array<int,2> GetBasisPairList(int p0) {
        return call_method< blitz::Array<int,2> >(py_self, "GetBasisPairList", p0);
    }

    blitz::Array<int,2> default_GetBasisPairList(int p0) {
        return CustomPotentialSphericalBase<2>::GetBasisPairList(p0);
    }

    PyObject* py_self;
};

struct CustomPotential_LaserLength_Z_2_Wrapper: CustomPotential_LaserLength_Z<2>
{
    CustomPotential_LaserLength_Z_2_Wrapper(PyObject* py_self_, const CustomPotential_LaserLength_Z<2>& p0):
        CustomPotential_LaserLength_Z<2>(p0), py_self(py_self_) {}

    CustomPotential_LaserLength_Z_2_Wrapper(PyObject* py_self_):
        CustomPotential_LaserLength_Z<2>(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CustomPotential_LaserLength_Z<2>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,2> p0, boost::shared_ptr<Wavefunction<2> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,2> p0, boost::shared_ptr<Wavefunction<2> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserLength_Z<2>::UpdatePotentialData(p0, p1, p2, p3);
    }

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotentialSphericalBase<2>::SetBasisPairs(p0, p1);
    }

    blitz::Array<int,2> GetBasisPairList(int p0) {
        return call_method< blitz::Array<int,2> >(py_self, "GetBasisPairList", p0);
    }

    blitz::Array<int,2> default_GetBasisPairList(int p0) {
        return CustomPotentialSphericalBase<2>::GetBasisPairList(p0);
    }

    PyObject* py_self;
};

struct CustomPotential_LaserLength_X_2_Wrapper: CustomPotential_LaserLength_X<2>
{
    CustomPotential_LaserLength_X_2_Wrapper(PyObject* py_self_, const CustomPotential_LaserLength_X<2>& p0):
        CustomPotential_LaserLength_X<2>(p0), py_self(py_self_) {}

    CustomPotential_LaserLength_X_2_Wrapper(PyObject* py_self_):
        CustomPotential_LaserLength_X<2>(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        CustomPotential_LaserLength_X<2>::ApplyConfigSection(p0);
    }

    void UpdatePotentialData(blitz::Array<std::complex<double>,2> p0, boost::shared_ptr<Wavefunction<2> > p1, std::complex<double> p2, double p3) {
        call_method< void >(py_self, "UpdatePotentialData", p0, p1, p2, p3);
    }

    void default_UpdatePotentialData(blitz::Array<std::complex<double>,2> p0, boost::shared_ptr<Wavefunction<2> > p1, std::complex<double> p2, double p3) {
        CustomPotential_LaserLength_X<2>::UpdatePotentialData(p0, p1, p2, p3);
    }

    void SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        call_method< void >(py_self, "SetBasisPairs", p0, p1);
    }

    void default_SetBasisPairs(int p0, const blitz::Array<int,2>& p1) {
        CustomPotentialSphericalBase<2>::SetBasisPairs(p0, p1);
    }

    blitz::Array<int,2> GetBasisPairList(int p0) {
        return call_method< blitz::Array<int,2> >(py_self, "GetBasisPairList", p0);
    }

    blitz::Array<int,2> default_GetBasisPairList(int p0) {
        return CustomPotentialSphericalBase<2>::GetBasisPairList(p0);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libpotential)
{
    class_< SphericalKineticEnergyEvaluator<2>, SphericalKineticEnergyEvaluator_2_Wrapper >("SphericalKineticEnergyEvaluator_2", init<  >())
        .def(init< const SphericalKineticEnergyEvaluator<2>& >())
        .def_readwrite("Mass", &SphericalKineticEnergyEvaluator<2>::Mass)
        .def("ApplyConfigSection", (void (SphericalKineticEnergyEvaluator<2>::*)(const ConfigSection&) )&SphericalKineticEnergyEvaluator<2>::ApplyConfigSection, (void (SphericalKineticEnergyEvaluator_2_Wrapper::*)(const ConfigSection&))&SphericalKineticEnergyEvaluator_2_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (SphericalKineticEnergyEvaluator<2>::*)(blitz::Array<std::complex<double>,2>, boost::shared_ptr<Wavefunction<2> >, std::complex<double>, double) )&SphericalKineticEnergyEvaluator<2>::UpdatePotentialData, (void (SphericalKineticEnergyEvaluator_2_Wrapper::*)(blitz::Array<std::complex<double>,2>, boost::shared_ptr<Wavefunction<2> >, std::complex<double>, double))&SphericalKineticEnergyEvaluator_2_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialSphericalBase<2>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialSphericalBase<2>::SetBasisPairs, (void (SphericalKineticEnergyEvaluator_2_Wrapper::*)(int, const blitz::Array<int,2>&))&SphericalKineticEnergyEvaluator_2_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialSphericalBase<2>::*)(int) )&CustomPotentialSphericalBase<2>::GetBasisPairList, (blitz::Array<int,2> (SphericalKineticEnergyEvaluator_2_Wrapper::*)(int))&SphericalKineticEnergyEvaluator_2_Wrapper::default_GetBasisPairList)
    ;

    class_< CustomPotential_LaserLength_Z<2>, CustomPotential_LaserLength_Z_2_Wrapper >("CustomPotential_LaserLength_Z_2", init<  >())
        .def(init< const CustomPotential_LaserLength_Z<2>& >())
        .def_readwrite("Scaling", &CustomPotential_LaserLength_Z<2>::Scaling)
        .def("ApplyConfigSection", (void (CustomPotential_LaserLength_Z<2>::*)(const ConfigSection&) )&CustomPotential_LaserLength_Z<2>::ApplyConfigSection, (void (CustomPotential_LaserLength_Z_2_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_Z_2_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength_Z<2>::*)(blitz::Array<std::complex<double>,2>, boost::shared_ptr<Wavefunction<2> >, std::complex<double>, double) )&CustomPotential_LaserLength_Z<2>::UpdatePotentialData, (void (CustomPotential_LaserLength_Z_2_Wrapper::*)(blitz::Array<std::complex<double>,2>, boost::shared_ptr<Wavefunction<2> >, std::complex<double>, double))&CustomPotential_LaserLength_Z_2_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialSphericalBase<2>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialSphericalBase<2>::SetBasisPairs, (void (CustomPotential_LaserLength_Z_2_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_Z_2_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialSphericalBase<2>::*)(int) )&CustomPotentialSphericalBase<2>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_Z_2_Wrapper::*)(int))&CustomPotential_LaserLength_Z_2_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength_Z<2>::Coefficient)
        .staticmethod("Coefficient")
    ;

    class_< CustomPotential_LaserLength_X<2>, CustomPotential_LaserLength_X_2_Wrapper >("CustomPotential_LaserLength_X_2", init<  >())
        .def(init< const CustomPotential_LaserLength_X<2>& >())
        .def_readwrite("Scaling", &CustomPotential_LaserLength_X<2>::Scaling)
        .def("ApplyConfigSection", (void (CustomPotential_LaserLength_X<2>::*)(const ConfigSection&) )&CustomPotential_LaserLength_X<2>::ApplyConfigSection, (void (CustomPotential_LaserLength_X_2_Wrapper::*)(const ConfigSection&))&CustomPotential_LaserLength_X_2_Wrapper::default_ApplyConfigSection)
        .def("UpdatePotentialData", (void (CustomPotential_LaserLength_X<2>::*)(blitz::Array<std::complex<double>,2>, boost::shared_ptr<Wavefunction<2> >, std::complex<double>, double) )&CustomPotential_LaserLength_X<2>::UpdatePotentialData, (void (CustomPotential_LaserLength_X_2_Wrapper::*)(blitz::Array<std::complex<double>,2>, boost::shared_ptr<Wavefunction<2> >, std::complex<double>, double))&CustomPotential_LaserLength_X_2_Wrapper::default_UpdatePotentialData)
        .def("SetBasisPairs", (void (CustomPotentialSphericalBase<2>::*)(int, const blitz::Array<int,2>&) )&CustomPotentialSphericalBase<2>::SetBasisPairs, (void (CustomPotential_LaserLength_X_2_Wrapper::*)(int, const blitz::Array<int,2>&))&CustomPotential_LaserLength_X_2_Wrapper::default_SetBasisPairs)
        .def("GetBasisPairList", (blitz::Array<int,2> (CustomPotentialSphericalBase<2>::*)(int) )&CustomPotentialSphericalBase<2>::GetBasisPairList, (blitz::Array<int,2> (CustomPotential_LaserLength_X_2_Wrapper::*)(int))&CustomPotential_LaserLength_X_2_Wrapper::default_GetBasisPairList)
        .def("Coefficient", &CustomPotential_LaserLength_X<2>::Coefficient)
        .staticmethod("Coefficient")
    ;

    class_< DynamicPotentialEvaluator<KineticEnergyPotential<2>,2> >("KineticEnergyPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<CoulombPotential<2>,2> >("CoulombPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<CoulombPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<CoulombPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2> >("SingleActiveElectronPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<SingleActiveElectronPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<OverlapPotential<2>,2> >("OverlapPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<OverlapPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<OverlapPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2> >("ComplexAbsorbingPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<2>,2>::CalculateExpectationValue)
    ;

}

