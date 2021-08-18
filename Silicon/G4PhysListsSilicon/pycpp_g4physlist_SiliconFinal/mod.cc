#include "G4Interfaces/PhysicsListPyExport.hh"
// #include "G4PhysListsSilicon/PhysicsListNIEL.hh"
#include "G4PhysListsSilicon/QGSP_INCLXX_HP_mod.hh"

PYTHON_MODULE
{
  PhysicsListPyExport::exportPhysList<QGSP_INCLXX_HP_mod>();
}
