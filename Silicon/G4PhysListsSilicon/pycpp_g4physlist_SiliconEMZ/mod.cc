#include "G4Interfaces/PhysicsListPyExport.hh"
#include "G4PhysListsSilicon/PhysicsListSilicon.hh"

PYTHON_MODULE
{
  PhysicsListPyExport::exportPhysList<PhysicsListSilicon>();
}
