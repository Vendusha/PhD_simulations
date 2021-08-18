#include "G4Interfaces/PhysicsListPyExport.hh"
#include "G4PhysListsSilicon/PhysicsList.hh"

PYTHON_MODULE
{
  PhysicsListPyExport::exportPhysList<PhysicsList>();
}
