#include "G4Interfaces/PhysicsListPyExport.hh"
// #include "G4PhysListsSilicon/PhysicsListNIEL.hh"
#include "G4PhysListsSilicon/QGSP_BERT_HP_mod.hh"

PYTHON_MODULE
{
  PhysicsListPyExport::exportPhysList<QGSP_BERT_HP_mod>();
}
