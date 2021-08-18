#include "G4PhysListsSilicon/PhysicsListSilicon.hh"
#include "G4PhysListsSilicon/PhysListEmStandardNR.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4OpticalPhysics.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BERT_HP.hh"
#include "G4Decay.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4Version.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"

PhysicsListSilicon :: PhysicsListSilicon() : G4VModularPhysicsList()
{
  // RegisterPhysics( new G4EmStandardPhysics(1));
  // RegisterPhysics( new G4EmStandardPhysicsSS(1));
  // RegisterPhysics( new G4EmExtraPhysics("extra EM"));
  RegisterPhysics( new PhysListEmStandardNR(1) );
  // RegisterPhysics( new G4HadronElasticPhysics(1));
  RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP("hadronQGSP_BIC_HP",true));
  RegisterPhysics( new G4DecayPhysics("decay",1) );
  // RegisterPhysics( new G4EmStandardPhysics_option4(1));
  // RegisterPhysics( new G4StoppingPhysics("stopping",1,true));
  // RegisterPhysics( new G4IonPhysics("ion"));
  // RegisterPhysics( new G4OpticalPhysics());

  /*  G4Scintillation * theScintillationProcess = new G4Scintillation("Scintillation");
  theScintillationProcess->SetScintillationYieldFactor(1.);
  theScintillationProcess->SetTrackSecondariesFirst(true);
  G4EmSaturation*emSaturation=G4LossTableManager::Instance()-> EmSaturation();
  theScintillationProcess->AddSaturation(emSaturation);
  */
}

PhysicsListSilicon::~PhysicsListSilicon()
{
}
