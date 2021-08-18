#include "G4PhysListsSilicon/PhysicsListNIEL.hh"
#include "G4PhysListsSilicon/PhysListEmStandardNR.hh"
#include "G4PhysListsSilicon/G4EmStandardPhysics_option4_mod.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4DecayPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4EmExtraPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonElasticPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"

// #include "G4OpticalPhysics.hh"
#include "G4PhysListsSilicon/QGSP_BIC_HP_mod.hh"
// #include "QGSP_BERT_HP.hh"
// #include "G4Decay.hh"
// #include "G4HadronElasticPhysics.hh"
// #include "G4EmStandardPhysicsSS.hh"
// #include "G4Version.hh"
// #include "G4HadronPhysicsQGSP_BIC_HP.hh"

PhysicsListNIEL :: PhysicsListNIEL() : G4VModularPhysicsList()
{
  RegisterPhysics( new G4EmStandardPhysics_option1());
  // RegisterPhysics( new G4EmStandardPhysics_option1());
  // RegisterPhysics( new G4EmStandardPhysics_option4_mod());
  // RegisterPhysics( new G4DecayPhysics("decay",1) );
  // RegisterPhysics( new G4RadioactiveDecayPhysics());

  // RegisterPhysics( new QGSP_BIC_HP_mod());
  // RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP());
  // RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP("hadronQGSP_BIC_HP",true));
  // RegisterPhysics( new G4EmExtraPhysics("extra EM"));
  // RegisterPhysics( new G4HadronElasticPhysicsHP(1));
  // RegisterPhysics( new G4StoppingPhysics("stopping",1,true));
  // RegisterPhysics( new G4IonPhysics());
  // RegisterPhysics( new G4IonElasticPhysics());
  // RegisterPhysics( new G4EmStandardPhysics(1));
  // RegisterPhysics( new G4EmStandardPhysicsSS(1));
  // RegisterPhysics( new PhysListEmStandardNR(1) );
  // RegisterPhysics( new G4OpticalPhysics());

  /*  G4Scintillation * theScintillationProcess = new G4Scintillation("Scintillation");
  theScintillationProcess->SetScintillationYieldFactor(1.);
  theScintillationProcess->SetTrackSecondariesFirst(true);
  G4EmSaturation*emSaturation=G4LossTableManager::Instance()-> EmSaturation();
  theScintillationProcess->AddSaturation(emSaturation);
  */
}

PhysicsListNIEL::~PhysicsListNIEL()
{
}
