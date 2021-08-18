#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>

//After you have modified this file as needed for your project you can remove this line and commit <NOCOMMIT>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

int main(int argc,char**argv) {

  //Open .griff file(s) specified on the command line:
  GriffDataReader dr(argc,argv);
  // GriffAnaUtils::SegmentIterator diode(&diode)

  //Extract and dump meta-data:
  dr.setup()->dump();

  //Check (as an example) that the simulated geometry is of the correct type:
  if (dr.setup()->geo().getName()!="G4GeoSilicon/GeoSilicon")
    return 1;
  //energy_init=0.027 keV for 0.2 MeV neutrons, 0.095 for 0.7 MeV neutrons, 2.861 for 3 MeV neutrons, 5.756 for 6 Mev neutrons, 9.7 for 10 MeV neutrons, 13.4 for 14 MeV neutrons, 17.4 for 18 MeV neutrons, 
  auto energy_init = 0.095;
  // auto energy_init= dr.setup()->gen().getParameterDouble("fixed_energy_eV")/1000000;
  // auto energy_init= dr.setup()->gen().getParameterDouble("energy_MeV");
  const std::string particle_type= dr.setup()->gen().getParameterString("particle_type");
  std::cout<< "The generator shooted " << particle_type <<" of the energy "<< energy_init << "MeV ."<<std::endl;
  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto nbins = 256;
  // auto nbins_dim=128;
  // auto n = 0;
  auto h_recoil_PKA_high_Z = hc.book1D("Recoil_spectra_high_Z",nbins, 0, energy_init*1000, "Recoil_spectra_high_Z");
  auto h_recoil_PKA_once = hc.book1D("Recoil_spectra_PKA",nbins, 0, energy_init*1000, "Recoil_spectra_PKA");
  auto h_recoil_PKA_low_Z = hc.book1D("Recoil_spectra_low_Z",nbins, 0, energy_init*1000, "Recoil_spectra_low_Z");
  // auto h_recoil_SKA = hc.book1D("Recoil_spectra_SKA",nbins, 0,energy_init*1000, "Recoil_spectra_SKA");
  auto h_Daughter_PKA_once = hc.bookCounts("PKA_particle_type", "PKA_particle_type");
  auto Si27_once = h_Daughter_PKA_once ->addCounter("Si27");
  auto Si28_once = h_Daughter_PKA_once ->addCounter("Si28");
  auto Si29_once = h_Daughter_PKA_once ->addCounter("Si29");
  auto Si30_once = h_Daughter_PKA_once ->addCounter("Si30");
  auto Si31_once = h_Daughter_PKA_once ->addCounter("Si31");
  auto Al25_once = h_Daughter_PKA_once ->addCounter("Al25");
  auto Al26_once = h_Daughter_PKA_once ->addCounter("Al26");
  auto Al27_once = h_Daughter_PKA_once ->addCounter("Al27");
  auto Al28_once = h_Daughter_PKA_once ->addCounter("Al28");
  auto Al29_once = h_Daughter_PKA_once ->addCounter("Al29");
  auto Al30_once = h_Daughter_PKA_once ->addCounter("Al30");
  auto Mg24_once = h_Daughter_PKA_once ->addCounter("Mg24");
  auto Mg25_once = h_Daughter_PKA_once ->addCounter("Mg25");
  auto Mg26_once = h_Daughter_PKA_once ->addCounter("Mg26");
  auto Mg27_once = h_Daughter_PKA_once ->addCounter("Mg27");
  auto Na23_once = h_Daughter_PKA_once ->addCounter("Na23");
  auto neutron_once = h_Daughter_PKA_once ->addCounter("neutron");
  auto proton_once = h_Daughter_PKA_once ->addCounter("proton");
  auto alpha_once = h_Daughter_PKA_once ->addCounter("alpha");
  auto gamma_once = h_Daughter_PKA_once ->addCounter("gamma");
  auto deuteron_once = h_Daughter_PKA_once ->addCounter("deuteron");
  auto others_once = h_Daughter_PKA_once ->addCounter("others");
  h_recoil_PKA_high_Z -> setXLabel("Energy [keV]");
  h_recoil_PKA_high_Z -> setYLabel("Frequency [PKA-1 keV-1]");
  h_recoil_PKA_low_Z -> setXLabel("Energy [keV]");
  h_recoil_PKA_low_Z -> setYLabel("Frequency [PKA-1 keV-1]");
  h_recoil_PKA_once -> setXLabel("Energy [keV]");
  h_recoil_PKA_once -> setYLabel("Frequency [PKA-1 keV-1]");
  //Loop over events and extract info via Griff interface:
  while (dr.loopEvents()) {
    // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
    // dr.runNumber(),dr.eventNumber(),dr.nTracks());
    while (auto neutron = primary_neutrons.next()){
      if (neutron->nDaughters()>0){
        for(uint32_t i=0; i<neutron->nDaughters(); ++i) {
          if(neutron->getDaughter(i)->atomicNumber()>0 && neutron->getDaughter(i)->pdgCode()!=2212){
          auto trk=neutron->getDaughter(i);
  if (trk->pdgName()=="Si28"){ Si28_once += 1;}
              else if (trk->pdgName()=="Si29"){
                Si29_once += 1;}
              else if (trk->pdgName()=="Si30"){
                Si30_once += 1;}
              else if (trk->pdgName()=="Si27"){
                Si27_once +=1 ;}
              else if (trk->pdgName()=="Si31"){
                Si31_once +=1 ;}
              else if (trk->pdgName()=="Al25"){
                Al25_once +=1 ;}
              else if (trk->pdgName()=="Al26"){
                Al26_once +=1 ;}
              else if (trk->pdgName()=="Al27"){
                Al27_once +=1 ;}
              else if (trk->pdgName()=="Al28"){
                Al28_once +=1 ;}
              else if (trk->pdgName()=="Al29"){
                Al29_once +=1 ;}
              else if (trk->pdgName()=="Al30"){
                Al30_once +=1 ;}
              else if (trk->pdgName()=="Mg24"){
                Mg24_once +=1 ;}
              else if (trk->pdgName()=="Mg25"){
                Mg25_once +=1 ;}
              else if (trk->pdgName()=="Mg26"){
                Mg26_once +=1 ;}
              else if (trk->pdgName()=="Mg27"){
                Mg27_once +=1 ;}
              else if (trk->pdgName()=="gamma"){
                gamma_once +=1 ;}
              else if (trk->pdgName()=="neutron"){
                neutron_once += 1;}
              else if (trk->pdgName()=="proton"){
                proton_once += 1;}
              else if (trk->pdgName()=="deuteron"){
                deuteron_once += 1;}
              else if (trk->pdgName()=="Na23"){
                Na23_once += 1;}
              else if (trk->pdgName()=="alpha"){
                alpha_once += 1;}
              else{
                // printf("Neutron daughter info %s %s %f keV atomic number %d \n",
          // trk->pdgNameCStr(),trk->creatorProcessCStr(),trk->startEKin()/Units::keV,trk->atomicNumber());
                // printf("Particle is %s /n",track->pdgNameCStr());
                others_once += 1;}
  if (trk->startEKin()*1000>0.021 && (trk->atomicNumber()==14 or trk->atomicNumber()==13 or trk->atomicNumber()==12)){
    h_recoil_PKA_high_Z ->fill(trk->startEKin()*1000); 
    // printf("Particle daughter info %s %s %f keV \n",
    // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()*1000);
  }else if(trk->startEKin()*1000>0.021 && (trk->atomicNumber()<8 and trk->atomicNumber()>0)){
    h_recoil_PKA_low_Z ->fill(trk->startEKin()*1000);}

  if(trk->startEKin()*1000>0.021 && trk->atomicNumber()>0){ h_recoil_PKA_once->fill(trk->startEKin()*1000); }

  break;
          }
      }
      }
    }
  }
  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  hc.saveToFile("silicon",true);

  return 0;
}


