#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>
#include <cstdlib>
#include <math.h>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

class DaughterFilter: public GriffAnaUtils::ITrackFilter {

public:
  bool filter(const GriffDataRead::Track*trk) const {

    auto par = trk->getParent();
    if (!par) return false;
    return par->isPrimary();
  }
};

int main(int argc,char**argv) {
  std::fstream Coulomb;
  Coulomb.open("Coulomb.txt", std::ofstream::out | std::ofstream::trunc);
  Coulomb.close();
  std::fstream Elastic;
  Elastic.open("Elastic.txt", std::ofstream::out | std::ofstream::trunc);
  Elastic.close();
  std::fstream Inelastic;
  Inelastic.open("Inelastic.txt", std::ofstream::out | std::ofstream::trunc);
  Inelastic.close();
  GriffDataReader dr(argc,argv);
  GriffAnaUtils::SegmentIterator diode(&dr);
  diode.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Diode"));

  // std::fstream Coulomb;
  Coulomb.open ("Coulomb.txt",std::fstream::app);
  // std::fstream Elastic;
  Elastic.open ("Elastic.txt",std::fstream::app);
  // std::fstream Inelastic;
  Inelastic.open ("Inelastic.txt",std::fstream::app);

  dr.setup()->dump();

  if (dr.setup()->geo().getName()!="G4GeoSilicon/GeoSilicon")
    return 1;
  ///////////////////Printing info//////////////////////
  // auto energy_init = 200;
  auto energy_init= dr.setup()->gen().getParameterDouble("energy_MeV");
  const std::string particle_type= dr.setup()->gen().getParameterString("particle_type");
  std::cout<< "The generator shooted " << particle_type <<" of the energy "<< energy_init << "MeV ."<<std::endl;
  // GriffAnaUtils::TrackIterator neutrons(&dr); 
  // neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112)); //PDGCode for neutrons
  auto descendantFilter = new GriffAnaUtils::TrackFilter_Descendant();
   // descendantFilter->setNegated();
  GriffAnaUtils::TrackIterator descendant(&dr);
  descendant.addFilter(descendantFilter);
  //////////////////////Filter primary particles////////////////////////////////////
  GriffAnaUtils::TrackIterator primary_particles(&dr);
  primary_particles.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto nbins = 8192;
  //hadElastic, protonInelastic, CoulombScat
  ////////////////Atomic number/////////////////////
  auto h_PKA_atomic_number = hc.book1D("Atomic_number_PKA",32, 0, 32, "Atomic_number__PKA");
  h_PKA_atomic_number ->setXLabel("Proton number Z");
  h_PKA_atomic_number ->setYLabel("Number of particles created");
  auto h_PKA_atomic_number_Coulomb = hc.book1D("Atomic_number_PKA_Coulomb",32, 0, 32, "Atomic_number_PKA_Coulomb");
  h_PKA_atomic_number_Coulomb ->setXLabel("Proton number Z");
  h_PKA_atomic_number_Coulomb ->setYLabel("Number of particles created");
  auto h_PKA_atomic_number_Elastic = hc.book1D("Atomic_number_PKA_Elastic",32, 0, 32, "Atomic_number_PKA_Elastic");
  h_PKA_atomic_number_Elastic ->setXLabel("Proton number Z");
  h_PKA_atomic_number_Elastic ->setYLabel("Number of particles created");
    auto h_PKA_atomic_number_Inelastic = hc.book1D("Atomic_number_PKA_Inelastic",32, 0, 32, "Atomic_number_PKA_Inelastic");
    h_PKA_atomic_number_Inelastic ->setXLabel("Proton number Z");
    h_PKA_atomic_number_Inelastic ->setYLabel("Number of particles created");
  ////////////////Recoils/////////////////////
 auto h_primary_recoil = hc.book1D("Recoil_spectra_PKA", nbins, 0, energy_init*1000, "Recoil_spectra_PKA");
 h_primary_recoil -> setXLabel("Energy [keV]");
 h_primary_recoil -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Coulomb = hc.book1D("Recoil_spectra_PKA_Coulomb", nbins, 0, energy_init*1000, "Recoil_spectra_Coulomb");
 h_primary_recoil_Coulomb -> setXLabel("Energy [keV]");
 h_primary_recoil_Coulomb -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Elastic = hc.book1D("Recoil_spectra_PKA_Elastic", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Elastic");
 h_primary_recoil_Elastic -> setXLabel("Energy [keV]");
 h_primary_recoil_Elastic -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Inelastic = hc.book1D("Recoil_spectra_PKA_Inelastic", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Inelastic");
 h_primary_recoil_Inelastic -> setXLabel("Energy [keV]");
 h_primary_recoil_Inelastic -> setYLabel("Frequency [PKA-1 keV-1]");
 ////////////////Recoils Silicon/////////////////////
 auto h_primary_recoil_Silicon = hc.book1D("Recoil_spectra_PKA_Silicon", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Silicon");
 h_primary_recoil_Silicon -> setXLabel("Energy [keV]");
 h_primary_recoil_Silicon -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Silicon_Coulomb = hc.book1D("Recoil_spectra_PKA_Silicon_Coulomb", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Silicon_Coulomb");
 h_primary_recoil_Silicon_Coulomb -> setXLabel("Energy [keV]");
 h_primary_recoil_Silicon_Coulomb -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Silicon_Elastic = hc.book1D("Recoil_spectra_PKA_Silicon_Elastic", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Silicon_Elastic");
 h_primary_recoil_Silicon_Elastic -> setXLabel("Energy [keV]");
 h_primary_recoil_Silicon_Elastic -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Silicon_Inelastic = hc.book1D("Recoil_spectra_PKA_Silicon_Inelastic", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Silicon_Inelastic");
 h_primary_recoil_Silicon_Inelastic -> setXLabel("Energy [keV]");
 h_primary_recoil_Silicon_Inelastic -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_alpha = hc.book1D("Recoil_spectra_PKA_alpha", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_alpha");
 h_primary_recoil_alpha -> setXLabel("Energy [keV]");
 h_primary_recoil_alpha -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_recoil_PKA_high_Z = hc.book1D("Recoil_spectra_high_Z", nbins,0, energy_init*1000, "Recoil_spectra_high_Z");
 h_recoil_PKA_high_Z -> setXLabel("Energy [keV]");
 h_recoil_PKA_high_Z -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_recoil_PKA_low_Z = hc.book1D("Recoil_spectra_low_Z", nbins, 0, energy_init*1000, "Recoil_spectra_low_Z");
 h_recoil_PKA_low_Z -> setXLabel("Energy [keV]");
 h_recoil_PKA_low_Z -> setYLabel("Frequency [PKA-1 keV-1]");
 auto new_Coulomb = 0;
 auto new_Elastic = 0;
 auto new_Inelastic = 0;
 unsigned int no_track = 0;
  while (dr.loopEvents()) {
    // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
                                      // dr.runNumber(),dr.eventNumber(),dr.nTracks());
    new_Coulomb = 1;
    new_Elastic = 1;
    new_Inelastic = 1;
    no_track=0;
 for (auto track = dr.trackBegin();track!=dr.trackEnd();++track) {
            no_track += 1;
   // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
          // dr.runNumber(),dr.eventNumber(),dr.nTracks());
     //////////////////////////////Daughters section////////////////////////////
         auto primary_recoil = track;
         if (primary_recoil->atomicNumber()==14 or primary_recoil->atomicNumber()==13 or primary_recoil->atomicNumber()==12){
           h_recoil_PKA_high_Z ->fill(primary_recoil->startEKin()*1000); }else if(primary_recoil->startEKin()*1000>0.021 && (primary_recoil->atomicNumber()<8 and primary_recoil->atomicNumber()>0)){
           h_recoil_PKA_low_Z ->fill(primary_recoil->startEKin()*1000);}
         if(primary_recoil->atomicNumber()>1){h_primary_recoil->fill(primary_recoil->startEKin()*1000); }
         h_PKA_atomic_number->fill(primary_recoil->atomicNumber());
         if(primary_recoil->atomicNumber()==2){h_primary_recoil_alpha->fill(primary_recoil->startEKin()*1000); }
         // std::cout<<primary_recoil->atomicMass()<<std::endl;
         if(primary_recoil->atomicNumber()==14){h_primary_recoil_Silicon->fill(primary_recoil->startEKin()*1000); }
         //////////division into processes///////////
         if(primary_recoil -> creatorProcess() == "CoulombScat"){
           h_PKA_atomic_number_Coulomb->fill( primary_recoil -> atomicNumber());
          h_primary_recoil_Coulomb->fill( primary_recoil -> startEKin()*1000);
          h_primary_recoil_Silicon_Coulomb->fill( primary_recoil -> startEKin()*1000);
         h_primary_recoil->fill( primary_recoil -> startEKin()*1000);
         if( new_Coulomb == 1 && new_Elastic == 1 && new_Inelastic==1 )
           { new_Coulomb = 0; }
         }else if (primary_recoil -> creatorProcess() == "protonInelastic" || primary_recoil -> creatorProcess() == "neutronInelastic"){
           h_PKA_atomic_number_Inelastic -> fill(primary_recoil->atomicNumber());
           if( primary_recoil->atomicNumber()>1 ){h_primary_recoil_Inelastic->fill(primary_recoil->startEKin()*1000); }
           if( primary_recoil->atomicNumber()>1 ){h_primary_recoil_Silicon_Inelastic->fill(primary_recoil->startEKin()*1000); }
           if(new_Inelastic==1){
             new_Inelastic =0;}
         }else if (primary_recoil->creatorProcess()=="hadElastic"){
           h_PKA_atomic_number_Elastic->fill(primary_recoil->atomicNumber());
           h_primary_recoil_Elastic->fill(primary_recoil->startEKin()*1000);
           h_primary_recoil_Silicon_Elastic->fill(primary_recoil->startEKin()*1000);
           if (new_Elastic==1){
             new_Elastic=0;}
        }
         if (no_track==dr.nTracks()){
           if(new_Elastic==0){Elastic<<dr.eventNumber()<<std::endl;
           }else if (new_Inelastic ==0){
             Inelastic<<dr.eventNumber()<<std::endl;}
         else{
           Coulomb<<dr.eventNumber()<<std::endl;}}
 }
  }
  Coulomb.close();
  Inelastic.close();
  Elastic.close();
   
  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  // hc.saveToFile("PKA",true);

  return 0;
}


