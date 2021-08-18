#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

int main(int argc,char**argv) {

  //Open .griff file(s) specified on the command line:
  GriffDataReader dr(argc,argv);
  GriffAnaUtils::SegmentIterator silicon(&dr);
  silicon.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Diode"));

  //Extract and dump meta-data:
  dr.setup()->dump();

  //Check (as an example) that the simulated geometry is of the correct type:
  // if (dr.setup()->geo().getName()!="G4GeoSilicon/GeoSiliconLGAD")
    // return 1;
  //energy_init=0.027 keV for 0.2 MeV neutrons, 0.095 for 0.7 MeV neutrons, 2.861 for 3 MeV neutrons, 5.756 for 6 Mev neutrons, 9.7 for 10 MeV neutrons, 13.4 for 14 MeV neutrons, 17.4 for 18 MeV neutrons, 
  auto energy_init = 120;
  // auto energy_init= dr.setup()->gen().getParameterDouble("fixed_energy_eV")/1000000;
  // auto energy_init= dr.setup()->gen().getParameterDouble("energy_MeV");
  const std::string particle_type= dr.setup()->gen().getParameterString("particle_type");
  std::cout<< "The generator shooted " << particle_type <<" of the energy "<< energy_init << "GeV ."<<std::endl;
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto nbins = 4096;
  auto h_evtcounts = hc.bookCounts("Type of particle","evtcounts");
  auto h_ionization = hc.book1D("Ionization_spectra",nbins, 0,0.67, "Ionization_spectra");
  auto h_spallation = hc.book1D("Spallation_spectra",nbins, 0, 26, "Spallation_spectra");
  auto h_edep = hc.book1D("Edep", nbins, 0, 0.67, "Edep");
  auto h_edep_big = hc.book1D("Edep_everything", nbins, 0, 26, "Edep_everything");
  h_ionization -> setXLabel("Ionization energy [MeV]");
  h_ionization -> setYLabel("Counts [-]");
  h_edep -> setXLabel("Ionization energy [MeV]");
  h_edep -> setYLabel("Counts [-]");
  auto MeV_0_1 = h_evtcounts->addCounter("MoreThan0_1MeV");
  auto MeV_0_2 = h_evtcounts->addCounter("E_0_2MeV");
  auto MeV_0_5 = h_evtcounts->addCounter("E_0_5MeV");
  auto MeV_1 = h_evtcounts->addCounter("E1MeV");
  auto MeV_5 = h_evtcounts->addCounter("E5MeV");
  auto MeV_10 = h_evtcounts->addCounter("E10MeV");
  auto MeV_15 = h_evtcounts->addCounter("E15MeV");
  auto MeV_20 = h_evtcounts->addCounter("E20MeV");
  while (dr.loopEvents()) {
    bool Spallation = 0;
    // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
    // dr.runNumber(), dr.eventNumber(), dr.nTracks()); while (auto seg = si_transmission.next()) {
    double edep = 0;
    double w = 1;
    while (auto segment = silicon.next()){
      edep += segment->eDep();
      w = segment->getTrack()->weight();
    }
    if (edep>20){++MeV_20;}
    if (edep>15){++MeV_15;}
    if (edep>10){++MeV_10;}
    if (edep>5){++MeV_5;}
    if (edep>1){++MeV_1;}
    if (edep>0.5){++MeV_0_5;}
    if (edep>0.2){++MeV_0_2;}
    if (edep>0.1){++MeV_0_1;}
    if (edep>0){
      // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
      // dr.runNumber(), dr.eventNumber(), dr.nTracks());
      // printf("Crazy deposition now %f.\n", edep);
      for (auto track = dr.trackBegin();track!=dr.trackEnd();++track){
        if (track->isPrimary()){
          printf("Primary particle info %s %s %f keV \n",
            track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()*1000);
          for(uint i=0; i<track->nDaughters(); ++i){
            auto trk = track->getDaughter(i);
    printf(" ===info %s %s %f keV \n",
           trk -> pdgNameCStr(),trk -> creatorProcessCStr(),trk -> startEKin()*1000);
            if (trk-> creatorProcess()=="protonInelastic"){Spallation = 1;}
      }

    if (edep>0){h_edep->fill(edep/Units::MeV,w);
      h_edep_big->fill(edep/Units::MeV,w);
      if (Spallation == 1){
        h_spallation->fill(edep/Units::MeV,w);
      }else{
        h_ionization->fill(edep/Units::MeV,w);
      }
    }

      // printf("info %s %s %f GeV atomic number %d \n",
             // trk->pdgNameCStr(),trk->creatorProcessCStr(),trk->startEKin()/Units::GeV,trk->atomicNumber());
        }
      }
      // printf("Particle is %s /n",track->pdgNameCStr());
    // printf("Particle daughter info %s %s %f keV \n",
    // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()*1000);
      }
  }
  hc.saveToFile("silicon",true);

  return 0;
}


