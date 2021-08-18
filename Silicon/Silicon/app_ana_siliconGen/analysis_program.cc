#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>


class DaughterFilter: public GriffAnaUtils::ITrackFilter {

public:
  bool filter(const GriffDataRead::Track*trk) const {

    auto par = trk->getParent();
    if (!par) return false;
    if (!par->isPrimary()) return false;
    //11 is an electron, one could alternatively
    //check if par->pdgName() equals "electron":
    // return par->pdgName() == particle_type;
    return par->pdgCode() == 1000140280;
  }
};
class SKAFilter: public GriffAnaUtils::ITrackFilter {

public:
  bool filter(const GriffDataRead::Track*trk) const {

    auto par = trk->getParent();
    if (!par) return false;
    // if (par->atomicNumber!=14) return false;
    auto grandpar = par->getParent();
    if (!grandpar) return false;
    if (!grandpar->isPrimary()) return false;
    //11 is an electron, one could alternatively
    //check if par->pdgName() equals "electron":
    // return par->pdgName() == particle_type;
    return par->atomicNumber() == 14;
  }
};

int main(int argc,char**argv) {

  //Open .griff file(s) specified on the command line:
  GriffDataReader dr(argc,argv);
  // GriffAnaUtils::SegmentIterator diode(&diode)

  //Extract and dump meta-data:
  dr.setup()->dump();

  //Check (as an example) that the simulated geometry is of the correct type:
  if (dr.setup()->geo().getName()!="G4GeoSilicon/GeoSilicon")
    return 1;
  // auto energy_init = 200;
  auto energy_init= dr.setup()->gen().getParameterDouble("fixed_energy_eV");
  // const std::string particle_type= dr.setup()->gen().getParameterString("particle_type");
  std::cout<< "The generator shooted Silicon particles of the energy "<< energy_init << "eV."<<std::endl;

  GriffAnaUtils::TrackIterator ti(&dr);
  ti.addFilter(new DaughterFilter);
  GriffAnaUtils::TrackIterator SKA(&dr);
  SKA.addFilter(new SKAFilter);
  GriffAnaUtils::SegmentIterator silicon(&dr);
  silicon.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Silicon"));


  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto nbins = 8192;
  auto nbins_dim=1024;
  // auto nbins_dim=128;
  // auto n = 0;
  auto h_primary = hc.book1D("Primary_silicon",nbins, 0, energy_init/1000, "Primary_silicon");
  auto h_recoil_PKA_high_Z = hc.book1D("Recoil_spectra_high_Z",nbins, 0, energy_init/1000, "Recoil_spectra_high_Z");
  auto h_recoil_PKA = hc.book1D("Recoil_spectra_all_PKA",nbins, 0, energy_init/1000, "Recoil_spectra_all_PKA");
  auto h_recoil_PKA_low_Z = hc.book1D("Recoil_spectra_low_Z",nbins, 0, energy_init/1000, "Recoil_spectra_low_Z");
  auto h_recoil_SKA = hc.book1D("Recoil_spectra_SKA",nbins, 0, energy_init/1000, "Recoil_spectra_SKA");
  auto h_recoil_Cascade = hc.book1D("Recoil_spectra_Cascade",nbins, 0, energy_init/1000, "Recoil_spectra_Cascade");
  // auto h_recoil_SKA = hc.book1D("Recoil_spectra_SKA",nbins, 0,energy_init*1000, "Recoil_spectra_SKA");
  auto h_Daughter_PKA = hc.bookCounts("PKA_particle_type", "PKA_particle_type");
  auto h_Daughter_SKA = hc.bookCounts("SKA_particle_type", "SKA_particle_type");
  auto h_Daughter_Cascade = hc.bookCounts("Cascade_particle_type", "Cascade_particle_type");
  auto h_particle_yz = hc.book2D("YZ_distribution_neutrons",nbins_dim,-10,10,nbins_dim,-10,10,"YZ_distribution_neutrons"); //two dimensional histogram in the direction of the detector, visualise the incoming neutrons //values are in centimeters
  h_particle_yz->setXLabel("y_position [cm]");
  h_particle_yz->setYLabel("z_position [cm]");

  auto Si27 = h_Daughter_PKA ->addCounter("Si27");
  auto Si27_SKA = h_Daughter_SKA ->addCounter("Si27");
  auto Si27_Cascade = h_Daughter_Cascade ->addCounter("Si27");
  auto Si28 = h_Daughter_PKA ->addCounter("Si28");
  auto Si28_SKA = h_Daughter_SKA ->addCounter("Si28");
  auto Si28_Cascade = h_Daughter_Cascade ->addCounter("Si28");
  auto Si29 = h_Daughter_PKA ->addCounter("Si29");
  auto Si29_SKA = h_Daughter_SKA ->addCounter("Si29");
  auto Si29_Cascade = h_Daughter_Cascade ->addCounter("Si29");
  auto Si30 = h_Daughter_PKA ->addCounter("Si30");
  auto Si30_SKA = h_Daughter_SKA ->addCounter("Si30");
  auto Si30_Cascade = h_Daughter_Cascade ->addCounter("Si30");
  auto Si31 = h_Daughter_PKA ->addCounter("Si31");
  auto Si31_SKA = h_Daughter_SKA ->addCounter("Si31");
  auto Si31_Cascade = h_Daughter_Cascade ->addCounter("Si31");
  auto Al25 = h_Daughter_PKA ->addCounter("Al25");
  auto Al25_SKA = h_Daughter_SKA ->addCounter("Al25");
  auto Al25_Cascade = h_Daughter_Cascade ->addCounter("Al25");
  auto Al26 = h_Daughter_PKA ->addCounter("Al26");
  auto Al26_SKA = h_Daughter_SKA ->addCounter("Al26");
  auto Al26_Cascade = h_Daughter_Cascade ->addCounter("Al26");
  auto Al27 = h_Daughter_PKA ->addCounter("Al27");
  auto Al27_SKA = h_Daughter_SKA ->addCounter("Al27");
  auto Al27_Cascade = h_Daughter_Cascade ->addCounter("Al27");
  auto Al28 = h_Daughter_PKA ->addCounter("Al28");
  auto Al28_SKA = h_Daughter_SKA ->addCounter("Al28");
  auto Al28_Cascade = h_Daughter_Cascade ->addCounter("Al28");
  auto Al29 = h_Daughter_PKA ->addCounter("Al29");
  auto Al29_SKA = h_Daughter_SKA ->addCounter("Al29");
  auto Al29_Cascade = h_Daughter_Cascade ->addCounter("Al29");
  auto Al30 = h_Daughter_PKA ->addCounter("Al30");
  auto Al30_SKA = h_Daughter_SKA ->addCounter("Al30");
  auto Al30_Cascade = h_Daughter_Cascade ->addCounter("Al30");
  auto Mg24 = h_Daughter_PKA ->addCounter("Mg24");
  auto Mg24_SKA = h_Daughter_SKA ->addCounter("Mg24");
  auto Mg24_Cascade = h_Daughter_Cascade ->addCounter("Mg24");
  auto Mg25 = h_Daughter_PKA ->addCounter("Mg25");
  auto Mg25_SKA = h_Daughter_SKA ->addCounter("Mg25");
  auto Mg25_Cascade = h_Daughter_Cascade ->addCounter("Mg25");
  auto Mg26 = h_Daughter_PKA ->addCounter("Mg26");
  auto Mg26_SKA = h_Daughter_SKA ->addCounter("Mg26");
  auto Mg26_Cascade = h_Daughter_Cascade ->addCounter("Mg26");
  auto Mg27 = h_Daughter_PKA ->addCounter("Mg27");
  auto Mg27_SKA = h_Daughter_SKA ->addCounter("Mg27");
  auto Mg27_Cascade = h_Daughter_Cascade ->addCounter("Mg27");
  auto Na23 = h_Daughter_PKA ->addCounter("Na23");
  auto Na23_SKA = h_Daughter_SKA ->addCounter("Na23");
  auto Na23_Cascade = h_Daughter_Cascade ->addCounter("Na23");
  auto neutron = h_Daughter_PKA ->addCounter("neutron");
  auto neutron_SKA = h_Daughter_SKA ->addCounter("neutron");
  auto neutron_Cascade = h_Daughter_Cascade ->addCounter("neutron");
  auto proton = h_Daughter_PKA ->addCounter("proton");
  auto proton_SKA = h_Daughter_SKA ->addCounter("proton");
  auto proton_Cascade = h_Daughter_Cascade ->addCounter("proton");
  auto alpha = h_Daughter_PKA ->addCounter("alpha");
  auto alpha_SKA = h_Daughter_SKA ->addCounter("alpha");
  auto alpha_Cascade = h_Daughter_Cascade ->addCounter("alpha");
  auto gamma = h_Daughter_PKA ->addCounter("gamma");
  auto gamma_SKA = h_Daughter_SKA ->addCounter("gamma");
  auto gamma_Cascade = h_Daughter_Cascade ->addCounter("gamma");
  auto deuteron = h_Daughter_PKA ->addCounter("deuteron");
  auto deuteron_SKA = h_Daughter_SKA ->addCounter("deuteron");
  auto deuteron_Cascade = h_Daughter_Cascade ->addCounter("deuteron");
  auto others = h_Daughter_PKA ->addCounter("others");
  auto others_SKA = h_Daughter_SKA ->addCounter("others");
  auto others_Cascade = h_Daughter_Cascade ->addCounter("others");
  h_recoil_PKA_high_Z -> setXLabel("Energy [keV]");
  h_recoil_PKA_high_Z -> setYLabel("Frequency [PKA-1 keV-1]");
  h_recoil_PKA_low_Z -> setXLabel("Energy [keV]");
  h_recoil_PKA_low_Z -> setYLabel("Frequency [PKA-1 keV-1]");
  h_recoil_SKA -> setXLabel("Energy [keV]");
  h_recoil_SKA -> setYLabel("Frequency [SKA-1 keV-1]");
  //Loop over events and extract info via Griff interface:
  while (dr.loopEvents()) {
    // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
                                      // dr.runNumber(),dr.eventNumber(),dr.nTracks());
    // double edep = 0;
    // while (auto segment = silicon.next()){
      // edep += segment->eDep();
    // }
    // printf("This particle deposited this energy. %f", edep);
    while (auto track = ti.next()){
      // printf("Particle info %s %s %f keV parent is %s \n",
             // track->pdgNameCStr(), track->creatorProcessCStr(), track->startEKin()/Units::keV, track->getParent()->pdgNameCStr());
      if (track->pdgName()=="Si28"){
                Si28+=1;
                // printf("Silicon Knocked Out! Number of babies %d PDGname is %s and PDGcode is %i. \n",track->nDaughters(), track->pdgNameCStr(), track->pdgCode());
                //             // printf("Silicon Knocked Out! Number of babies is %d and it happened at x %f y %f z %f Parent is %s.\n",track->nDaughters(), track->firstStep()->preGlobalX()/Units::cm, track->firstStep()->preGlobalY()/Units::cm, track->firstStep()->preGlobalZ()/Units::cm, track->getParent()->pdgNameCStr());
               }
              else if (track->pdgName()=="Si29"){
                Si29+=1;}
              else if (track->pdgName()=="Si30"){
                Si30+=1;}
              else if (track->pdgName()=="Si27"){
                Si27+=1;}
              else if (track->pdgName()=="Si31"){
                Si31+=1;}
              else if (track->pdgName()=="Al25"){
                Al25+=1;}
              else if (track->pdgName()=="Al26"){
                Al26+=1;}
              else if (track->pdgName()=="Al27"){
                Al27+=1;}
              else if (track->pdgName()=="Al28"){
                Al28+=1;}
              else if (track->pdgName()=="Al29"){
                Al29+=1;}
              else if (track->pdgName()=="Al30"){
                Al30+=1;}
              else if (track->pdgName()=="Mg24"){
                Mg24+=1;}
              else if (track->pdgName()=="Mg25"){
                Mg25+=1;}
              else if (track->pdgName()=="Mg26"){
                Mg26+=1;}
              else if (track->pdgName()=="Mg27"){
                Mg27+=1;}
              else if (track->pdgName()=="gamma"){
                gamma+=1;
              }
              else if (track->pdgName()=="neutron"){
                neutron+=1;
                // printf("Particle info %s %s %f keV parent is %s \n",
                       // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()/Units::keV,track->getParent()->pdgNameCStr());
              }
              else if (track->pdgName()=="proton"){
                proton+=1;}
              else if (track->pdgName()=="deuteron"){
                deuteron+=1;}
              else if (track->pdgName()=="Na23"){
                Na23+=1;}
              else if (track->pdgName()=="alpha"){
                alpha+=1;}
              else{
                // printf("Neutron daughter info %s %s %f keV atomic number %d \n",
                       // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()/Units::keV,track->atomicNumber());
                // printf("Particle is %s /n",track->pdgNameCStr());
                others+=1;}
      if (track -> atomicNumber()==14 or track -> atomicNumber()==13 or track -> atomicNumber()==12){
        h_recoil_PKA_high_Z ->fill(track->startEKin()/Units::keV); 
        // printf("Particle daughter info %s %s %f keV \n",
        // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()*1000);
      }else if(track->atomicNumber()<8 and track->atomicNumber()>0){
        h_recoil_PKA_low_Z ->fill(track->startEKin()/Units::keV);}
      if(track->atomicNumber()>0){ h_recoil_PKA->fill(track->startEKin()/Units::keV);
      }}
    //////////////////////////LOOKING AT THE SKA SPECTRA////////////////////////
 while (auto track = SKA.next()){
   if (track->pdgName()=="Si28"){
                Si28_SKA+=1;
               }
              else if (track->pdgName()=="Si29"){
                Si29_SKA+=1;}
              else if (track->pdgName()=="Si30"){
                Si30_SKA+=1;}
              else if (track->pdgName()=="Si27"){
                Si27_SKA+=1;}
              else if (track->pdgName()=="Si31"){
                Si31_SKA+=1;}
              else if (track->pdgName()=="Al25"){
                Al25_SKA+=1;}
              else if (track->pdgName()=="Al26"){
                Al26_SKA+=1;}
              else if (track->pdgName()=="Al27"){
                Al27_SKA+=1;}
              else if (track->pdgName()=="Al28"){
                Al28_SKA+=1;}
              else if (track->pdgName()=="Al29"){
                Al29_SKA+=1;}
              else if (track->pdgName()=="Al30"){
                Al30_SKA+=1;}
              else if (track->pdgName()=="Mg24"){
                Mg24_SKA+=1;}
              else if (track->pdgName()=="Mg25"){
                Mg25_SKA+=1;}
              else if (track->pdgName()=="Mg26"){
                Mg26_SKA+=1;}
              else if (track->pdgName()=="Mg27"){
                Mg27_SKA+=1;}
              else if (track->pdgName()=="gamma"){
                gamma_SKA+=1;}
              else if (track->pdgName()=="neutron"){
                neutron_SKA+=1;}
              else if (track->pdgName()=="proton"){
                proton_SKA+=1;}
              else if (track->pdgName()=="deuteron"){
                deuteron_SKA+=1;}
              else if (track->pdgName()=="Na23"){
                Na23_SKA+=1;}
              else if (track->pdgName()=="alpha"){
                alpha_SKA+=1;}
              else{
                // printf("Neutron daughter info %s %s %f keV atomic number %d \n",
                       // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()/Units::keV,track->atomicNumber());
                // printf("Particle is %s /n",track->pdgNameCStr());
                others_SKA+=1;}
        h_recoil_SKA ->fill(track->startEKin()*1000); 
    }

 for (auto track = dr.trackBegin();track!=dr.trackEnd();++track) {
   if (track->isPrimary()){
     h_primary->fill(track->startEKin()*1000);
   }
   if (!track->isPrimary()){
   if (track->pdgName()=="Si28"){
                Si28_Cascade+=1;
               }
              else if (track->pdgName()=="Si29"){
                Si29_Cascade+=1;}
              else if (track->pdgName()=="Si30"){
                Si30_Cascade+=1;}
              else if (track->pdgName()=="Si27"){
                Si27_Cascade+=1;}
              else if (track->pdgName()=="Si31"){
                Si31_Cascade+=1;}
              else if (track->pdgName()=="Al25"){
                Al25_Cascade+=1;}
              else if (track->pdgName()=="Al26"){
                Al26_Cascade+=1;}
              else if (track->pdgName()=="Al27"){
                Al27_Cascade+=1;}
              else if (track->pdgName()=="Al28"){
                Al28_Cascade+=1;}
              else if (track->pdgName()=="Al29"){
                Al29_Cascade+=1;}
              else if (track->pdgName()=="Al30"){
                Al30_Cascade+=1;}
              else if (track->pdgName()=="Mg24"){
                Mg24_Cascade+=1;}
              else if (track->pdgName()=="Mg25"){
                Mg25_Cascade+=1;}
              else if (track->pdgName()=="Mg26"){
                Mg26_Cascade+=1;}
              else if (track->pdgName()=="Mg27"){
                Mg27_Cascade+=1;}
              else if (track->pdgName()=="gamma"){
                gamma_Cascade+=1;}
              else if (track->pdgName()=="neutron"){
                neutron_Cascade+=1;}
              else if (track->pdgName()=="proton"){
                proton_Cascade+=1;}
              else if (track->pdgName()=="deuteron"){
                deuteron_Cascade+=1;}
              else if (track->pdgName()=="Na23"){
                Na23_Cascade+=1;}
              else if (track->pdgName()=="alpha"){
                alpha_Cascade+=1;}
              else{
                // printf("Neutron daughter info %s %s %f keV atomic number %d \n",
                       // track->pdgNameCStr(),track->creatorProcessCStr(),track->startEKin()/Units::keV,track->atomicNumber());
                // printf("Particle is %s /n",track->pdgNameCStr());
                others_Cascade+=1;}
   h_recoil_Cascade ->fill(track->startEKin()/Units::keV); 
   if (track->atomicNumber()==14){
     h_particle_yz->fill(track->firstStep()->preGlobalY()/Units::um, track->firstStep()->preGlobalZ()/Units::um, track->weight());}
   }
 }
    }
    //       // std::cout<<trk->getDaughter()->pdgNameCstr()
    //     }
        // if (seg->volumeName()=="Diode" && seg->endAtVolumeBoundary()) {
          // auto s = seg->firstStep();
          // double x(s->preLocalX()/Units::mm),y(s->preLocalY()/Units::mm);
          // if (x*x+y*y>0.01) {
            // h2d_box_hitmap->fill(x,y);
            // h1d_box_hitradius->fill(sqrt(x*x+y*y));
          // }
        // }
      // }
    // }
  // }

  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  hc.saveToFile("silicon",true);

  return 0;
}


