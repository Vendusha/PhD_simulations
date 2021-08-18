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

  GriffDataReader dr(argc,argv);
  GriffAnaUtils::SegmentIterator diode(&dr);
  diode.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Diode"));


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
  ////////////////////////////////Filter primary recoils//////////////////////////////////////
  // GriffAnaUtils::TrackIterator primary_recoil(&dr);
  // primary_recoil.addFilter(new DaughterFilter);
  //////////////////////Filter primary particles////////////////////////////////////
  GriffAnaUtils::TrackIterator primary_particles(&dr);
  primary_particles.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto nbins = 8192;
  auto nbins_dim=1024;
  auto time_max=2;
  auto time_nbins=8;
  // auto time_bin=time_max/time_nbin;
  auto h_recoil_average = hc.book1D("Average_recoils",1300, 0, 1300, "Average_recoils");
  auto h_interaction_distance = hc.book1D("Interaction_distance",512, 0, 800, "Interaction_distance");
  auto h_distribution_time = hc.book1D("Distribution_in_time",time_nbins, 0, time_max, "Distribution_in_time");// this describes how many particles are in the system at each time 
  h_distribution_time -> setXLabel("Time [us]");
  h_distribution_time -> setYLabel("Particle count [-]");
  auto h_interaction_time = hc.book1D("Time_of_creation",time_nbins, 0, time_max, "Time_of_creation");// this describes how many particles are in the system at each time 
  h_interaction_time -> setXLabel("Time [us]");
  h_interaction_time -> setYLabel("Particle count [-]");
  auto h_interaction_distance_zoom = hc.book1D("Interaction_distance_zoom",512, 0,3, "Interaction_distance_zoom");
  auto h_interaction_distance_x = hc.book1D("Interaction_distance_x",512, 0,3, "Interaction_distance_x");
  auto h_interaction_distance_y = hc.book1D("Interaction_distance_y",512, 0,3, "Interaction_distance_y");
  auto h_interaction_distance_z = hc.book1D("Interaction_distance_z",512, 0,3, "Interaction_distance_z");
  h_interaction_distance_zoom -> setXLabel("Distance [um]");
  h_interaction_distance_zoom -> setYLabel("Reaction count [-]");
  h_interaction_distance_x -> setXLabel("Distance [um]");
  h_interaction_distance_x -> setYLabel("Reaction count [-]");
  h_interaction_distance_y -> setXLabel("Distance [um]");
  h_interaction_distance_y -> setXLabel("Distance [um]");
  h_interaction_distance_z -> setYLabel("Reaction count [-]");
  h_interaction_distance_z -> setYLabel("Reaction count [-]");
  h_interaction_distance -> setXLabel("Distance [um]");
  h_interaction_distance -> setYLabel("Reaction count [-]");
  auto h_cascade_atomic_number = hc.book1D("Atomic_number_concentration",32, 0, 32, "Atomic_number_concentration");
  auto h_PKA_atomic_number = hc.book1D("Atomic_number_concentration_PKA",32, 0, 32, "Atomic_number_concentration_PKA");
  auto h_recoil_average_all = hc.book1D("Average_silicon_recoils_no_gamma",1300, 0, 1300, "Average_silicon_recoils_no_gamma");
  h_recoil_average_all -> setXLabel("Number of produced silicon recoils");
  h_recoil_average_all -> setYLabel("Probability");
  auto h_recoil_average_silicon = hc.book1D("Average_silicon_recoils",1300, 0, 1300, "Average_silicon_recoils");
  h_recoil_average_silicon -> setXLabel("Number of produced silicon recoils");
  h_recoil_average_silicon -> setYLabel("Probability");
  auto h_primary_particle = hc.book1D("Primary_particles", nbins, energy_init*1/2*1000, energy_init*3/2*1000, "Primary_"+particle_type);
 h_primary_particle -> setXLabel("Energy [keV]");
 h_primary_particle -> setYLabel("Counts []");
 auto h_primary_recoil = hc.book1D("Recoil_spectra_PKA", nbins, 0, energy_init*1000, "Recoil_spectra_PKA");
 h_primary_recoil -> setXLabel("Energy [keV]");
 h_primary_recoil -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_NIEL = hc.book1D("PKA_NIEL", nbins, 0, 10, "PKA_NIEL");
 h_NIEL -> setXLabel("Energy [keV]");
 h_NIEL -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_no_thresh = hc.book1D("Recoil_spectra_PKA_no_thresh", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_no_thresh");
 h_primary_recoil_no_thresh -> setXLabel("Energy [keV]");
 h_primary_recoil_no_thresh -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Silicon = hc.book1D("Recoil_spectra_PKA_Silicon", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Silicon");
 h_primary_recoil_Silicon -> setXLabel("Energy [keV]");
 h_primary_recoil_Silicon -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_Silicon_no_thresh = hc.book1D("Recoil_spectra_PKA_Silicon_no_thresh", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_Silicon_no_thresh");
 h_primary_recoil_Silicon_no_thresh -> setXLabel("Energy [keV]");
 h_primary_recoil_Silicon_no_thresh -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_alpha = hc.book1D("Recoil_spectra_PKA_alpha", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_alpha");
 h_primary_recoil_alpha -> setXLabel("Energy [keV]");
 h_primary_recoil_alpha -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_alpha_no_thresh = hc.book1D("Recoil_spectra_PKA_alpha_no_thresh", nbins, 0, energy_init*1000, "Recoil_spectra_PKA_alpha_no_thresh");
 h_primary_recoil_alpha_no_thresh -> setXLabel("Energy [keV]");
 h_primary_recoil_alpha_no_thresh -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_primary_recoil_type = hc.bookCounts("PKA_type", "PKA_type");
 auto h_recoil_PKA_high_Z = hc.book1D("Recoil_spectra_high_Z", nbins,0, energy_init*1000, "Recoil_spectra_high_Z");
 h_recoil_PKA_high_Z -> setXLabel("Energy [keV]");
 h_recoil_PKA_high_Z -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_recoil_PKA_low_Z = hc.book1D("Recoil_spectra_low_Z", nbins, 0, energy_init*1000, "Recoil_spectra_low_Z");
 h_recoil_PKA_low_Z -> setXLabel("Energy [keV]");
 h_recoil_PKA_low_Z -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_recoil_Cascade = hc.book1D("Recoil_spectra_Cascade", nbins, 0, energy_init*1000, "Recoil_spectra_Cascade");
 auto h_Cascade_type = hc.bookCounts("Cascade_type", "Cascade_type");
 h_recoil_Cascade -> setXLabel("Energy [keV]");
 h_recoil_Cascade -> setYLabel("Frequency [PKA-1 keV-1]");
 auto h_particle_xy = hc.book2D("XY_distribution_neutrons", 8192,-0.7,0.7,8192,-0.7,0.7,"XY_distribution_neutrons"); //two dimensional histogram in the direction of the detector, visualise the incoming neutrons //values are in centimeters
 h_particle_xy->setXLabel("x_position [mm]");
 h_particle_xy->setYLabel("y_position [mm]");
 auto h_particle_yz = hc.book2D("YZ_distribution_neutrons", nbins_dim,-0.5,0.5,nbins_dim,-0.1,0.1,"YZ_distribution_neutrons"); //two dimensional histogram in the direction of the detector, visualise the incoming neutrons //values are in centimeters
  h_particle_yz->setXLabel("y_position [mm]");
  h_particle_yz->setYLabel("z_position [mm]");

  auto h_interaction_xz = hc.book2D("Point_of_interaction_xz", nbins_dim,-0.01,0.01,nbins_dim,-0.01,0.01,"Point_of_interaction_xz");
  auto h_interaction_xy = hc.book2D("Point_of_interaction_xy", nbins_dim,-0.5*Units::um,0.5*Units::um,nbins_dim,-0.5*Units::um,0.5*Units::um,"Point_of_interaction_xy");
  h_interaction_xz->setXLabel("x_position [um]");
  h_interaction_xz->setYLabel("z_position [um]");
  h_interaction_xy->setXLabel("x_position [um]");
  h_interaction_xy->setYLabel("y_position [um]");


  auto h_edep = hc.book1D("Energy_deposition",nbins, 0, energy_init*1000, "Energy_deposition");

  auto Si27 = h_primary_recoil_type ->addCounter("Si27");
  auto Si27_Cascade = h_Cascade_type ->addCounter("Si27");
  auto Si28 = h_primary_recoil_type ->addCounter("Si28");
  auto Si28_Cascade = h_Cascade_type ->addCounter("Si28");
  auto Si29 = h_primary_recoil_type ->addCounter("Si29");
  auto Si29_Cascade = h_Cascade_type ->addCounter("Si29");
  auto Si30 = h_primary_recoil_type ->addCounter("Si30");
  auto Si30_Cascade = h_Cascade_type ->addCounter("Si30");
  auto Si31 = h_primary_recoil_type ->addCounter("Si31");
  auto Si31_Cascade = h_Cascade_type ->addCounter("Si31");
  auto Al25 = h_primary_recoil_type ->addCounter("Al25");
  auto Al25_Cascade = h_Cascade_type ->addCounter("Al25");
  auto Al26 = h_primary_recoil_type ->addCounter("Al26");
  auto Al26_Cascade = h_Cascade_type ->addCounter("Al26");
  auto Al27 = h_primary_recoil_type ->addCounter("Al27");
  auto Al27_Cascade = h_Cascade_type ->addCounter("Al27");
  auto Al28 = h_primary_recoil_type ->addCounter("Al28");
  auto Al28_Cascade = h_Cascade_type ->addCounter("Al28");
  auto Al29 = h_primary_recoil_type ->addCounter("Al29");
  auto Al29_Cascade = h_Cascade_type ->addCounter("Al29");
  auto Al30 = h_primary_recoil_type ->addCounter("Al30");
  auto Al30_Cascade = h_Cascade_type ->addCounter("Al30");
  auto Mg24 = h_primary_recoil_type ->addCounter("Mg24");
  auto Mg24_Cascade = h_Cascade_type ->addCounter("Mg24");
  auto Mg25 = h_primary_recoil_type ->addCounter("Mg25");
  auto Mg25_Cascade = h_Cascade_type ->addCounter("Mg25");
  auto Mg26 = h_primary_recoil_type ->addCounter("Mg26");
  auto Mg26_Cascade = h_Cascade_type ->addCounter("Mg26");
  auto Mg27 = h_primary_recoil_type ->addCounter("Mg27");
  auto Mg27_Cascade = h_Cascade_type ->addCounter("Mg27");
  auto Na23 = h_primary_recoil_type ->addCounter("Na23");
  auto Na23_Cascade = h_Cascade_type ->addCounter("Na23");
  auto neutron = h_primary_recoil_type ->addCounter("neutron");
  auto neutron_Cascade = h_Cascade_type ->addCounter("neutron");
  auto proton = h_primary_recoil_type ->addCounter("proton");
  auto proton_Cascade = h_Cascade_type ->addCounter("proton");
  auto alpha = h_primary_recoil_type ->addCounter("alpha");
  auto alpha_Cascade = h_Cascade_type ->addCounter("alpha");
  auto gamma = h_primary_recoil_type ->addCounter("gamma");
  auto gamma_Cascade = h_Cascade_type ->addCounter("gamma");
  auto deuteron = h_primary_recoil_type ->addCounter("deuteron");
  auto deuteron_Cascade = h_Cascade_type ->addCounter("deuteron");
  auto others = h_primary_recoil_type ->addCounter("others");
  auto others_Cascade = h_Cascade_type ->addCounter("others");
  int first_daughter=1;

  // auto time_of_creation_parent = 0.0;
  //Loop over events and extract info via Griff interface:
  while (dr.loopEvents()) {
    // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
                                      // dr.runNumber(),dr.eventNumber(),dr.nTracks());
    auto edep = 0;
    while (auto segment = diode.next()) {
      edep+=segment->eDep();
    }
    if (edep>0){
      h_edep->fill(edep/Units::keV);}
 for (auto track = dr.trackBegin();track!=dr.trackEnd();++track) {
   if (track->isPrimary()){
     h_primary_particle ->fill(track->startEKin()*1000);
     //////////////////////////////Daughters section////////////////////////////
     if (track->nDaughters()>0){
       for(uint32_t i=0; i<track->nDaughters(); ++i) {
         auto primary_recoil = track->getDaughter(i);
         if (primary_recoil->startEKin()*1000>0.021 && (primary_recoil->atomicNumber()==14 or primary_recoil->atomicNumber()==13 or primary_recoil->atomicNumber()==12)){
           h_recoil_PKA_high_Z ->fill(primary_recoil->startEKin()*1000); }else if(primary_recoil->startEKin()*1000>0.021 && (primary_recoil->atomicNumber()<8 and primary_recoil->atomicNumber()>0)){
           h_recoil_PKA_low_Z ->fill(primary_recoil->startEKin()*1000);}
         if(primary_recoil->startEKin()*1000>0.021 && primary_recoil->atomicNumber()>1){h_primary_recoil->fill(primary_recoil->startEKin()*1000); }
         h_PKA_atomic_number->fill(primary_recoil->atomicNumber());
         if(primary_recoil->atomicNumber()>0){h_primary_recoil_no_thresh->fill(primary_recoil->startEKin()*1000); }
         h_PKA_atomic_number->fill(primary_recoil->atomicNumber());
         if(primary_recoil->startEKin()*1000>0.021 && primary_recoil->atomicNumber()==2){h_primary_recoil_alpha->fill(primary_recoil->startEKin()*1000); }
         if(primary_recoil->atomicNumber()==2){h_primary_recoil_alpha_no_thresh->fill(primary_recoil->startEKin()*1000); }
         // std::cout<<primary_recoil->atomicMass()<<std::endl;
         if(primary_recoil->startEKin()*1000>0.021 && primary_recoil->atomicNumber()==14){h_primary_recoil_Silicon->fill(primary_recoil->startEKin()*1000); }

         if(primary_recoil->atomicNumber()==14){h_primary_recoil_Silicon_no_thresh->fill(primary_recoil->startEKin()*1000); }

      if (primary_recoil->pdgName()=="Si28"){
                Si28+=1;
               }
              else if (primary_recoil->pdgName()=="Si29"){
                Si29+=1;}
              else if (primary_recoil->pdgName()=="Si30"){
                Si30+=1;}
              else if (primary_recoil->pdgName()=="Si27"){
                Si27+=1;}
              else if (primary_recoil->pdgName()=="Si31"){
                Si31+=1;}
              else if (primary_recoil->pdgName()=="Al25"){
                Al25+=1;}
              else if (primary_recoil->pdgName()=="Al26"){
                Al26+=1;}
              else if (primary_recoil->pdgName()=="Al27"){
                Al27+=1;}
              else if (primary_recoil->pdgName()=="Al28"){
                Al28+=1;}
              else if (primary_recoil->pdgName()=="Al29"){
                Al29+=1;}
              else if (primary_recoil->pdgName()=="Al30"){
                Al30+=1;}
              else if (primary_recoil->pdgName()=="Mg24"){
                Mg24+=1;}
              else if (primary_recoil->pdgName()=="Mg25"){
                Mg25+=1;}
              else if (primary_recoil->pdgName()=="Mg26"){
                Mg26+=1;}
              else if (primary_recoil->pdgName()=="Mg27"){
                Mg27+=1;}
              else if (primary_recoil->pdgName()=="gamma"){
                gamma+=1;
              }
              else if (primary_recoil->pdgName()=="neutron"){
                neutron+=1;
              }
              else if (primary_recoil->pdgName()=="proton"){
                proton+=1;}
              else if (primary_recoil->pdgName()=="deuteron"){
                deuteron+=1;}
              else if (primary_recoil->pdgName()=="Na23"){
                Na23+=1;}
              else if (primary_recoil->pdgName()=="alpha"){
                alpha+=1;}
              else{
                // printf("Neutron daughter info %s %s %f keV atomic number %d \n",
                       // primary_recoil->pdgNameCStr(),primary_recoil->creatorProcessCStr(),primary_recoil->startEKin()/Units::keV,primary_recoil->atomicNumber());
                // printf("Particle is %s /n",primary_recoil->pdgNameCStr());
                others+=1;}
       }
     }
       auto x_0 = 0.0;
       auto y_0 = 0.0;
       auto z_0 =0.0;
       auto time_0=0.0;
       first_daughter = 1;

     auto descendant_count = 0;
     auto descendant_count_no_gamma = 0;
     auto descendant_count_silicon = 0;
     auto edep_NIEL=0.0;
       descendantFilter->setAncestor(track);
       // printf("Primary particle. %s %s %f keV atomic number %d number of daughters %d \n",
       // track->pdgNameCStr(), track->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->nDaughters());
       // if (track->nDaughters()>0){
       while (auto particle = descendant.next()){
         if (!particle->isPrimary()){
           if (particle->getParent()->isPrimary()&& particle->nDaughters()>0 and first_daughter == 1){
           // new_event=1;
             time_0=particle->firstStep()->preTime()/Units::picosecond;
             x_0= particle->firstStep()->preGlobalX()/Units::um;
             y_0= particle->firstStep()->preGlobalY()/Units::um;
             z_0= particle->firstStep()->preGlobalZ()/Units::um;
             first_daughter = 0;
           }

         auto seg0 = particle->segmentBegin();
         auto segL = particle->segmentEnd();
         auto seg = seg0;
         for (;seg!=segL;++seg){
           if (seg->volumeName()=="Diode"){
             auto step0 = seg0->stepBegin();
             auto step = step0;
             auto stepE = seg0->stepEnd();
             for (;step!=stepE;++step){
               // edep_print = step->eDepNonIonising()/Units::eV;
               edep_NIEL = edep_NIEL+step->eDepNonIonising();
               // if (step->eDepNonIonising()>0){
               // std::cout<<edep_NIEL<<std::endl;
               // std::cout<<edep_print<<std::endl;
               }}}
         // h_NIEL->fill(edep_NIEL/Units::keV);
        

       // if (particle->firstStep()->preGlobalX()/Units::mm>0 && particle->firstStep()->preGlobalX()/Units::mm<0.001){
         h_particle_yz->fill(particle->firstStep()->preGlobalY()/Units::mm, particle->firstStep()->preGlobalZ()/Units::mm, particle->weight());
     // }
       if(!particle->isPrimary()){
         //////////////////////////The tree hierarchy print///////////////////////////
         //////////////////////////End of the tree hierarchy print
         auto time_of_creation=particle->firstStep()->preTime()/Units::picosecond-time_0;
         auto time_of_decay=particle->lastStep()->postTime()/Units::picosecond;
         // for(auto t = time_of_creation;t!=time_of_decay;t=t+time_bin){
           // h_distribution_time->fill(t);
         // }
         h_distribution_time->fill(time_of_decay);
         h_interaction_time->fill(time_of_creation);
       descendant_count+=1;
       auto x_now= (particle->firstStep()->preGlobalX()/Units::um)-x_0;
       auto y_now= (particle->firstStep()->preGlobalY()/Units::um)-y_0;
       auto z_now= (particle->firstStep()->preGlobalZ()/Units::um)-z_0;
       auto r = sqrt (pow(x_now,2)+pow(y_now,2)+pow(z_now,2));
       if (r>0){
         h_interaction_distance->fill(r);
         h_interaction_distance_zoom->fill(r);
         h_interaction_distance_x->fill(x_now);
         h_interaction_distance_y->fill(y_now);
         h_interaction_distance_z->fill(z_now);
         
         h_interaction_xz->fill(x_now,z_now);
         // if (z_now<=1){
           h_interaction_xy->fill(x_now,y_now);
       // }

       }
       h_cascade_atomic_number->fill(particle->atomicNumber());
       if (particle->atomicNumber()>1){descendant_count_no_gamma+=1;}
       if (particle->atomicNumber()==14){descendant_count_silicon+=1;}
       h_particle_xy->fill(particle->firstStep()->preGlobalX()/Units::mm, particle->firstStep()->preGlobalY()/Units::mm, particle->weight());

       // printf("Particle %s %s %f keV atomic number %d number of daughters %d \n",
              // particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->nDaughters());
     }
       }
     if (edep_NIEL>0){
       h_NIEL->fill(edep_NIEL/Units::keV);
         // std::cout<<edep_NIEL/Units::keV<<std::endl;
     }

     if (descendant_count>0){
       h_recoil_average ->fill(descendant_count);
       h_recoil_average_all ->fill(descendant_count_no_gamma);}
     if (descendant_count_silicon>0){
       h_recoil_average_silicon ->fill(descendant_count_silicon);}
     
   }
   }
   ///////////////////End of Primary particles//////////////////
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
   if (track->startEKin()*1000>0.021){
     h_recoil_Cascade ->fill(track->startEKin()*1000); 
   }
      }
 }
  }
  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  hc.saveToFile("silicon",true);

  return 0;
}


