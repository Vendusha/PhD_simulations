#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

int main(int argc,char**argv) {
  
  const char *arr[]={"0_1.txt","0_2.txt", "0_3.txt", "0_4.txt", "0_5.txt", "0_6.txt", "0_7.txt", "0_8.txt", "0_9.txt", "1_0.txt", "1_1.txt", "1_2.txt", "1_3.txt", "1_4.txt", "1_5.txt"};
  const char *arr_all[]={"0_1_all.txt","0_2_all.txt", "0_3_all.txt", "0_4_all.txt", "0_5_all.txt", "0_6_all.txt", "0_7_all.txt", "0_8_all.txt", "0_9_all.txt", "1_0_all.txt", "1_1_all.txt", "1_2_all.txt", "1_3_all.txt", "1_4_all.txt","1_5_all.txt"};
  std::ofstream ofs;
    for (int i=0; i<15;i++){
  ofs.open(arr[i], std::ofstream::out | std::ofstream::trunc);
  ofs.close();
    ofs.open(arr_all[i], std::ofstream::out | std::ofstream::trunc);
    ofs.close();}
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
  auto descendantFilter = new GriffAnaUtils::TrackFilter_Descendant();
   // descendantFilter->setNegated();
  GriffAnaUtils::TrackIterator descendant(&dr);
  descendant.addFilter(descendantFilter);
  //////////////////////Filter primary particles////////////////////////////////////
  GriffAnaUtils::TrackIterator primary_particles(&dr);
  primary_particles.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto new_event = 0;
  int event_no=1;
  int first_daughter=1;
  //Loop over events and extract info via Griff interface:
  while (dr.loopEvents()) {
 for (auto track = dr.trackBegin();track!=dr.trackEnd();++track) {
     auto descendant_count = 0;
     // auto x =track->firstStep()->preGlobalX();
     // auto y =track->firstStep()->preGlobalY();

     if (track->isPrimary()){
       auto x_0 = 0.0;
       auto y_0 = 0.0;
       auto z_0 =0.0;
       auto time_0=0.0;
       first_daughter = 1;

           if (track->nDaughters()>0){
             if (new_event !=0){
               ++event_no;
             for (int i=0; i<15;i++){
               std::fstream myfile;
               myfile.open (arr[i],std::fstream::app);
               myfile <<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<event_no<<std::endl;
               myfile.close();
               myfile.open (arr_all[i],std::fstream::app);
               myfile <<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<event_no<<std::endl;
               myfile.close();}
             new_event=0;}

             descendantFilter->setAncestor(track);

             // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
             //                                                            dr.runNumber(),dr.eventNumber(),dr.nTracks());
             // printf("Primary particle. %s %s %f keV atomic number %d number of daughters %d \n",
             //        track->pdgNameCStr(), track->creatorProcessCStr(), track->startEKin()/Units::keV, track->atomicNumber(), track->nDaughters());

           while (auto particle = descendant.next()){
             if (!particle->isPrimary()){
               if (particle->getParent()->isPrimary()&& particle->nDaughters()>0 and first_daughter == 1){
                 new_event=1;
                 time_0=particle->firstStep()->preTime()/Units::picosecond;
                 x_0= particle->firstStep()->preGlobalX()/Units::um;
                 y_0= particle->firstStep()->preGlobalY()/Units::um;
                 z_0= particle->firstStep()->preGlobalZ()/Units::um;
                 first_daughter = 0;
               }

               ///////////////////debugging segment
         //   if (particle->getParent()->isPrimary()){
         //      printf("- %s %s %f keV atomic number %d time %f parent %s %f keV  %d daugters death %f \n ",
         //             particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0, particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV, particle->nDaughters(), particle->lastStep()->postTime()/Units::picosecond-time_0);
         // }else if (particle->getParent()->getParent()->isPrimary()){
         //   printf("-- %s %s %f keV atomic number %d time %f parent %s %f keV %d daughters death %f \n",
         //          particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0, particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV,particle->nDaughters(), particle->lastStep()->postTime()/Units::picosecond-time_0);
         // }else if (particle->getParent()->getParent()->getParent()->isPrimary()){
         //   printf("--- %s %s %f keV atomic number %d  time %f parent %s %f keV  %d daughters death %f \n",
         //          particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0,  particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV, particle->nDaughters(), particle->lastStep()->postTime()/Units::picosecond-time_0);}else if (particle->getParent()->getParent()->getParent()->getParent()->isPrimary()){
         //     printf("---- %s %s %f keV atomic number %d  time %f parent %s %f keV %d daughters death %f \n",
         //            particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0,  particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV, particle->nDaughters(),particle->lastStep()->postTime()/Units::picosecond-time_0);}
           ///////////////////end of debugging segment
             }

           // segment to find the the time_0
           // auto seg0 = particle->segmentBegin();
           // auto segL = particle->segmentEnd();
           // auto seg = seg0;
           // for (;seg!=segL;++seg){
           //   if (seg->volumeName()=="Diode"){
           //     time_0=seg->stepBegin()->preTime()/Units::picosecond;
           //   }}
             if(!particle->isPrimary() && particle->atomicNumber()>1){

         //////////////////////////The tree hierarchy print///////////////////////////
         auto time_particle=particle->firstStep()->preTime()/Units::picosecond;
         auto time_frame=time_particle-time_0;
         std::fstream myfile;
         if (time_frame<0.1){
           myfile.open ("0_1.txt",std::fstream::app);
         }else if (time_frame<0.2){
           myfile.open ("0_2.txt", std::fstream::app);
         }else if (time_frame<0.3){
           myfile.open ("0_3.txt", std::fstream::app);
         }else if (time_frame<0.4){
           myfile.open ("0_4.txt", std::fstream::app);
         }else if (time_frame<0.5){
           myfile.open ("0_5.txt", std::fstream::app);
         }else if (time_frame<0.6){
           myfile.open ("0_6.txt", std::fstream::app);
         }else if (time_frame<0.7){
           myfile.open ("0_7.txt", std::fstream::app);
         }else if (time_frame<0.8){
           myfile.open ("0_8.txt", std::fstream::app);
         }else if (time_frame<0.9){
           myfile.open ("0_9.txt", std::fstream::app);
         }else if (time_frame<1.0){
           myfile.open ("1_0.txt", std::fstream::app);
         }else if (time_frame<1.1){
           myfile.open ("1_1.txt", std::fstream::app);
         }else if (time_frame<1.2){
           myfile.open ("1_2.txt", std::fstream::app);
         }else if (time_frame<1.3){
           myfile.open ("1_3.txt", std::fstream::app);
         }else if (time_frame<1.4){
           myfile.open ("1_4.txt", std::fstream::app);
         }else if (time_frame<1.5){
           myfile.open ("1_5.txt", std::fstream::app);}
         // }
       descendant_count+=1;
         auto x= (particle->firstStep()->preGlobalX()/Units::um)-x_0;
         auto y= (particle->firstStep()->preGlobalY()/Units::um)-y_0;
         auto z= (particle->firstStep()->preGlobalZ()/Units::um)-z_0;
         auto element=particle->atomicNumber();
         auto energy = particle->firstStep()->preEKin()/Units::keV;
         // auto time_death = particle->lastStep()->postTime()/Units::um;
           myfile <<element<<" "<<x<<" "<<y<<" "<<z<<" "<<energy<<std::endl;
           auto seg = particle->segmentBegin();
           auto segL = particle->segmentEnd();
           for (;seg!=segL;++seg){
           if (seg->volumeName()=="Diode"){
           auto step = seg->stepBegin();
           auto stepE= seg->stepEnd();
           for (;step!=stepE;++step){
             auto current_time=step->preTime()/Units::picosecond-time_0;
             // std::cout<<current_time/Units::picosecond<<std::endl;
         std::fstream allfile;
         if (current_time<0.1){
           allfile.open ("0_1_all.txt",std::fstream::app);
         }else if (current_time<0.2){
           allfile.open ("0_2_all.txt", std::fstream::app);
         }else if (current_time<0.3){
           allfile.open ("0_3_all.txt", std::fstream::app);
         }else if (current_time<0.4){
           allfile.open ("0_4_all.txt", std::fstream::app);
         }else if (current_time<0.5){
           allfile.open ("0_5_all.txt", std::fstream::app);
         }else if (current_time<0.6){
           allfile.open ("0_6_all.txt", std::fstream::app);
         }else if (current_time<0.7){
           allfile.open ("0_7_all.txt", std::fstream::app);
         }else if (current_time<0.8){
           allfile.open ("0_8_all.txt", std::fstream::app);
         }else if (current_time<0.9){
           allfile.open ("0_9_all.txt", std::fstream::app);
         }else if (current_time<1.0){
           allfile.open ("1_0_all.txt", std::fstream::app);
         }else if (current_time<1.1){
           allfile.open ("1_1_all.txt", std::fstream::app);
         }else if (current_time<1.2){
           allfile.open ("1_2_all.txt", std::fstream::app);
         }else if (current_time<1.3){
           allfile.open ("1_3_all.txt", std::fstream::app);
         }else if (current_time<1.4){
           allfile.open ("1_4_all.txt", std::fstream::app);
         }else if (current_time<1.5){
           allfile.open ("1_5_all.txt", std::fstream::app);}
         x= (step->preGlobalX()/Units::um)-x_0;
         y= (step->preGlobalY()/Units::um)-y_0;
         z= (step->preGlobalZ()/Units::um)-z_0;
         energy =step->preEKin()/Units::keV;
         allfile <<element<<" "<<x<<" "<<y<<" "<<z<<" "<<energy<<std::endl;
         allfile.close();
           }}}
       myfile.close();
         }

       // printf("Particle %s %s %f keV atomic number %d number of daughters %d \n",
              // particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->nDaughters());
     }
       }
     
   }
   }
 }
  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  hc.saveToFile("silicon_gif",true);

  return 0;
}


