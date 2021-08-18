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
  
  std::ofstream ofs;
  ofs.open("Events.txt", std::ofstream::out | std::ofstream::trunc);
  ofs.close();
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
       first_daughter = 1;

           if (track->nDaughters()>0){
             if (new_event !=0){
               ++event_no;
               std::fstream myfile;
               myfile.open ("Events.txt",std::fstream::app);
               myfile <<"0"<<" "<<"0"<<" "<<"0"<<" "<<"0"<<" "<<event_no<<std::endl;
               myfile.close();
             new_event=0;}

             descendantFilter->setAncestor(track);


           while (auto particle = descendant.next()){
             if (!particle->isPrimary()){
               if (particle->getParent()->isPrimary()&& particle->nDaughters()>0 and first_daughter == 1){
                 new_event=1;
                 x_0= particle->firstStep()->preGlobalX()/Units::um;
                 y_0= particle->firstStep()->preGlobalY()/Units::um;
                 z_0= particle->firstStep()->preGlobalZ()/Units::um;
                 first_daughter = 0;
               }

             if(!particle->isPrimary() && particle->atomicNumber()>1){

         //////////////////////////The tree hierarchy print///////////////////////////
         std::fstream myfile;
         myfile.open ("Events.txt", std::fstream::app);
       descendant_count+=1;
         auto x= (particle->firstStep()->preGlobalX()/Units::um)-x_0;
         auto y= (particle->firstStep()->preGlobalY()/Units::um)-y_0;
         auto z= (particle->firstStep()->preGlobalZ()/Units::um)-z_0;
         auto element=particle->atomicNumber();
         auto PKA_atomic=particle->atomicNumber();
         if (particle->getParent()->isPrimary()){
           PKA_atomic=particle->atomicNumber();
             }else if (particle->getParent()->getParent()->isPrimary()){
           PKA_atomic=particle->getParent()->atomicNumber();
             }else if (particle->getParent()->getParent()->getParent()->isPrimary()){
           PKA_atomic=particle->getParent()->getParent()->atomicNumber();
             }else if (particle->getParent()->getParent()->getParent()->getParent()->isPrimary()){
           PKA_atomic=particle->getParent()->getParent()->getParent()->atomicNumber();
             }else{
           PKA_atomic=14;}
         auto energy = particle->firstStep()->preEKin()/Units::keV;
         // auto time_death = particle->lastStep()->postTime()/Units::um;
         myfile <<element<<" "<<x<<" "<<y<<" "<<z<<" "<<energy<<" "<<PKA_atomic<<std::endl;
       myfile.close();
         }

       // printf("Particle %s %s %f keV atomic number %d number of daughters %d \n",
              // particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->nDaughters());
     }
       }
           }
     
   }
   }
 }
  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  hc.saveToFile("silicon_gif",true);

  return 0;
}


