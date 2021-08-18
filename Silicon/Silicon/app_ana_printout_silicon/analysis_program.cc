#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griffon for more info.

int main(int argc,char**argv) {
  
  GriffDataReader dr(argc,argv);
  GriffAnaUtils::SegmentIterator diode(&dr);
  diode.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Diode"));


  dr.setup()->dump();

  if (dr.setup()->geo().getName()!="G4GeoSilicon/GeoSilicon")
    return 1;
  ///////////////////Printing info//////////////////////
  // auto energy_init = 200;
  // auto energy_init= dr.setup()->gen().getParameterDouble("energy_MeV");
  // const std::string particle_type= dr.setup()->gen().getParameterString("particle_type");
  // std::cout<< "The generator shooted " << particle_type <<" of the energy "<< energy_init << "MeV ."<<std::endl;
  auto descendantFilter = new GriffAnaUtils::TrackFilter_Descendant();
   // descendantFilter->setNegated();
  GriffAnaUtils::TrackIterator descendant(&dr);
  descendant.addFilter(descendantFilter);
  //////////////////////Filter primary particles////////////////////////////////////
  GriffAnaUtils::TrackIterator primary_particles(&dr);
  primary_particles.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  int event_no=1;
  //Loop over events and extract info via Griff interface:
  while (dr.loopEvents()) {
    for (auto track = dr.trackBegin();track!=dr.trackEnd();++track) {
      // printf("Found event %i/%i [ntrks=%i] ===========================================\n",
                                                                                         // dr.runNumber(),dr.eventNumber(),dr.nTracks());


     if (track->isPrimary()){
       if (track->nDaughters()>0){
         auto time_0=track->firstStep()->preTime()/Units::picosecond;
         auto x_0= track->firstStep()->preGlobalX()/Units::um;
         auto y_0= track->firstStep()->preGlobalY()/Units::um;
         auto z_0= track->firstStep()->preGlobalZ()/Units::um;
         ++event_no;
             descendantFilter->setAncestor(track);
             printf("Found event %i/%i [ntrks=%i] ===========================================\n",
                                                                        dr.runNumber(),dr.eventNumber(),dr.nTracks());
             printf("Primary particle. %s %s %f keV atomic number %d number of daughters %d \n",
                    track->pdgNameCStr(), track->creatorProcessCStr(), track->startEKin()/Units::keV, track->atomicNumber(), track->nDaughters());

           while (auto particle = descendant.next()){
             if (!particle->isPrimary()){
               auto x_print= particle->firstStep()->preGlobalX()/Units::um;
               auto y_print= particle->firstStep()->preGlobalY()/Units::um;
               auto z_print= particle->firstStep()->preGlobalZ()/Units::um;
               ///////////////////debugging segment
           if (particle->getParent()->isPrimary()){
              printf("- %s %s %f keV atomic number %d time %f parent %s %f keV  %d daugters death %f x: %f, y:%f, z:%f\n ",
                     particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0, particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV, particle->nDaughters(), particle->lastStep()->postTime()/Units::picosecond-time_0, x_0,y_0,z_0);
         }else if (particle->getParent()->getParent()->isPrimary()){
             printf("-- %s %s %f keV atomic number %d time %f parent %s %f keV %d daughters death %f  x: %f, y: %f, z: %f\n",
                    particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0, particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV,particle->nDaughters(), particle->lastStep()->postTime()/Units::picosecond-time_0,x_print-x_0,y_print-y_0,z_print-z_0);
         }else if (particle->getParent()->getParent()->getParent()->isPrimary()){
           printf("--- %s %s %f keV atomic number %d  time %f parent %s %f keV  %d daughters death %f \n",
                  particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0,  particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV, particle->nDaughters(), particle->lastStep()->postTime()/Units::picosecond-time_0);}else if (particle->getParent()->getParent()->getParent()->getParent()->isPrimary()){
             printf("---- %s %s %f keV atomic number %d  time %f parent %s %f keV %d daughters death %f \n",
                    particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->firstStep()->preTime()/Units::picosecond-time_0,  particle->getParent()->pdgNameCStr(), particle->getParent()->startEKin()/Units::keV, particle->nDaughters(),particle->lastStep()->postTime()/Units::picosecond-time_0);}
           ///////////////////end of debugging segment
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


