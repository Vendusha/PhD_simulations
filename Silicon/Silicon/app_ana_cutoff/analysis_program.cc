#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <iomanip>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.
std::string remove_extension(const std::string& filename) {
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot); 
}


int main(int argc,char**argv) {
  GriffDataReader dr(argc-1,argv);
  GriffAnaUtils::SegmentIterator diode(&dr);
  diode.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Diode"));
  float energy_cutoff=std::stof(argv[2]);
  // std::cout << "Choose energy cutoff in keV: "<<std::endl;
  // std::cin >> energy_cutoff;
  // std::cout << energy_cutoff<<std::endl;
  std::stringstream stream;
  stream << std::fixed << std::setprecision(3) << energy_cutoff;
  std::string s = stream.str();
  auto filename = remove_extension(dr.setup()->gen().getParameterString("input_file"))+"_"+s;
  std::ofstream ofs;
  ofs.open(filename+".txt", std::ofstream::out | std::ofstream::trunc);
  ofs.close();

   dr.setup()->dump();

  if (dr.setup()->geo().getName()!="G4GeoSilicon/GeoSilicon")
    return 1;
  ///////////////////Printing info//////////////////////
  auto descendantFilter = new GriffAnaUtils::TrackFilter_Descendant();
  GriffAnaUtils::TrackIterator descendant(&dr);
  descendant.addFilter(descendantFilter);
  auto descendantFilterSecondGen = new GriffAnaUtils::TrackFilter_Descendant();
  GriffAnaUtils::TrackIterator descendant2(&dr);
  descendant2.addFilter(descendantFilterSecondGen);
  //////////////////////Filter primary particles////////////////////////////////////
  GriffAnaUtils::TrackIterator primary_particles(&dr);
  primary_particles.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  // SimpleHists::HistCollection hc;
  //Loop over events and extract info via Griff interface:
  // auto energy_cutoff = 0.000; //try to write it in keV
  std::fstream myfile;
  while (dr.loopEvents()) {
    myfile.open(filename+".txt", std::fstream::app);
 for (auto track = dr.trackBegin(); track!=dr.trackEnd(); ++track) {
         auto energy = track -> firstStep()->preEKin()/Units::keV;
         if (track -> isPrimary() && energy > energy_cutoff ) {
         auto x = (track -> firstStep() -> preGlobalX()/Units::um);
         auto y = (track -> firstStep() -> preGlobalY()/Units::um);
         auto z = (track -> firstStep() -> preGlobalZ()/Units::um);
         auto x_end = (track -> lastStep() -> postGlobalX()/Units::um);
         auto y_end = (track -> lastStep() -> postGlobalY()/Units::um);
         auto z_end = (track -> lastStep() -> postGlobalZ()/Units::um);
         auto element = track -> atomicNumber();
         printf("%s %s %f keV, daughters %d x:%f y:%f z:%f \n",
                track -> pdgNameCStr(), track -> creatorProcessCStr(), track -> startEKin()/Units::keV,track -> nDaughters(),x,y,z);
         myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
       if (track -> nDaughters() > 0){
       for(uint i=0; i<track->nDaughters(); ++i) {
         auto particle = track -> getDaughter(i);
         energy = particle -> firstStep()->preEKin()/Units::keV;
         if (energy > energy_cutoff ) {
            x = (particle -> firstStep() -> preGlobalX()/Units::um);
            y = (particle -> firstStep() -> preGlobalY()/Units::um);
            z = (particle -> firstStep() -> preGlobalZ()/Units::um);
            x_end = (particle -> lastStep() -> postGlobalX()/Units::um);
            y_end = (particle -> lastStep() -> postGlobalY()/Units::um);
            z_end = (particle -> lastStep() -> postGlobalZ()/Units::um);
            element = particle -> atomicNumber();
            printf("-%s %s %f keV, daughters %d x:%f y:%f z:%f \n",
                   particle -> pdgNameCStr(), particle -> creatorProcessCStr(), particle -> startEKin()/Units::keV,particle -> nDaughters(),x,y,z);
           myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
         if (particle->nDaughters()>0){
           for(uint ii=0; ii<particle->nDaughters(); ++ii) {
             auto particle2 = particle -> getDaughter(ii);
             energy = particle2 -> firstStep()->preEKin()/Units::keV;
             if (energy > energy_cutoff ) {
                x = (particle2 -> firstStep() -> preGlobalX()/Units::um);
                y = (particle2 -> firstStep() -> preGlobalY()/Units::um);
                z = (particle2 -> firstStep() -> preGlobalZ()/Units::um);
                x_end = (particle2 -> lastStep() -> postGlobalX()/Units::um);
                y_end = (particle2 -> lastStep() -> postGlobalY()/Units::um);
                z_end = (particle2 -> lastStep() -> postGlobalZ()/Units::um);
                element = particle2 -> atomicNumber();
                printf("--%s %s %f keV, daughters %d x:%f y:%f z:%f \n",
                       particle2 -> pdgNameCStr(), particle2 -> creatorProcessCStr(), particle2 -> startEKin()/Units::keV,particle2 -> nDaughters(),x,y,z);
               myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
               if (particle2->nDaughters()>0){
               for(uint iii=0; iii<particle2->nDaughters(); ++iii) {
                 auto particle3 = particle2 -> getDaughter(iii);
                 energy = particle3 -> firstStep()->preEKin()/Units::keV;
                 if (energy > energy_cutoff ) {
                   x = (particle3 -> firstStep() -> preGlobalX()/Units::um);
                   y = (particle3 -> firstStep() -> preGlobalY()/Units::um);
                   z = (particle3 -> firstStep() -> preGlobalZ()/Units::um);
                   x_end = (particle3 -> lastStep() -> postGlobalX()/Units::um);
                   y_end = (particle3 -> lastStep() -> postGlobalY()/Units::um);
                   z_end = (particle3 -> lastStep() -> postGlobalZ()/Units::um);
                   element = particle3 -> atomicNumber();
                   printf("---%s %s %f keV, daughters %d x:%f y:%f z:%f \n",
                          particle3 -> pdgNameCStr(), particle3 -> creatorProcessCStr(), particle3 -> startEKin()/Units::keV,particle3 -> nDaughters(),x,y,z);
                   myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                   if (particle3->nDaughters()>0){
                   for(uint iiii=0; iiii<particle3->nDaughters(); ++iiii) {
                     auto particle4 = particle3 -> getDaughter(iiii);
                     energy = particle4 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       x = (particle4 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle4 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle4 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle4 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle4 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle4 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle4 -> atomicNumber();
                       printf("----%s %s %f keV, daughters %d x:%f y:%f z:%f \n",
                              particle4 -> pdgNameCStr(), particle4 -> creatorProcessCStr(), particle4 -> startEKin()/Units::keV,particle4 -> nDaughters(),x,y,z);
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                    if (particle4->nDaughters()>0){
                   for(uint iiiii=0; iiiii<particle4->nDaughters(); ++iiiii) {
                     auto particle5 = particle4 -> getDaughter(iiiii);
                     energy = particle5 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       printf("-----%s %s %f keV atomic number %d number of daughters %d \n",
                              particle5 -> pdgNameCStr(), particle5 -> creatorProcessCStr(), particle5 -> startEKin()/Units::keV, particle5 -> atomicNumber(), particle5 -> nDaughters());
                       x = (particle5 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle5 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle5 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle5 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle5 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle5 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle5 -> atomicNumber();
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                    if (particle5->nDaughters()>0){
                   for(uint iiiiii=0; iiiiii<particle5->nDaughters(); ++iiiiii) {
                     auto particle6 = particle5 -> getDaughter(iiiiii);
                     energy = particle6 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       printf("------%s %s %f keV atomic number %d number of daughters %d \n",
                              particle6 -> pdgNameCStr(), particle6 -> creatorProcessCStr(), particle6 -> startEKin()/Units::keV, particle6 -> atomicNumber(), particle6 -> nDaughters());
                       x = (particle6 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle6 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle6 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle6 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle6 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle6 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle6 -> atomicNumber();
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                       if (particle6->nDaughters()>0){
                   for(uint iiiiiii=0; iiiiiii<particle6->nDaughters(); ++iiiiiii) {
                     auto particle7 = particle6 -> getDaughter(iiiiiii);
                     energy = particle7 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       printf("-------%s %s %f keV atomic number %d number of daughters %d \n",
                              particle7 -> pdgNameCStr(), particle7 -> creatorProcessCStr(), particle7 -> startEKin()/Units::keV, particle7 -> atomicNumber(), particle7 -> nDaughters());
                       x = (particle7 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle7 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle7 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle7 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle7 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle7 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle7 -> atomicNumber();
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                       if (particle7->nDaughters()>0){
                   for(uint iiiiiiii=0; iiiiiiii<particle7->nDaughters(); ++iiiiiiii) {
                     auto particle8 = particle7 -> getDaughter(iiiiiiii);
                     energy = particle8 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       printf("--------%s %s %f keV atomic number %d number of daughters %d \n",
                              particle8 -> pdgNameCStr(), particle8 -> creatorProcessCStr(), particle8 -> startEKin()/Units::keV, particle8 -> atomicNumber(), particle8 -> nDaughters());
                       x = (particle8 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle8 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle8 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle8 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle8 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle8 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle8 -> atomicNumber();
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                       if (particle8->nDaughters()>0){
                   for(uint iiiiiiiii=0; iiiiiiiii<particle8->nDaughters(); ++iiiiiiiii) {
                     auto particle9 = particle8 -> getDaughter(iiiiiiiii);
                     energy = particle9 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       printf("---------%s %s %f keV atomic number %d number of daughters %d \n",
                              particle9 -> pdgNameCStr(), particle9 -> creatorProcessCStr(), particle9 -> startEKin()/Units::keV, particle9 -> atomicNumber(), particle9 -> nDaughters());
                       x = (particle9 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle9 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle9 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle9 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle9 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle9 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle9 -> atomicNumber();
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;
                       if (particle9->nDaughters()>0){
                   for(uint iiiiiiiiii=0; iiiiiiiiii<particle9->nDaughters(); ++iiiiiiiiii) {
                     auto particle10 = particle9 -> getDaughter(iiiiiiiiii);
                     energy = particle10 -> firstStep()->preEKin()/Units::keV;
                     if (energy > energy_cutoff ) {
                       printf("----------%s %s %f keV atomic number %d number of daughters %d \n",
                              particle10 -> pdgNameCStr(), particle10 -> creatorProcessCStr(), particle10 -> startEKin()/Units::keV, particle10 -> atomicNumber(), particle10 -> nDaughters());
                       x = (particle10 -> firstStep() -> preGlobalX()/Units::um);
                       y = (particle10 -> firstStep() -> preGlobalY()/Units::um);
                       z = (particle10 -> firstStep() -> preGlobalZ()/Units::um);
                       x_end = (particle10 -> lastStep() -> postGlobalX()/Units::um);
                       y_end = (particle10 -> lastStep() -> postGlobalY()/Units::um);
                       z_end = (particle10 -> lastStep() -> postGlobalZ()/Units::um);
                       element = particle10 -> atomicNumber();
                       myfile << element <<" "<< x <<" "<< y <<" "<< z <<" "<< energy <<" "<< x_end <<" "<< y_end <<" "<< z_end <<" "<<std::endl;

                     }}}
                     }}}
                     }}}
                     }}}
                     }}}
                     }}}
                     }}}
                 }}}
             } // if energy_cutoff 3nd generation
           } // if energy_cutoff 2nd generation
         }//end of looping in daughters of the second generation
       } //end of looping in daughters of the first generation
       } //if number of daughters of first generation >0
         } //if track is primary
         }
   }
 }
  myfile.close();
  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  // hc.saveToFile("silicon_gif",true);

  return 0;
}


