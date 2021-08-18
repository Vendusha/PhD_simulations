#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Utils/ArrayMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <iostream>

//After you have modified this file as needed for your project you can remove this line and commit <NOCOMMIT>

//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

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
  auto energy_init= dr.setup()->gen().getParameterDouble("energy_MeV");
  // auto energy_init= dr.setup()->gen().getParameterDouble("fixed_energy_eV");
  const std::string particle_type= dr.setup()->gen().getParameterString("particle_type");
  std::cout<< "The generator shooted " << particle_type <<" of the energy "<< energy_init << "MeV ."<<std::endl;
  // GriffAnaUtils::TrackIterator neutrons(&dr); 
  // neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112)); //PDGCode for neutrons
  auto descendantFilter = new GriffAnaUtils::TrackFilter_Descendant();
   // descendantFilter->setNegated();
  GriffAnaUtils::TrackIterator descendant(&dr);
  descendant.addFilter(descendantFilter);
  // GriffAnaUtils::TrackIterator ti(&dr);
  // ti.addFilter(new DaughterFilter);

  //Book histograms (see https://confluence.esss.lu.se/display/DG/SimpleHists for more info):
  SimpleHists::HistCollection hc;
  auto nbins = 8192;
  auto nbins_dim=1024;
  // auto h_recoil_average = hc.book1D("Recoil_spectra_Cascade",nbins, 0, energy_init/1000, "Recoil_spectra_Cascade");
  auto h_recoil_Cascade = hc.book1D("Recoil_spectra_Cascade",nbins, 0, energy_init/1000, "Recoil_spectra_Cascade");
  auto h_Daughter_Cascade = hc.bookCounts("Cascade_particle_type", "Cascade_particle_type");
  auto h_particle_yz = hc.book2D("YZ_distribution_neutrons",nbins_dim,-10,10,nbins_dim,-10,10,"YZ_distribution_neutrons"); //two dimensional histogram in the direction of the detector, visualise the incoming neutrons //values are in centimeters
  h_particle_yz->setXLabel("y_position [cm]");
  h_particle_yz->setYLabel("z_position [cm]");

  auto Si28_Cascade = h_Daughter_Cascade ->addCounter("Si28");
  auto Si29_Cascade = h_Daughter_Cascade ->addCounter("Si29");
  auto Si30_Cascade = h_Daughter_Cascade ->addCounter("Si30");
  auto others_Cascade = h_Daughter_Cascade ->addCounter("others");
  h_recoil_Cascade -> setXLabel("Energy [keV]");
  h_recoil_Cascade -> setYLabel("Frequency [PKA-1 keV-1]");
  //Loop over events and extract info via Griff interface:
  while (dr.loopEvents()) {
    printf("Found event %i/%i [ntrks=%i] ===========================================\n",
                                      dr.runNumber(),dr.eventNumber(),dr.nTracks());
 for (auto track = dr.trackBegin();track!=dr.trackEnd();++track) {
   if (track->isPrimary()){
       descendantFilter->setAncestor(track);
     while (auto particle = descendant.next()){
       printf("Particle %s %s %f keV atomic number %d number of daughters %d \n",
              particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->nDaughters());
   if (!particle->isPrimary()){
     if (particle->getParent()->isPrimary()){
     printf("- %s %s %f keV atomic number %d  parent %s \n",
            particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber(), particle->getParent()->pdgNameCStr());
     }else if (particle->getParent()->getParent()->isPrimary()){
     printf("--- %s %s %f keV atomic number %d  \n",
     particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber());
     }else if (particle->getParent()->getParent()->getParent()->isPrimary()){
     printf("----- %s %s %f keV atomic number %d  \n",
     particle->pdgNameCStr(), particle->creatorProcessCStr(), particle->startEKin()/Units::keV, particle->atomicNumber());
     }
   }
   if (particle->pdgName()=="Si28"){
                Si28_Cascade+=1;
               }
              else if (particle->pdgName()=="Si29"){
                Si29_Cascade+=1;}
              else if (particle->pdgName()=="Si30"){
                Si30_Cascade+=1;}
              else{
                // printf("Neutron daughter info %s %s %f keV atomic number %d \n",
                       // particle->pdgNameCStr(),particle->creatorProcessCStr(),particle->startEKin()/Units::keV,particle->atomicNumber());
                // printf("Particle is %s /n",particle->pdgNameCStr());
                others_Cascade+=1;}
   if (particle->startEKin()*1000>0.021){
     h_recoil_Cascade ->fill(particle->startEKin()*1000); 
   }
   if (particle->atomicNumber()==14){
     h_particle_yz->fill(particle->firstStep()->preGlobalY()/Units::um, particle->firstStep()->preGlobalZ()/Units::um, particle->weight());}
   }
 }}

  }

  //Save histograms to a file which can be browsed with ess_simplehists_browse:
  hc.saveToFile("silicon",true);

  return 0;
}


