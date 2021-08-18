//After you have modified this file as needed for your project you can remove this line and commit <NOCOMMIT>

/////////////////////////////////////////
// Declaration of our geometry module: //
/////////////////////////////////////////

#include "G4Interfaces/GeoConstructPyExport.hh"

class GeoSilicon : public G4Interfaces::GeoConstructBase
{
public:
  GeoSilicon();
  virtual ~GeoSilicon(){}
  virtual G4VPhysicalVolume* Construct();
  //(add more member functions and data here if needed)
};

PYTHON_MODULE { GeoConstructPyExport::exportGeo<GeoSilicon>("GeoSilicon"); }

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

#include "G4Box.hh"
#include "G4Orb.hh"

//As an example, we put a spherical sample in front of a box representing a
//detector tube (using the Griff analysis to look for positions where neutrons
//enter the detector).

GeoSilicon::GeoSilicon()
  : GeoConstructBase("G4GeoSilicon/GeoSilicon")
{
  //Free parameters, tunable from python and command-line, printed in log-files,
  //stored in Griff-files. The values specified below will be the default
  //values. Avoid the temptation to instead hard-code "interesting" parameters
  //in the Construct() method:

  addParameterDouble("diode_posz_mm",0.0);
  addParameterDouble("diode_thickness_mm",0.1);
  addParameterDouble("diode_size_mm",1.0);
  addParameterString("material_diode","G4_Si");
  // addParameterString("material_lab",
                     // "IdealGas:formula=0.7*Ar+0.3*CO2:temp_kelvin=293.15:pressure_atm=2.0");

  //Note: It is possible and easy (but done done here to keep the example
  //      simple) to impose constraints on parameter ranges, etc., making sure
  //      the geometry will only be build with sensible parameters. For
  //      instance, one should not put the sample radius larger than the
  //      detector-sample distance.
}

G4VPhysicalVolume* GeoSilicon::Construct()
{
  //Parameters (converting to G4 units immediately as is best practice):
  const double diode_posz = getParameterDouble("diode_posz_mm")*Units::mm;
  const double diode_size = getParameterDouble("diode_size_mm")*Units::mm;
  const double diode_thickness_um = getParameterDouble("diode_thickness_mm")*Units::mm;
  // const double det_size = getParameterDouble("detector_size_cm")*Units::cm;
  // const double det_depth = 1.0*Units::cm;//ok to hardcode non-interesting parameters
  // const double det_sample_dist = getParameterDouble("detector_sample_dist_cm")*Units::cm;
  auto mat_diode = getParameterMaterial("material_diode");
  // auto mat_lab = getParameterMaterial("material_lab");
  auto mat_vacuum = getMaterial("Vacuum");

  //Notice that materials are created above via the NamedMaterialProvider. Avoid
  //the temptation to instead create your own G4Material instances in the code,
  //since it invites bugs and reduces flexibility, reusability and readability!
  //
  //Find more info at: https://confluence.esss.lu.se/display/DG/NamedMaterials

  //World volume (must be big enough for the sample and detector to fit inside):
  // const double dz_world = 1.001 * (std::abs<double>(_posz)+sample_radius+det_depth+det_sample_dist);
  // const double dxy_world = 1.001 * std::max<double>(0.5*det_size,sample_radius);
  auto worldvols = place(new G4Box("World", 0.001*Units::m, 0.001*Units::m, 0.001*Units::m),mat_vacuum,0);
  
  auto lvWorld = worldvols.logvol;

  //Sample:
  place(new G4Box("Diode",0.5*diode_size, 0.5*diode_size, 0.5*diode_thickness_um), mat_diode,0,0,diode_posz,lvWorld,YELLOW);

  //Detector:
  // place(new G4Box("Detector",0.5*det_size,0.5*det_size,0.5*det_depth),mat_det,
        // 0,0,sample_posz+det_sample_dist+0.5*det_depth,lvWorld,RED);

  //Note that we added visualisation colours for volumes above. Other colours
  //include BLUE, GREEN, ORANGE, CYAN, PURPLE, SILVER, GOLD, WHITE, BROWN, ...

  return worldvols.physvol;
}
