import G4CustomPyGen
import Core.Units as Units
import math
import Core.FindData
import numpy as np
import random
import G4GriffGen.GriffGen as Gen
def readAmBe():
    datafile = Core.FindData("SFGenerators","Neutron_CDF_Temp.txt")
    datafile_1 = Core.FindData("SFGenerators","Neutron_Energy_Temp.txt")
    with open(datafile) as f:
        neutronCDF=np.array(f.read().split(),dtype=np.float64)
    with open(datafile_1) as f_1:
        neutronEnergy=np.array(f_1.read().split(),dtype=np.float64)
    return (neutronCDF,neutronEnergy)

class NeutronGen(G4CustomPyGen.GenBase): 
#this is neutron generator
    def declare_parameters(self):
        self.addParameterDouble("SourceDistance",30) #distance of the source in cm
        self.addParameterDouble("BeamWidth",2) #width of the beam in cm
        self.addParameterDouble("BeamHeight",2) #height of the beam in cm
        self.addParameterDouble("NeutronEnergy",1000) #NeutronEnergy in meV
        #dimensions of the source
    def init_generator(self,gun):
        gun.set_type('neutron')
    def generate_event(self,gun):
        gun.set_direction(0,0,1)
        gun.set_type('neutron')
        gun.set_energy(max(0,self.randGauss(0.001*self.NeutronEnergy*Units.meV,self.NeutronEnergy*Units.meV)))

        xdir=self.rand(-self.BeamWidth/2,self.BeamWidth/2)*Units.cm
        ydir=self.rand(-self.BeamHeight/2,self.BeamHeight/2)*Units.cm
        zdir=self.SourceDistance*Units.cm
        gun.set_position(xdir,ydir,-zdir)
class ThermalizedGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterString("energy_histogram",
               "BeamMonitor/Thermalized.shist:ThermalizedNeutronZoom:eV")
        self.addParameterDouble("cone_opening_deg",45.0, 5.0, 180.0 )
    def init_generator(self,gun):
        gun.set_type('neutron')
        self._uzcut = math.cos(self.getParameterDouble("cone_opening_deg")*math.pi/180)
        gun.set_position(112*Units.mm,0,0)
        self._esampler = self.create_hist_sampler(self.energy_histogram)
    def generate_event(self,gun):
        gun.set_energy(self._esampler())
        while True:
            ux,uy,uz = self.randIsotropic()
            if uz>self._uzcut:
                break
        gun.set_direction(uz,ux,uy)
class ThermalizedGammaGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterString("energy_histogram",
               "BeamMonitor/Thermalized.shist:GammasFromThermalizer:MeV")
        self.addParameterDouble("cone_opening_deg",45.0, 5.0, 180.0 )
    def init_generator(self,gun):
        gun.set_type('gamma')
        self._uzcut = math.cos(self.getParameterDouble("cone_opening_deg")*math.pi/180)
        gun.set_position(112*Units.mm,0,0)
        self._esampler = self.create_hist_sampler(self.energy_histogram)
    def generate_event(self,gun):
        gun.set_energy(self._esampler())
        while True:
            ux,uy,uz = self.randIsotropic()
            if uz>self._uzcut:
                break
        gun.set_direction(uz,ux,uy)

class AmBeGen(G4CustomPyGen.GenBase): 
#this is an AmBeGen taken from an external code
    def declare_parameters(self):
        #dimensions of the source
        self.addParameterDouble("Radius_Source_cm", 1.5) 
        self.addParameterDouble("Length_cm", 4)
    def init_generator(self,gun):
        self.neutronCDF,self.neutronEnergy=readAmBe()
        gun.set_position(0,0,0)
        gun.set_type('neutron')

    def generate_event(self,gun):
        gun.set_type('neutron')
        gun.set_random_direction()
        xdir=self.rand(-1,1)*self.Radius_Source_cm*Units.cm
        ydir=self.rand(-1,1)*self.Radius_Source_cm*Units.cm
        zdir=self.rand(-0.5,0.5)*self.Length_cm*Units.cm
        gun.set_position(xdir,ydir,zdir)
        prob=self.rand()
        idx=np.searchsorted(self.neutronCDF,prob,side="left")
        if idx==757:
            idx_high=idx
            idx=idx-1
        else:
            idx_high=idx+1
        split = (prob - self.neutronCDF[idx])/(self.neutronCDF[idx_high] - self.neutronCDF[idx]);
        energy = (self.neutronEnergy[idx]+split*(self.neutronEnergy[idx_high] - self.neutronEnergy[idx]))*Units.MeV;
        gun.set_energy(max(0,energy))
        gun.fire()
        if (energy >0.5*Units.MeV and energy <6*Units.MeV): # if the neutron is within the range of 0.5-6, it emmits, 4.4 MeV Gamma. If it is 0.5-1.9, it emmits 3.21483 and additional 4.4
            gun.set_random_direction()
            gun.set_type('gamma')
            gun.set_energy(max(0,self.randGauss(0.022*Units.MeV,4.43803*Units.MeV))) #first parameter is the Gaussian smearing, 0.44 is chosen more or less randomly (0.5% of the peak)
            #gun.set_energy(4.43803*Units.MeV)
            gun.fire()
            if energy<1.9*Units.MeV:
                gun.set_random_direction()
                #gun.set_energy(max(0,self.randGauss(0.015*Units.MeV,3.21483*Units.MeV))) #first parameter is the Gaussian smearing, 0.03 is chosen as 0.5 percent of the peak
                gun.set_energy(3.21483*Units.MeV)
                gun.fire()

class IsoAngle(G4CustomPyGen.GenBase): 
#this is a single energy neutron generator with an isotropic cone 
    def declare_parameters(self):
        self.addParameterDouble("cone_opening_deg",30.0, 5.0, 180.0 ) 
    def init_generator(self,gun):
        self._uzcut = math.cos(self.getParameterDouble("cone_opening_deg")*math.pi/180)
        gun.set_position(0,0,0)
        gun.set_energy(25*Units.meV)
        gun.set_type('neutron')
    def generate_event(self,gun):
        while True:
            ux,uy,uz = self.randIsotropic()
            if uz>self._uzcut: 
                break
        gun.set_direction(-uz,-ux,-uy)
        gun.set_position(0,0,0)
        gun.set_energy(25*Units.meV)

class SinglePath(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("energy_MeV",1.6) # 120 GeV
        self.addParameterString("particle_type",'pi+')
        # self.addParameterDouble("GaussianCenter" ,22.02) # center of the Gaussian in ms
        # self.addParameterDouble("GaussianSpread" ,0.11) # % of the Gaussian FWHM in absolute numbers
    def init_generator(self,gun):
        gun.set_type(self.particle_type)
        gun.set_energy((self.energy_MeV)*Units.MeV)
        gun.set_position (0,0,-20*Units.cm)
    def generate_event(self,gun):
        gun.set_direction(0,0,1)

#####################################################

class energyGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterDouble("energy_MeV", 0.2)
        self.addParameterString("particle_type",'proton')
    def init_generator(self,gun):
        gun.set_type(self.particle_type)
        gun.set_energy((self.energy_MeV)*Units.MeV)
        # gun.set_position()
    def generate_event(self,gun):
        gun.set_direction(0,0,1)
        # x=self.rand(-0.05,0.05)*Units.cm
        # y=self.rand(-0.05,0.05)*Units.cm
        gun.set_position(0,0,-0.1*Units.cm)
        # gun.set_position (0,0,-20*Units.cm)
        # gun.set_position(0.1*Units.cm,0.1*Units.cm,-0.1*Units.cm)

     	# gun.set_random_direction()

#####################################################
        
def launch(geo):
    import G4Launcher
    launcher = G4Launcher()
    #launcher.addParameterString('event_gen','SinglePath')
    launcher.addParameterString('event_gen','else') 
    launcher.addParameterString('griff_file','else.txt') 
    #geometry:
    launcher.setGeo(geo)

    #generator
    if launcher.getParameterString('event_gen')=='IsoAngle':
        gen = IsoAngle()
    elif launcher.getParameterString('event_gen')=='energyGenPKA':
        gen=energyGen()
    elif launcher.getParameterString('event_gen')=='GriffGen':
        gen=Gen.create()
        gen.input_file=launcher.getParameterString("griff_file")
        gen.primary_only = False
    elif launcher.getParameterString('event_gen')=='energyGen':
        gen=energyGen()
    elif launcher.getParameterString('event_gen')=='SinglePath':
        gen=SinglePath()
    elif launcher.getParameterString('event_gen')=='ThermalizedGammaGen':
        gen=ThermalizedGammaGen()
    elif launcher.getParameterString('event_gen')=='AmBeGen':
        gen=AmBeGen()
    elif launcher.getParameterString('event_gen')=='ThermalizedGen':
        gen=ThermalizedGen()
    elif launcher.getParameterString('event_gen')=="SiliconGen":
        import G4StdGenerators.FlexGenDefaultSpherical
        gen = G4StdGenerators.FlexGenDefaultSpherical.create()
    else:
        import G4StdGenerators.FlexGenDefaultSpherical
        gen = G4StdGenerators.FlexGenDefaultSpherical.create()
        # gen.particleName = 'neutron'
        # print(gen.__dict__)
        gen.pdgCode = 2112
        gen.fixed_energy_eV = 200e3
    launcher.setGen(gen)

    #filter:
    if launcher.getParameterString('event_gen')=='energyGenPKA':
        import G4CollectFilters.StepFilterPKA as F
        launcher.setFilter(F.create())


    #general stuff:
    launcher.cmd_preinit('/process/eLoss/StepFunction 0.1 0.001 um')
    launcher.cmd_preinit('/process/eLoss/minKinEnergy 1 eV')
    #launcher.cmd_postinit('/run/setCutForAGivenParticle proton 0 um')
    launcher.setOutput('silicon','REDUCED')


    #launch:
    launcher.go()
