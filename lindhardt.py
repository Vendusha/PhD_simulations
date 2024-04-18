from scipy.interpolate import interp1d
import numpy as np
hbar    = 1.0546e-27
speed_c = 3.e10
mass_e  = 9.11e-28
esuChg  = 4.8032e-10
alpha   = 7.29735e-3
mass_p  = 1.67e-24
amu     = 1.6605E-24
mass_si = 28.086*amu
mev2erg = 1.6022e-6
erg2mev = 624151.0
cm2barn = 1.0e24
# silicon properties needed for 
Na    = 6.022e23
mSi   = 28.086
rhoSi = 2.33
Ntarg = (Na*rhoSi)/mSi  # target density (number per cm3)
Nthick = 0.01
rhoSi = 2.33
Nthick = 0.01           # silicon target thickness (cm)
def Sig_Rec_WM(K_proj,trec_min):
  z_proj = 1.0
  z_targ = 14.0
  k_proj = K_proj * mev2erg             #----------- KE of projectile -------------#
  m_proj = mass_p
  e_proj = k_proj + (m_proj*speed_c**2)
  m_targ = mass_si
  e_targ = mass_si*speed_c**2
  p_proj = np.sqrt((e_proj/speed_c)**2-(m_proj*speed_c)**2)
  m_12_sq = (m_proj**2) + (m_targ**2) + 2.0*(e_proj/speed_c**2)*(e_targ/speed_c**2)
  m_12   = np.sqrt(m_12_sq)
  mu_rel = m_proj*m_targ/m_12
  beta_r = (p_proj*speed_c)/e_proj
  p_proj_pr = (p_proj*m_targ)/m_12
  Tmax = (2.0*m_targ*p_proj**2)/m_12_sq
  Tmin = trec_min*mev2erg              #-------- KE needed to displace atom ------#
  a0  = hbar**2 / (mass_e * esuChg**2)
  cTF = 0.88534
  aTF = (cTF * a0) / pow(z_targ,(1./3.))
  As = (hbar/(2.0*p_proj_pr*aTF))**2 * (1.13 + (3.76 * ((alpha*z_targ*z_proj)/(beta_r))**2 ))
  #
  a1 = (z_proj*z_targ*esuChg**2)**2
  a2 = (p_proj_pr*beta_r*speed_c)**2
  a3 = Tmax*As
  # print(Tmax)
  # Exclude values greater than kinematically allowed (i.e. > Tmax)
  if Tmin < Tmax:
    numer = np.pi*(a1/a2)*Tmax*(Tmax-Tmin)
    denom = (a3+Tmin)*(a3+Tmax)
    myvalue = numer/denom
  else:
    myvalue = 0.0

  return (myvalue)


def centers_from_borders(borders):
    """From histogram, calculates the centers of everything"""
    return borders[:-1] + np.diff(borders) / 2


def niel_alpha(KEalpha):
  KE_alpha = KEalpha   # alpha projectile energy (MeV)
  xvals=[]
  yvals=[]
  filo = open("NIELdata/NIEL_Alpha_Burke.dat","r")
  filo.readline()
  filo.readline()
  filo.readline()
  for line in filo:
    line = line.split()
    xval = float(line[0])
    xval = 4.0 * xval  # alpha NIEL data given per nucleon, I want alpha KE so x 4
    yval = float(line[1])
    yval = yval * rhoSi * Nthick
    # convert also to 100um silicon!!!!!!!!!!!!!!!!!!!!!!
    xvals.append(xval)
    yvals.append(yval)
  ynew = interp1d(xvals,yvals,kind='cubic')
  return(ynew(KE_alpha))

# https://www.bnl.gov/isd/documents/93416.pdf from here
def niel_si_lindhard_classic(e_Si):
  e_Si=e_Si*1.e6 # convert to eV
  k = 0.1462
  eps = 2.147e-5*e_Si
  g = 3.4008*eps**(1/6)+0.40244*eps**(3/4)+eps
  E_NIEL = e_Si*(1/(1+k*g))
  return(E_NIEL*1e-6)
   # absolute value in MeV
def NRT_dpa(e_Si): # input energy in MeV
  E_NIEL=niel_si_lindhard_classic(e_Si)*1e6 # convert to eV
  E_d = 21
  if E_NIEL>E_d and E_NIEL<2*E_d/0.8:
    return(1)
  elif E_NIEL>2*E_d/0.8:
    return(0.8*E_NIEL/(2*E_d))
  else:
    return(0)
def niel_si_lindhard_updated(e_Si): 
  k = 0.1462
  eps = 2.147e-5*e_Si
  g = 0.90656*eps**(1/6)+1.6812*eps**(3/4)+0.74422*eps
  E_NIEL = e_Si/(1+k*g)
  return E_NIEL # absolute value in MeV


# def niel_si(trec):
#   E_rec = trec*1.e6 # silicon target recoil energy (in eV)
#   Epsd  = 2.147e-5 * E_rec
#   gEpsd = Epsd + 0.40244*pow(Epsd,(3./4.)) + 3.4008*pow(Epsd,(1./6.))
#   Kd = 0.1462
#   denom = 1.0 + (Kd*gEpsd)
#   return(E_rec*1.e-6/denom) # absolute value in MeV
#   #return(1.0/denom) # NIEL fraction of recoil energy
