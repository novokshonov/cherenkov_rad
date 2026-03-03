import numpy as np
from scipy.optimize import root_scalar




def chr_angle(n, _gamma=100):
  # Cherenkov angle in crystal - see Fig1 in schemes.pdf
  # n - refractive index of media
  # _gamma - Lorentz factor
  beta = np.sqrt(1 - 1/_gamma**2) # reduced velocity
  theta_ch = np.arccos(1 / (n * beta))
  return theta_ch



def chr_angle_vac(n1, _theta_cr=np.pi/4, _gamma=100, _n2=1):
  # Cherenkov radiation angle after crossing boundary crystal/vacuum - see Fig. 1 in schemes.pdf
  # n1 - refractive index of the media
  # n2 - refractive index outside the media
  # _theta_cr - crystal rotation angle (see Fig. 1)
  # _gamma - Lorentz factor
  beta = np.sqrt(1 - 1/_gamma**2) # reduced velocity
  theta_ch = np.arccos(1 / (n * beta)) # ChR angle (Fig. 1)
  alpha2 = np.arcsin(np.sin(theta_ch - _theta_cr) * n1) # ChR to crystal normal angle (Fig. 1)
  theta_ch_vac = alpha2 + _theta_cr
  return theta_ch_vac



def refr_ind(lmbd):
  # Sellmeier equation - dependence of refractive index on wavelength of fused silica crystal
  # lmbd - wavelength in MICROMETERS !!!!!! 
  eps = 1 + 0.6961663 * lmbd**2 / (lmbd**2 - 0.0684043**2) + 0.4079426 * lmbd**2 / (lmbd**2 - 0.1162414**2) + 0.8974794 * lmbd**2 / (lmbd**2 - 9.896161**2) 
  n = np.sqrt(eps) # refractive index
  return n



def wvl_peak_pos(_theta_cr=0.3 _theta_obs=0.97, _gamma=100, _n2=1):
  # The function finds position of the spectral peak [nm] at given theta_cr and _theta_obs
  # _theta_cr and _theta_obs are shown in Fig.1 (schemes.pdf)
  # _n2 are refractive index ouside the crystal
  # We soppose the crystal is out of Fused Silica -> Sellmeyer Equation is taken below in min_func2
  beta = np.sqrt(1 - 1/_gamma**2) # Reduced velocity in crystal

  # The minimisation funtion. What is n1(lambda) at given _theta_cr and _theta_obs
  min_f1 = lambda n1: n1 * np.sin(np.arccos(1/(n1*beta)) - theta_cr) - _n2 * np.sin(theta_obs - theta_cr)
  # Here the n1(lambda) is found
  res1 = root_scalar(min_f1, bracket=[1.01, 1.99]) # The limits are also according to fused silica crystal

  # Now knowing the n1(lambda) we are finding the lambda corresponding to it
  min_f2 = lambda lmbd: refr_ind(lmbd) - res1.root
  # The same minimalisation function is used
  res2 = root_scalar(min_f2, bracket=[0.15, 4]) # Limints are lambda in micrometers
  return res2.root
