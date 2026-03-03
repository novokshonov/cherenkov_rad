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

def peak_lmbd_on_cr_angle(theta_cr, _n2=1, _obs_angle=0.97, _gamma=100):
  beta = np.sqrt(1 - 1/_gamma**2)
  min_func1 = lambda n1: np.arcsin(_n2 * np.sin(_obs_angle - theta_cr) / n1) - np.arccos(1 / (n1 * beta)) - theta_cr
  sol1 = root_scalar(min_func1(n1), bracket=[0.0, 2], method='brentq')
  min_func2 = lambda lmbd: 1 + 0.6961663 * lmbd**2 / (lmbd**2 - 0.0684043**2) + 0.4079426 * lmbd**2 / (lmbd**2 - 0.1162414**2) + 0.8974794 * lmbd**2 / (lmbd**2 - 9.896161**2) - sol**2
  sol2 = root_scalar(min_func2(lmbd), bracket=[0.1, 2], method='brentq')
  return sol2
