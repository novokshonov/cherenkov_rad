import numpy as np

def chr_angle(n, _gamma=100):
  # n - refractive index of media
  # _gamma - Lorentz factor
  beta = np.sqrt(1 - 1/_gamma**2) # reduced velocity
  theta_ch = np.arccos(1 / (n * beta))
  return theta_ch

def chr_angle_vac(n1, _theta_cr=np.pi/4, _gamma=100, _n2=1):
  # chr_angle_vac - Cherenkov radiation angle after it leaves the target/crystal, crossed the medium boundary.
  # n1 - refractive index of the media
  # n2 - refractive index outside the media
  # _theta_cr - angel between the particle velocity and the normal to crystal surface (0 deg in normal case - velocity is perpendicular to crystal) [rad]
  # _gamma - Lorentz factor
  beta = np.sqrt(1 - 1/_gamma**2) # reduced velocity
  theta_ch = np.arccos(1 / (n * beta)) # ChR angle - angle between the particle velocity ang the ChR
  alpha2 = np.arcsin(np.sin(theta_ch - _theta_cr) * n1) # Angle between the ChR after crossing the boundary and the normal to the crystal/target surface [rad]
  chr_angle_vac = alpha2 + _theta_cr # Angle between the particle velocity and the ChR after crossing the boundary [rad]
  return chr_angle_vac

def refr_ind(lmbd):
  # Sellmeier equation - dependence of refractive index on wavelength of fused silica crystal
  # lmbd - wavelength [um]
  eps = 1 + 0.6961663 * lmbd**2 / (lmbd**2 - 0.0684043**2) + 0.4079426 * lmbd**2 / (lmbd**2 - 0.1162414**2) + 0.8974794 * lmbd**2 / (lmbd**2 - 9.896161**2) 
  n = np.sqrt(eps) # refractive index
  return n
