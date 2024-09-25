#! /usr/bin/python3
import numpy as np

# reference conditions for the upper mantle
# depth_ref  = 30e3 # m
# visc_ref   = 1e21
# temp_ref   = 495
# depth_ref  = 20e3 # m
# visc_ref   = 1e20
# temp_ref   = 495
# strain_ref = 2.6e-16
midmantle_viscosity_jump = 25
# depth_ref  = 0e3 # m
visc_ref   = 1e19
temp_ref   = 600
depth_ref  = 0e3 # m




strain = 9.e-13
cr_yr = 0.1 #m/yr
yr = 365*24*60*60
cr = cr_yr/yr #m/s
shear_d_ref = cr/strain #m
shear_d = np.array([shear_d_ref/2, shear_d_ref, shear_d_ref*2])

strain_ref = cr/shear_d
print("Shear zone thickness = ", shear_d/1.e3, "km")
print("strain rate = ", strain_ref, "s-1")


# rheological parameters (e.g. Hirth and Kohlstedt, 2003)
Edisl = 125e3; Ediff = 0
Vdisl = 0;  Vdiff = 0; Vdiff_lowermant = 0
n = 4;        R = 8.314
adiabat = 0.3; # K/km
temp_ref   = temp_ref + (depth_ref * 1e-3 * adiabat)
press_ref  = 2900. * 9.81 * depth_ref
# press_ref = 0


# compute dislocation creep prefactor using:
# eta = 1/2*Aeq^(-1/n) * strain_rate^((1-n)/n) * exp((E+PV)/nRT)
# Aeq   = ((2 *eta)^(-n)) / (strain_rate^(n-1) * exp(-(E+PV)/RT))
# Adisl = ((2*visc_ref)**(-1. * n)) / (  strain_ref**(n-1) * np.exp(-(Edisl+press_ref*Vdisl)/(R*temp_ref))  )
for i in range(len(strain_ref)):
    Adisl = (0.5*np.exp(Edisl/(n*R*temp_ref))*strain_ref[i]**((1-n)/n)/visc_ref)**(n)
    visc_check_disl = (1/2)*Adisl**(-1/n) * strain_ref[i]**((1-n)/n) * np.exp((Edisl+press_ref*Vdisl)/(n*R*temp_ref))
    print("Strain rate= %e - Dislocation prefactor = %e. Check: %e = %e" % (strain_ref[i], Adisl,visc_ref,visc_check_disl))

