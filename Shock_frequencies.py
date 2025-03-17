import numpy as np
import astropy.constants as const
from scipy.optimize import fsolve
import astropy.units as u
import plasmapy 
import plasmapy.particles as pp

me = (const.m_e.value) * (u.kg) # Electron mass in CGS (g)
e = 4.8032e-10  *u.cm**(3/2) * u.g**(1/2) * u.s**(-1) # Electron charge in CGS (esu)
c = 3e10 *u.cm/u.s


v0 = (1e7 * u.cm/u.s).to( u.m/u.s)
Te = 15 * u.eV
Ti = Te
ne = (7e19 /u.cm**3).to(1/u.m**3)



He = pp.Particle( 2, mass_numb=4, Z= 2)
N2 = pp.Particle( 7, mass_numb=14, Z= 7)
Proton =pp.proton
Ions = He #Choose here what kind of ions you have. use also
# Ions = CustomParticle( user_A *proton.mass, charge =  user_Z * proton.charge.value * u.C, symbol ="customHe+")
 
electrons = pp.electron
Z = Ions.atomic_number
mi = (Ions.mass.value) * (u.kg)
gamma = 5/3  



#Input velocity in simulation (reflective wall)
v = v0

def equation(v):
    v = v * u.m/u.s
    
    Mach_q = ((v**2 * mi ).to(u.J) / ( gamma * (Z+1) * Ti.to(u.J) ))
    
    # eta calculation
    numerator_eta = (gamma + 1) * Mach_q
    denominator_eta = (gamma - 1) * Mach_q + 2
    eta = numerator_eta / denominator_eta
    
    # v0 equation
    lhs = v0
    rhs = v * (1 - 1 / eta)
    return lhs - rhs

#v1 is the upstream value of velocity in SHOCK frame. It allows me to calculate the true Mach number.

v1_solution = fsolve(equation, (v.to(u.m/u.s)).value) 


v1 =  v1_solution[0] * u.m/u.s
#Calculating jump conditions, 2 refers to downstream in shock frame ref.

Mach_q = ((v1**2 * mi ).to(u.J) / ( gamma * (Z+1) * Ti.to(u.J) ))
print(Mach_q**0.5)
Mach_q = ((v**2 * mi ).to(u.J) / ( gamma * (Z+1) * Ti.to(u.J) ))
theta = ((gamma-1)*Mach_q + 2)*(2*gamma*Mach_q - gamma +1 )/(gamma+1)**2/Mach_q
eta = (gamma +1)*Mach_q/ ( 2 + (gamma-1)*Mach_q)


ne2 = ne*eta
Te2 = Te*theta
v2 = v1/eta
ni = ne/Z
ni2 = ne2/Z



print(v0)
print(v1)


#Mean free path calculation of ions at the shock front colliding with the colder background
nu_ii = plasmapy.formulary.SingleParticleCollisionFrequencies(Ions, Ions, v_drift=v1, T_b =Te, n_b = ni, Coulomb_log=10)
nu_ie = plasmapy.formulary.SingleParticleCollisionFrequencies(Ions, electrons, v_drift=v1, T_b =Te, n_b = ni, Coulomb_log=10)

nu = nu_ii.Lorentz_collision_frequency + nu_ie.Lorentz_collision_frequency
v_thermal = plasmapy.formulary.thermal_speed( Te2, Ions)
vth = np.sqrt(Te2.to(u.J)/me.to(u.kg)/1836/2).to(u.m/u.s)  #just to check

tau = (1/nu).to(u.s)
print(v_thermal, vth)
mfp_ii = (v_thermal/nu).to(u.m)
print('mean free path: ' + f'{mfp_ii:1.3E}')


unit_length = mfp_ii
unit_time = tau

print('units.length = ' f'{unit_length:1.3E}')
print('units.time = ' f'{unit_time:1.3E}')
print('units.number_density = '  f'{ne:1.3E}') # num/m^3
print('units.temperature    = ' f'{Te:1.3E}')



