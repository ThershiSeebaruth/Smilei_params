#-----------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI -- REFLECTIVE WALL CASE
# ----------------------------------------------------------------------

import math as m
import numpy as np
#
q  = 4.80326e-10       #statCoulomb
c  = 3.*1e10           #cm/s
me = 9.10938356*1e-28  #g
mi = 1836  #me/me = 1, mi/me = 1836


cost = 1.380649e-23 #/ 8.617333262e-5
## INPUT ##



# Eon plasma frequency in 1/s
v0 = 1.5e7 #velocity in cm/s --> 300 km/s 
v0_c = v0/c

## ------- PISTON PLASMA -------

TeP_eV = 45 # temperature of electrons in eV
TiP_eV = 45 # temperature of ions in eV 
TeP_mec2 = TeP_eV/(511*1e3) # temperature nondimensionalization
TiP_mec2 = TiP_eV/(511*1e3) # temperature nondimensionalization

ZP = 2  #fully ionized He -- see FLASH sim.

cost = 1.380649e-23 #/ 8.617333262e-5

#ne = 10 mbar / (kb* 290 K)  = ion density -->
#eon density for reference= ion * ZP * density jump condition (max 4)
p = 35 #mbar
p = p * 100 #Pascal
T_chamber = 290 #K
n0 = p/cost/T_chamber/1e6 * ZP * 4 

neP = n0/4 #This actually has a profile, varies from 1e17 to 1e20 cm-3 
niP = neP/ZP/4

neP_n0 = neP/n0 #density nondimensionalization
niP_n0 = niP/n0 #density nondimensionalization

miP = mi*4 



#--------Parameters for Input file----------------------------------------------------------------------
wpe = m.sqrt(4.*m.pi*n0*q**2/me)            # Eon plasma frequency in 1/s
c_wpe = c/wpe                               # Electron skin depth is reference length for Smilei), in cm.
inv_ld2 = (1/7.43e2)**2*(n0/TeP_eV)*(1+ZP)  # 1/lambda_D^2
deb_len = 1/np.sqrt(inv_ld2)                # Debye length

#------------- Case 1 -----------------------------------------------------------------------------------



k = 2.8  ## Should be < 3
dx = k*deb_len/c_wpe


#--------------------------------------------------------------------------------------------------------

cell_per_patch = 8 #min 6
num_patches = 2**11 #power of 2 (prioritize large number of patches?)

#--------------------------------------------------------------------------------------------------------

Nx= num_patches*cell_per_patch 
Lx = dx*Nx 
Nppc = int(k/2*250) #350 #?
#---------------------------------------------------------------------------------------------------------

# Normalization time in ns
dt = 0.9*dx  
Tsim_ns = 2 #Simulation time in ns
t_sim = int(Tsim_ns*wpe/1e9) #nondimensionalization 
Nt = t_sim/dt
debug_time = t_sim/100

t_sim = 50*dt #Testing parameter	
debug_time = 48


Main(
    geometry = "1Dcartesian",
    
    number_of_patches = [num_patches],
    
    interpolation_order = 2,

    timestep = dt,  
    
    simulation_time = t_sim,
    
    #time_fields_frozen = 5*t_sim,
    
    cell_length = [dx],  
    grid_length = [Lx],   
    
    EM_boundary_conditions = [ ["reflective", "silver-muller"],
    ],
    
    solve_poisson = True, #True or False?
    print_every  = int(t_sim/50.), 
    reference_angular_frequency_SI = wpe,
)

Species(
    name = "eonP",
    position_initialization = "random",
    momentum_initialization = "mj",
    particles_per_cell= Nppc,
    mass = 1.0,
    charge =-1.0,
    number_density = trapezoidal(neP_n0, xvacuum=2*dx, xplateau=16378*dx, xslope1=2*dx, xslope2=2*dx  ),
    mean_velocity = [-v0_c, 0., 0.], #need to change the values of the velocity 1
    temperature = [TeP_mec2], 
    boundary_conditions = [[ "reflective", "remove"],
    ],
)

Species(
    name = "ionP",
    position_initialization = "random",
    momentum_initialization = "mj",
    particles_per_cell= Nppc,
    mass = miP,
    charge = ZP,
    number_density = trapezoidal(niP_n0,  xvacuum=2*dx, xplateau=16378*dx, xslope1=2*dx, xslope2=2*dx  ),
    mean_velocity = [-v0_c, 0., 0.], #need to change the values of the velocity 2
    temperature = [TeP_mec2], 
    boundary_conditions = [["reflective", "remove"],
    ],
)



Collisions( #this corresponds to collisions between electrons and electrons 
    species1 = ["ionP"],
    species2 = ["ionP"],
    coulomb_log = 0,
    debug_every = debug_time*1,
)

Collisions( #this corresponds to collisions between electrons and electrons 
    species1 = ["eonP"],
    species2 = ["eonP"],
    coulomb_log = 0,
    debug_every = debug_time*1,
)

Collisions( #this corresponds to collisions between electrons and electrons 
    species1 = ["eonP"],
    species2 = ["ionP"],
    coulomb_log = 0,
    debug_every = debug_time*1,
)



DiagScalar(
	every=debug_time,
	vars = ["Utot", "Ukin", "Uelm", "Uexp", "Ubal", "Ukin_bnd", "Uelm_bnd", "Ukin_ionP", "Ukin_eonP"],
	precision =10,
)


# ---------- Ion and Electron number densities -------------
DiagParticleBinning( #(0)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)

DiagParticleBinning( #(1)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)



# ------------- Ion and Electron phase spaces -------------------
DiagParticleBinning( #(2)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
                 ["vx",   -0.5,   0.5, 500] #1024 is just a ssufficiently large number? 
        ]
)

DiagParticleBinning( #(3)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
                 ["vx",   -5*v0_c,    5*v0_c, 500] #1024 is just a ssufficiently large number? 
        ]
)



DiagProbe(
    every = debug_time,
    origin = [0],
    corners = [[Lx]],
    number = [Nx],

    fields = ['Ex','Ey']


)



DiagParticleBinning( #(4)
        deposited_quantity = "weight_vx_px",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx) ],
        ]
)

DiagParticleBinning( #(5)
        deposited_quantity = "weight_vx_px",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx,  int(Lx/dx)],
        ]
)

DiagParticleBinning( #(6)
        deposited_quantity = "weight_vx_py",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx) ],
        ]
)

DiagParticleBinning( #(7)
        deposited_quantity = "weight_vx_py",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)

DiagParticleBinning( #(8)
        deposited_quantity = "weight_vx_pz",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx) ],
        ]
)

DiagParticleBinning( #(9)
        deposited_quantity = "weight_vx_pz",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)




DiagParticleBinning( #(10)
        deposited_quantity = "weight_vy_py",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx) ],
        ]
)

DiagParticleBinning( #(11)
        deposited_quantity = "weight_vy_py",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)

DiagParticleBinning( #(12)
        deposited_quantity = "weight_vy_pz",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx) ],
        ]
)

DiagParticleBinning( #(13)
        deposited_quantity = "weight_vy_pz",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)


DiagParticleBinning( #(14)
        deposited_quantity = "weight_vz_pz",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx) ],
        ]
)

DiagParticleBinning( #(15)
        deposited_quantity = "weight_vz_pz",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, int(Lx/dx)],
        ]
)

