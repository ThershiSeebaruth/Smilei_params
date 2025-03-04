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


scale = 1
# Eon plasma frequency in 1/s
v0 = 1.0e7/np.sqrt(scale) #velocity in cm/s --> 300 km/s 
v0_c = v0/c

## ------- PISTON PLASMA -------

TeP_eV = 15/scale # temperature of electrons in eV
TiP_eV = 15/scale # temperature of ions in eV 
TeP_mec2 = TeP_eV/(511*1e3) # temperature nondimensionalization
TiP_mec2 = TiP_eV/(511*1e3) # temperature nondimensionalization

ZP = 2  #fully ionized He -- see FLASH sim.

cost = 1.380649e-23 #/ 8.617333262e-5

#ne = 10 mbar / (kb* 290 K)  = ion density -->
#eon density for reference= ion * ZP * density jump condition (max 4)
#p = 35 #mbar
#p = p * 100 #Pascal
#T_chamber = 290 #K
#n0 = p/cost/T_chamber/1e6 * ZP * 4 
n0 = 4e18 #downstream value
neP = n0/4 #This actually has a profile, varies from 1e17 to 1e20 cm-3 
niP = neP/ZP

neP_n0 = neP/n0 #density nondimensionalization
niP_n0 = niP/n0 #density nondimensionalization

miP = mi*4 



#--------Parameters for Input file----------------------------------------------------------------------
wpe = m.sqrt(4.*m.pi*n0*q**2/me)            # Eon plasma frequency in 1/s
c_wpe = c/wpe                               # Electron skin depth is reference length for Smilei), in cm.
#inv_ld2 = (1/7.43e2)**2*(n0/TeP_eV)*(1+ZP)  # 1/lambda_D^2
#deb_len = 1/np.sqrt(inv_ld2)                # Debye length

#------------- Case 1 -----------------------------------------------------------------------------------



#k = 30  ## Should be < 3
#dx = k*deb_len/c_wpe
dx = 0.143*6

#--------------------------------------------------------------------------------------------------------

cell_per_patch = 6 #min 6
num_patches = 2**9 #power of 2 (prioritize large number of patches?)

#--------------------------------------------------------------------------------------------------------

Nx= num_patches*cell_per_patch 
Lx = dx*Nx 
Nppc =20000 #int(k/2*250) #350 #?
#---------------------------------------------------------------------------------------------------------

# Normalization time in ns
dt = 0.95*dx  
Tsim_ns = 4 #Simulation time in ns
t_sim = int(Tsim_ns*wpe/1e9) #nondimensionalization 
Nt = t_sim/dt
debug_time =int( Nt/100/Tsim_ns*2)

#t_sim = 50*dt #Testing parameter	
#debug_time = 1


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
    
    solve_poisson = False, #True or False?
    print_every  = int(Nt/100), 
    reference_angular_frequency_SI = wpe,
)

Species(
    name = "eonP",
    position_initialization = "random",
    momentum_initialization = "mj",
    particles_per_cell= Nppc,
    mass = 1.0,
    charge =-1.0,
    number_density = trapezoidal(neP_n0, xvacuum=2*dx, xplateau=3038*dx, xslope1=30*dx, xslope2=2*dx  ),
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
    number_density = trapezoidal(niP_n0, xvacuum=2*dx, xplateau=3038*dx, xslope1=30*dx, xslope2=2*dx  ),
    mean_velocity = [-v0_c, 0., 0.], #need to change the values of the velocity 1
    temperature = [TeP_mec2],
    boundary_conditions = [[ "reflective", "remove"],
    ],
)

Species(
    name = "ion_track",
    position_initialization = "random",
    momentum_initialization = "mj",
    particles_per_cell= int(Nppc/100),
    mass = miP,
    charge = ZP,
    number_density = trapezoidal(niP_n0, xvacuum=2*dx, xplateau=3038*dx, xslope1=30*dx, xslope2=2*dx  ),
    mean_velocity = [-v0_c, 0., 0.], #need to change the values of the velocity 1
    temperature = [TeP_mec2],
    boundary_conditions = [[ "reflective", "remove"],
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

Collisions( #this corresponds to collisions between electrons and electrons 
    species1 = ["eonP"],
    species2 = ["ion_track"],
    coulomb_log = 0,
    debug_every = debug_time*1,
)
Collisions( #this corresponds to collisions between electrons and electrons 
    species1 = ["ionP"],
    species2 = ["ion_track"],
    coulomb_log = 0,
    debug_every = debug_time*1,
)
Collisions( #this corresponds to collisions between electrons and electrons 
    species1 = ["ion_track"],
    species2 = ["ion_track"],
    coulomb_log = 0,
    debug_every = debug_time*1,
)

Checkpoints(
    #restart_dir ="/ccc/scratch/cont003/gen7678/seebarut/Hua_e18",
    dump_minutes = 1420,
    exit_after_dump = True,
    keep_n_dumps = 1,
)



DiagScalar(
	every=debug_time,
	vars = ["Utot", "Ukin", "Uelm", "Uexp", "Ubal", "Ukin_bnd", "Uelm_bnd", "Ukin_ionP", "Ukin_eonP"],
	precision =10,
)

debug_x = int(Lx/dx/4)
# ---------- Ion and Electron number densities -------------
DiagParticleBinning( #(0)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)

DiagParticleBinning( #(1)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)



# ------------- Ion and Electron phase spaces -------------------
DiagParticleBinning( #(2)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
                 ["vx",   -0.5,   0.5, 400] #1024 is just a ssufficiently large number? 
        ]
)

DiagParticleBinning( #(3)
        deposited_quantity = "weight",
        every =debug_time, 
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x],
                 ["vx",   -3*v0_c,    3*v0_c, 400] #1024 is just a ssufficiently large number? 
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
        axes = [ ["x",    0,    Lx, debug_x ],
        ]
)

DiagParticleBinning( #(5)
        deposited_quantity = "weight_vx_px",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx,  debug_x],
        ]
)

DiagParticleBinning( #(6)
        deposited_quantity = "weight_vx_py",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x ],
        ]
)

DiagParticleBinning( #(7)
        deposited_quantity = "weight_vx_py",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)

DiagParticleBinning( #(8)
        deposited_quantity = "weight_vx_pz",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x ],
        ]
)

DiagParticleBinning( #(9)
        deposited_quantity = "weight_vx_pz",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)




DiagParticleBinning( #(10)
        deposited_quantity = "weight_vy_py",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x ],
        ]
)

DiagParticleBinning( #(11)
        deposited_quantity = "weight_vy_py",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)

DiagParticleBinning( #(12)
        deposited_quantity = "weight_vy_pz",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x ],
        ]
)

DiagParticleBinning( #(13)
        deposited_quantity = "weight_vy_pz",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)


DiagParticleBinning( #(14)
        deposited_quantity = "weight_vz_pz",
        every =debug_time,
        species = ["ionP"],
        axes = [ ["x",    0,    Lx, debug_x ],
        ]
)

DiagParticleBinning( #(15)
        deposited_quantity = "weight_vz_pz",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
        ]
)

DiagParticleBinning( #(16)
        deposited_quantity = "weight",
        every =debug_time,
        species = ["eonP"],
        axes = [ ["x",    0,    Lx, debug_x],
	["ekin", 0, 100, 1000, "logscale", "edge_inclusive"]
        ]
)

