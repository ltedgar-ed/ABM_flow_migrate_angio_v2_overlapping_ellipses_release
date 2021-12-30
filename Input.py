#------------------------------------------------------------------------------
# The Input class contains all the input parameters, can be imported into any module
# that requires these parameters using import *
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import sys

from math import pi

from random_seed_numbers import *

#------------------------------------------------------------------------------
## Input parameters

i_on_linux = False

# Uncomment if running from bash script in the command line
#i_run = int(sys.argv[1])
i_run = 16

# Uncomment if running as stand alone in Spyder
i_branch_alpha = 0.5


i_bifurcation_rule = 1

i_rand_init = True

# Random seed generated from rand.randint(1,1000000000)
i_rseed = rseeds[i_run-1]

# Type of simulation to run
#i_sim_type = "demo vess"
#i_sim_type = "A branch"
i_sim_type = "ideal cap bed"
#i_sim_type = "PD curve"

i_fsi_on = False
i_Pext = 0*12.96
i_fsi_linear_PD = False

# Time step size
i_EndTime = 96
i_poe_dt = 0.5

# Number of time steps
i_Nt = int(i_EndTime/i_poe_dt)

## POE force balance parameters
# Ellipse properties
i_poe_Ain = 5.
i_poe_Bin = 2.5
i_poe_Aout = 1.5*i_poe_Ain
i_poe_Bout = 1.5*i_poe_Bin

# Dampening paramter (kg/hr)
i_poe_eta = 3.

# Force paramters (N*)
i_poe_kmig = 3.*i_poe_eta
i_poe_krep = 1.*i_poe_eta
i_poe_katt = (1/3)*i_poe_eta
#i_poe_katt = 0

i_find_acc_free = True
i_del_vel_min = 1e-2
    
# Dynamic viscosity of blood (Pa*-hr)
i_mu = 1.26 *1e-5

# EC migration speed (um/hr)
i_EC_speed = 3. 

# Width of each cell (um)                                             
i_cell_size = 5. 


# Create the animation of the simulation
i_plot_animation = False                                                          # Turn on the MATLAB plotting at the end of the simuation


# Inlet flow condition
i_inlet_flow_on = False


# Polarization re - alignment weights: w1 + w2 + w3 = 1
i_w2 = 1.0                                                                        # Polarization weight - flow component
i_w3 = 0.0                                                                        # Polarization weight - random walk component
i_w1 = 1 - i_w2 - i_w3                                                            # Polarization weight - persistence component

i_pol_subpop = True
i_subpop = 0.5
i_polshift_label = "_polsubpop_" + str(i_subpop)

# -------------------------
# demo vess model parameters
if (i_sim_type == "demo vess"):
    # Length of each vessel segment/length of each cell (m)
    i_vess_length = 2*i_cell_size
    
    # Inlet pressure (Pa*)
    i_Pin = 100*12.96     
    
    # Outlet pressure (Pa*)
    i_Pout = 0*12.96                             
    
    # Initial number of cells in vessel
    i_num_cell = 4                                                                    
       
    # Number of vessel segments
    i_Nseg = 5 
        
    # Number of nodes
    i_Nn = i_Nseg + 1
    
# -------------------------
# A branch model parameters
if (i_sim_type == "A branch"):
    # Length of each vessel segment/length of each cell (m)
    i_vess_length = 2*i_cell_size
    
    # Inlet pressure (Pa*) 
    i_Pin = 100*12.96       
    
    # Outlet pressure (Pa)
    i_Pout = 0*12.96                             
    
    # Initial number of cells in vessel
    i_num_cell = 4                                                                    
    
    # Number of nodes
    i_Nn = 40   
       
    # Number of vessel segments                                                               
    i_Nseg = 40

    # Flag for Dirchlet cell boundary conditions
    i_A_branch_dirchlet = False                                                                     


# -------------------------
# Ideal capillary bed parameters
if (i_sim_type == "ideal cap bed"):
    # Length of each vessel segment/length of each cell (m)
    i_vess_length = 10.
    
    # Width of each cell (m)                                             
    i_cell_size = 5.     
    
    # Pressure boundary conditions [Part, Pvein] (Pa)
    i_Part = 100*12.96
    i_Pvein = 0*12.96
    
    # Number of honeycombs long and high
    i_network_num_hc_long = 1
    i_network_num_hc_high = 1
    
    # Number of cells in the and capillaries
    i_num_cell = 4
    
    # Name of the MATLAB file containing the network geometry
    i_cap_bed_name = 'ideal_cap_bed_3x3'
    
    i_in_cap_bed_filename = i_cap_bed_name + '.mat'
    
    i_yes_dirichlet_BCs = False
    i_yes_flip_periodic_BCs = False
    

# -------------------------
# pressure-diameter curve parameters
if (i_sim_type == "PD curve"):
    # Length of each vessel segment/length of each cell (m)
    i_vess_length = 2*i_cell_size
    
    # Inlet pressure (Pa*)
    i_Pin = 100*12.96     
    
    # Outlet pressure (Pa*)
    i_Pout = i_Pin                            
    
    # Initial number of cells in vessel
    i_num_cell = 4                                                                    
       
    # Number of vessel segments
    i_Nseg = 5 
        
    # Number of nodes
    i_Nn = i_Nseg + 1
    
    i_w2 = 0.
    i_w3 = 0.
    i_w1 = 1.
    
    i_poe_kmig = 0.
    
    
# -------------------------
# Output file name
# Output file for demo vess simulation
if (i_sim_type == "demo vess"):
    i_out_filename = "ABM_output_demo_vess" + ("_Nt_" + str(i_Nt)) + ("_Ncell_" + str(i_num_cell))
    
# Output file for A branch simulation
if (i_sim_type == "A branch"):
    i_out_filename = "ABM_output_A_branch" + ("_Nt_" + str(i_Nt)) + ("_Ncell_" + str(i_num_cell))
  
    # Output file for A branch simulation
if (i_sim_type == "ideal cap bed"):
    i_out_filename = "ABM_output_" + i_cap_bed_name + ("_Nt_" + str(i_Nt)) + ("_" + str(i_num_cell) + "cells_")
  
if (i_sim_type == "PD curve"):
    i_out_filename = "ABM_output_PD_curve" + ("_Nt_" + str(i_Nt)) + ("_Ncell_" + str(i_num_cell)) + ("_Ptm_" + str(int(round((i_Pin - i_Pext)/12.96,0))))
    
# Add force transmission parameter info to the file name    
if (i_poe_krep != 0.):
    i_out_filename += "_kext_" + "{0:.2f}".format(i_poe_krep)
        
if (i_poe_katt != 0.):
    i_out_filename += "_kcoh_" + "{0:.2f}".format(i_poe_katt)

if (i_poe_kmig != 0.):
    i_out_filename += "_kmig_" + "{0:.2f}".format(i_poe_kmig)   

if (i_bifurcation_rule == 2):
    i_out_filename += "_BF2"

if (i_pol_subpop == True):
    i_out_filename += i_polshift_label
    
# Add run number to the output file name
i_out_filename += ("_run" + str(i_run))

# Linux version
if (i_on_linux == True):
    i_out_filename += ".mat"