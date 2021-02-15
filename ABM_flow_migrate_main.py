#------------------------------------------------------------------------------
# An agent-based model of the flow-coupled migration of vascular endothelial cells
# during developmental angiogenesis
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import random as rand

from Input import *
from Vessel import *

import Geom
import Network
import Flow
import Output


#------------------------------------------------------------------------------
## Initialise the simulation
# Set the random seed
rand.seed(i_rseed)
print(i_rseed)

# Create the container class for the network geometry
geom = Geom.Geom()

# Construct the network
Geom.create_network(geom)

# Solve for initial flow
Flow.solve_for_flow(geom)

# Print out the time step
print("Migration time step {0} out of {1}...".format(0, i_Nt))

# Initial polarity alignment
Network.realign_polarity(geom)

# Calculate the initial force balance on cells
Network.cell_randomise_force(geom)

# Perform initial cell migration step
Network.cell_migration(geom, False, False, 0)
    
# Intial intercalation step
Network.intercalate(geom, True)

# Initialise the ouput file
out = Output.OutFile(i_out_filename, geom)

# Save the vessel and cell states for output
out.save_vessels(geom)
out.save_cells(geom)


#------------------------------------------------------------------------------
## Step through time
for t in range(i_Nt):
    # Set the current time
    geom.time = t+1
    
    # Print out the time step
    print("Migration time step {0} out of {1}...".format(geom.time, i_Nt))
    
    # Realign polarity
    Network.realign_polarity(geom)

    # Calculate force balance on cells
    Network.cell_force_balance(geom)
    
    # Perform cell migration
    Network.cell_migration(geom, True, False, 0)
   
    # Intercalate cell to correct spacing after migration
    Network.intercalate(geom, False)
    
    # Check cell positions and make sure nothing has out of range
    #Network.check_cells(geom)
    
    # Solve for flow in the new network configuration
    Flow.solve_for_flow(geom)
    
    # Save the vessel and cell states for output
    out.save_vessels(geom)
    out.save_cells(geom)
    
    
#------------------------------------------------------------------------------   
## Write the output files  
out.write_output_file()    


#------------------------------------------------------------------------------
## Terminate the simulation
print("Done!")


#------------------------------------------------------------------------------