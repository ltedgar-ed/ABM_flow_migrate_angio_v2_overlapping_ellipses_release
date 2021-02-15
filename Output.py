#------------------------------------------------------------------------------
# The Output class handles all printing, output, writing to files, etc
# *** NOTE: The Output class requires a folder inside the working directory,
# 'Simulation Data,' which it uses to store the output files
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import scipy.io as sio
import numpy as np
#import matlab.engine
import DGraph

from Input import *

import os
#------------------------------------------------------------------------------
# The out_file class stores all the information for writing the output file
class OutFile:
    #---
    # OutFile Constructor
    def __init__(self, filename, geom):
        # Name of the ouput file
        self.filename = filename
        
        # Container of nodal positions
        self.nodes = geom.nodes.copy()

        self.dgraph = DGraph.DGraph(len(self.nodes))        
        
        # Number of vessels
        v = len(geom.vessels)
        
        # Container of vessel connectivity
        self.vess_conn = np.zeros((v, 2))

        for vess in geom.vessels:
            self.vess_conn[vess.ID-1][0] = vess.n0
            self.vess_conn[vess.ID-1][1] = vess.n1
            self.dgraph.addEdge(vess.n0,vess.n1)
        
        # Container of nodal pressures
        self.nodal_pressures = np.zeros((len(geom.nodes), i_Nt+1))
        
        # Container of vessel diameter
        self.vess_diameter = np.zeros((v, i_Nt+1))
        
        # Container of vessel flow
        self.vess_flow = np.zeros((v, i_Nt+1))
        
        # Container of vessel wall shear stress
        self.vess_WSS = np.zeros((v, i_Nt+1))
        
        # Container of vessel cell number
        self.vess_num_cells = np.zeros((v, i_Nt+1))
        
        # Container of the state of ECs at each time point
        self.cells = np.zeros((i_Nt+1,),dtype=np.object)
        
        # Container for the boundary conditions
        self.BCs = geom.BCs
        
        # Containter for the bifurcation nodes
        self.bifur_nodes = geom.bifur_nodes
        
        # Container for bifurcation status
        self.bifur_status = np.zeros((len(geom.bifur), i_Nt+1))
        
        # Initialize the containers
        for vess in geom.vessels:
            self.nodal_pressures[vess.n0][0] = vess.P0
            self.nodal_pressures[vess.n1][0] = vess.P1
                
            self.vess_diameter[vess.ID-1][0] = vess.D
            self.vess_flow[vess.ID-1][0] = vess.Q
            self.vess_WSS[vess.ID-1][0] = vess.tau
            self.vess_num_cells[vess.ID-1][0] = vess.num_cells
                
        for i in range(len(geom.bifur)):
            self.bifur_status[i][0] = geom.bifur[i]

        #self.save_cells(geom)
    
        self.input = {}
        self.save_input_params()
        
        self.interc_iters = np.zeros((1, i_Nt+1))
        
        
        
        
    #---
    # Save input parameters
    def save_input_params(self):
#        self.input['alpha'] = i_branch_alpha
        self.input['Nt'] = i_Nt
        self.input['Ain'] = i_poe_Ain
        self.input['Bin'] = i_poe_Bin
        self.input['Aout'] = i_poe_Aout
        self.input['Bout'] = i_poe_Bout
        self.input['eta'] = i_poe_eta
        self.input['krep'] = i_poe_krep
        self.input['katt'] = i_poe_katt
        self.input['kmig'] = i_poe_kmig
        self.input['dt'] = i_poe_dt
        self.input['mu'] = i_mu
        self.input['cell size'] = i_cell_size
        self.input['vess length'] = i_vess_length
        self.input['num cell'] = i_num_cell                                                                  
        #self.input['Nseg'] = i_Nseg 
        #self.input['Nn'] = i_Nn
        self.input['w1'] = i_w1
        self.input['w2'] = i_w2
        self.input['w3'] = i_w3
        self.input['Pext'] = i_Pext
        
    
    #---    
    # Save the vessel state to the containers
    def save_vessels(self, geom):
        time_step = geom.time
        
        for vess in geom.vessels:
            self.nodal_pressures[vess.n0][time_step] = vess.P0
            self.nodal_pressures[vess.n1][time_step] = vess.P1
                
            self.vess_diameter[vess.ID-1][time_step] = vess.D
            self.vess_flow[vess.ID-1][time_step] = vess.Q
            self.vess_WSS[vess.ID-1][time_step] = vess.tau
            self.vess_num_cells[vess.ID-1][time_step] = vess.num_cells
            

        for i in range(len(geom.bifur)):
            self.bifur_status[i][time_step] = geom.bifur[i]
            
        self.interc_iters[0][time_step] = geom.total_interc_iters
    
    
    #---        
    # Save the cell state to the containers
    def save_cells(self, geom):
        time_step = geom.time
        
        num_cells_total = 0
        
        for vess in geom.vessels:
            for cell in vess.cells:
                num_cells_total += 1
                
        curr_cell_state = np.zeros((num_cells_total, 13))
        
        cell_counter = 0
        
        for vess in geom.vessels:
            for cell in vess.cells:
                curr_cell_state[cell_counter][0] = cell.ID
                curr_cell_state[cell_counter][1] = cell.vessID
                curr_cell_state[cell_counter][2] = cell.xi
                curr_cell_state[cell_counter][3] = cell.zeta
                curr_cell_state[cell_counter][4] = cell.rho
                curr_cell_state[cell_counter][5] = cell.polarity.x
                curr_cell_state[cell_counter][6] = cell.polarity.y
                curr_cell_state[cell_counter][7] = cell.mig_net_force.x
                curr_cell_state[cell_counter][8] = cell.mig_net_force.y
                curr_cell_state[cell_counter][9] = cell.mig_vel.x
                curr_cell_state[cell_counter][10] = cell.mig_vel.y
                curr_cell_state[cell_counter][11] = cell.pos_overlap
                curr_cell_state[cell_counter][12] = cell.neg_overlap
                cell_counter += 1
        
        self.cells[time_step] = curr_cell_state
        
    
    #---    
    # Write the final output file once the simulation has ended    
    def write_output_file(self):        
        # Save the output of the simulation as a MATLAB .mat file in the Simulation Data directory
        curr_path = os.getcwd()
        #save_path = curr_path + "\Simulation Data"
        
        if (i_on_linux == False):
            save_path = curr_path + "\Simulation Data"
        else:
            save_path = curr_path + "/SimulationData"
            
        os.chdir(save_path)        
        sio.savemat(self.filename, {'nodes' : self.nodes, 'vess_conn' : self.vess_conn, 'nodal_pressures' : self.nodal_pressures, 'vess_diameter' : self.vess_diameter, 'vess_flow' : self.vess_flow, 'vess_WSS' : self.vess_WSS, 'vess_num_cells' : self.vess_num_cells, 'cells' : self.cells, 'PBCs' : self.BCs, 'bifur_nodes' : self.bifur_nodes, 'bifur' : self.bifur_status, 'input' : self.input, 'interc_iters' : self.interc_iters})
        os.chdir(curr_path)
        
        # Open the MATLAB plotting script to visualize the results
        if (i_plot_animation == True):
            print("Plotting animation...")
            
            mateng = matlab.engine.start_matlab()
            
            if (i_sim_type == "A branch"):
                mateng.plot_ABM_output_A_branch(self.filename, nargout=0)
            elif (i_sim_type == "demo vess"):
                mateng.plot_ABM_output_demo_vess(self.filename, nargout=0)
                
            mateng.quit()

        
        
#------------------------------------------------------------------------------
# Print the nodal positions to the console
def print_nodes(geom):
    # Print the nodal array     
    for node in geom.nodes:
        print("({0}, {1})".format(round(node[0], 1), round(node[1], 1)))


#------------------------------------------------------------------------------
# Print the vessels status to the console
def print_vessels(geom):
    # Print the vessels to the console
    print()
    
    for vess in geom.vessels:
        print(vess)
        

#------------------------------------------------------------------------------
# Print info on vessel neighbours to the console
def print_vess_neighbours(geom):
    # Print info on neighbours to console for verification            
    for vess in geom.vessels:
        print("Vessel {0} has upstream neighbors {1} and downstream neighbors {2}".format(vess.ID, vess.neigh0, vess.neigh1))


#------------------------------------------------------------------------------  
# Print information on the cells to the console
def print_cells(geom):
    # Print the cells to the console
    print()
    for vess in geom.vessels:
        for cell in vess.cells:
            print(cell)

#------------------------------------------------------------------------------
# Print the number of cells in each vessel to the consoles
def print_num_cells(geom):
    # Print the number of cells in each vessel to the console
    print()
    for vess in geom.vessels:
        print("Vessel {0} has {1} cells with cell number set as {2}".format(vess.ID, len(vess.cells), vess.num_cells))

    
def print_cells_net_force(geom):
    
    for vess in geom.vessels:
        for cell in vess.cells:
            print("Net force acting on cell {0} is ({1}, {2}, {3})".format(cell.ID, cell.net_force.x, cell.net_force.y, cell.net_force.z))
    
    
def calc_mean_diam(geom):
    
    mean_diam = 0.
    
    for vess in geom.vessels:
        mean_diam += vess.D
        
    mean_diam /= len(geom.vessels)
    
    print("Mean vessel diameter {0} um.".format(mean_diam))
    
    