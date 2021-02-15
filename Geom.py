#------------------------------------------------------------------------------
# The Geom class contains the methods needed to construct the vascular network
# as a collection of nodes, vessel segments, and boundary conditions
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import scipy.io as sio

from math import cos
from math import sin
from math import radians

from Vessel import *
from LTETools import *
from Input import *
from ECell import *

#import matlab.engine

#------------------------------------------------------------------------------
class Geom:
    # Geom constructor 
    def __init__(self):
        # Current time
        self.time = 0
        
        # List of nodes (initially empty)       
        self.nodes = []
        
        # List of vessel segments (intially empty)
        self.vessels = []     
        
        # List of nodal degree (initially empty)                                                              
        self.node_degree = []                                                                

        # List of boundary conditions (intially empty)
        self.BCs = []                                                                        

        # List of bifurcation nodes (intially empty)
        self.bifur_nodes = []                                                                

        # List of bifurcation state (initially empty)
        self.bifur = []       

        self.total_cells = 0
        
        self.total_interc_iters = 0


    #---
    # Count the number of cells in the domain
    def count_cells(self):
        total_cell_num = 0
        
        for vess in self.vessels:
            total_cell_num += vess.num_cells
                
        self.total_cells = total_cell_num

    
    #---
    # Check that all cell positions are in range
    def check_cells(self):
        
        for vess in self.vessels:
            for cell in vess.cells:
                assert isinstance(cell.xi, float), 'xi is wrong type in cell {0}'.format(cell.ID)
                assert isinstance(cell.zeta, float), 'zeta is wrong type in cell {0}'.format(cell.ID)
                assert isinstance(cell.rho, float), 'rho is wrong type in cell {0}'.format(cell.ID)
                
                assert (cell.xi >= 0.) and (cell.xi <= 1.0), 'xi out of range in cell {0}'.format(cell.ID)
                assert (cell.zeta >= 0.) and (cell.zeta <= 1.0), 'zeta out of range  in cell {0}'.format(cell.ID)
                assert (cell.rho >= 0.), 'rho out of range in cell {0}'.format(cell.ID)
    
    
    #---
    # Calculate mean velocity (by cell)            
    def calc_vel_norm(self):
    
        mean_vel = 0.
    
        for vess in self.vessels:
            for cell in vess.cells:
                mean_vel += cell.vel.length()
    
        return mean_vel/self.total_cells
    
    
    #---
    # Calculate mean acceleration (by cell)             
    def calc_acc_norm(self):
    
        mean_acc = 0.
    
        for vess in self.vessels:
            for cell in vess.cells:
                mean_acc += cell.acc.length()
    
        return mean_acc/self.total_cells
    
    
    #---
    # Calculate max acceleration (by cell)             
    def calc_max_vel(self):
    
        max_vel = 0.
    
        for vess in self.vessels:
            for cell in vess.cells:
                if (cell.vel.length() > max_vel):
                    max_vel = cell.vel.length()
                    
        return max_vel
    
    
    #---
    # Calculate max acceleration (by cell)             
    def calc_max_acc(self):
    
        max_acc = 0.
    
        for vess in self.vessels:
            for cell in vess.cells:
                if (cell.acc.length() > max_acc):
                    max_acc = cell.acc.length()
                    
        return max_acc
    
    
    #---
    # Calculate mean change in vessel diameter (by vessel)  
    def calc_diam_norm(self):
        
        mean_dD = 0.
        
        for vess in self.vessels:
            mean_dD += vess.D - vess.Dold
                
        return mean_dD/len(self.vessels)    
    
    
#------------------------------------------------------------------------------
# Create the vessel network
def create_network(geom):

    if (i_sim_type == "A branch"):
            create_A_branch_network(geom)
        
    if (i_sim_type == "demo vess") or (i_sim_type == "PD curve"):
        create_demo_vess(geom)
    
    if (i_sim_type == "ideal cap bed"):
        load_ideal_cap_bed(geom)
        
    geom.count_cells()
    

#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_demo_vess(geom):
    
    # Create the first vessel section (lower left)
    geom.nodes.append((0, 0))

    for i in range(i_Nseg):
        node = (i_vess_length*(i+1), 0)
        geom.nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, i_num_cell, i_mu)
        geom.vessels.append(vess)
    
    # Print the nodal array     
    #Output.print_nodes(geom.nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
    
    find_vessel_neighbours(geom, BC_nodes)
    
#    # Assign the Dirchlet cell number condition to the outlet vessel
#    if (i_A_branch_dirchlet == True):
#        vessels[19].dirichlet = vessels[19].i_num_cells
    
    # Calculate the degree of each node
    for node in geom.nodes:
        geom.node_degree.append(0)
    
    for vess in geom.vessels:
        geom.node_degree[vess.n0] += 1
        geom.node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(geom.node_degree)):
        if (geom.node_degree[i] == 3):
            geom.bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(geom.bifur_nodes)):
        geom.bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in geom.vessels:
        vess.calc_length(geom.nodes[vess.n0], geom.nodes[vess.n1])
        vess.update_diameter()
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    geom.BCs.append((0, i_Pin))
    geom.BCs.append((i_Nn-1, i_Pout))
    

#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_A_branch_network(geom):    
    # Create the first vessel section (lower left)
    geom.nodes.append((0, 0))

    for i in range(5):
        node = (0, i_vess_length*(i+1))
        geom.nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, i_num_cell, i_mu)
        geom.vessels.append(vess)
    
    # Create the second vessel section (proximal horizontal branch)
    for i in range(6, 16):
        node = (i_vess_length*(i-5), 50.)
        geom.nodes.append(node)
        
        vess = Vessel(i, i-1, i, i_num_cell, i_mu)
        geom.vessels.append(vess)
    
    # Create the third vessel section (lower right)
    for i in range(16, 21):
        node = (100., 50. - i_vess_length*(i-15))
        geom.nodes.append(node)
        
        vess = Vessel(i, i-1, i, i_num_cell, i_mu)
        geom.vessels.append(vess)
    
    # Create the fourth vessel section (upper left)
    for i in range(21, 26):
        node = (0, 50. + i_vess_length*(i-20))
        geom.nodes.append(node)
        
        if (i == 21):
            vess = Vessel(i, 5, i, i_num_cell, i_mu)
        else:
            vess = Vessel(i, i-1, i, i_num_cell, i_mu)
            
        geom.vessels.append(vess)
    
    # Create the fifth vessel segment (distal horizontal branch)
    for i in range(26, 36):
        node = (i_vess_length*(i-25), 100.)
        geom.nodes.append(node)
        
        vess = Vessel(i, i-1, i, i_num_cell, i_mu)
        geom.vessels.append(vess)
    
    # Create the sixth vessel segment (lower right)    
    for i in range(36, 41):
        node = (100., 100. - i_vess_length*(i-35))
        if (i != 40):
            geom.nodes.append(node)
        
        if (i == 40):
            vess = Vessel(i, i-1, 15, i_num_cell, i_mu)
        else:
            vess = Vessel(i, i-1, i, i_num_cell, i_mu)
    
        geom.vessels.append(vess)
    
    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
    
    
    find_vessel_neighbours(geom, BC_nodes)
    
    # Assign the Dirchlet cell number condition to the outlet vessel
    if (i_A_branch_dirchlet == True):
        geom.vessels[19].dirichlet = geom.vessels[19].num_cells
    
    # Calculate the degree of each node
    for node in geom.nodes:
        geom.node_degree.append(0)
    
    for vess in geom.vessels:
        geom.node_degree[vess.n0] += 1
        geom.node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(geom.node_degree)):
        if (geom.node_degree[i] == 3):
            geom.bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(geom.bifur_nodes)):
        geom.bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in geom.vessels:
        vess.calc_length(geom.nodes[vess.n0], geom.nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    geom.BCs.append((0, i_Pin))
    geom.BCs.append((20, i_Pout))


#------------------------------------------------------------------------------
# Create the vessel network (load the ideal capillary bed from the MATLAB file)
def load_ideal_cap_bed(geom):
    # Generate the ideal capillary bed using the MATLAB script
#    mateng = matlab.engine.start_matlab()
#    mateng.generate_ideal_cap_bed(i_network_num_hc_long, i_network_num_hc_high, i_num_cell, nargout=0)
#    mateng.quit()
    
    # Load the network geometry from the MATLAB file
    ideal_cap_bed = sio.loadmat(i_in_cap_bed_filename, mat_dtype=True)
    vess_seg = ideal_cap_bed['vess_seg']

    filtered_vess_seg = []                                                      # Container of filtered segments
    rp = 8                                                                      # Rounding precision
        
    # Filter out any repeated vessel segments from the raw data
    for seg in vess_seg:
         
        yes_append = True
        
        # Optimize this if possible, must be a better way to do this
        for fvess_seg in filtered_vess_seg:
            
            if (round(seg[0], rp) == round(fvess_seg[0], rp) and round(seg[1], rp) == round(fvess_seg[1], rp) and round(seg[2], rp) == round(fvess_seg[2], rp) and round(seg[3], rp) == round(fvess_seg[3], rp)):
                yes_append = False
                
        if (yes_append == True):
            filtered_vess_seg.append(seg)   
    
    # Create the sub-vessel segments for the ABM based on the input segment size    
    for seg in filtered_vess_seg:
        x0 = float(seg[0])
        y0 = float(seg[1])
        x1 = float(seg[2])
        y1 = float(seg[3])
        
        seg_vect = Vect(x1 - x0, y1 - y0, 0.)
        seg_length = seg_vect.length()
        seg_vect.unit()
        
        num_sub_vess = int(round(seg_length/i_vess_length))
        seg_vess_length = seg_length/num_sub_vess
        
        for i in range(num_sub_vess):
            node1 = (round(x0 + i*(seg_vess_length*seg_vect.x), rp), round(y0 + i*(seg_vess_length*seg_vect.y), rp))
            node2 = (round(x0 + (i+1)*(seg_vess_length*seg_vect.x), rp), round(y0 + (i+1)*(seg_vess_length*seg_vect.y), rp))
            
            if (node1 not in geom.nodes):
                geom.nodes.append(node1)
                
            if (node2 not in geom.nodes):
                geom.nodes.append(node2)
                
            index1 = geom.nodes.index(node1)
            index2 = geom.nodes.index(node2)
            
            vess = Vessel(len(geom.vessels)+1, index1, index2, i_num_cell, i_mu)
                            
            vess.calc_length(node1, node2)
            vess.update_diameter()
            vess.update_conductance()
            
            geom.vessels.append(vess)
    
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in geom.vessels:        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
            
    # Give geom.vessels information about their neighbours
    BC_nodes = []
    
    find_vessel_neighbours(geom, BC_nodes)
    
    # Calculate the degree of each node
    for node in geom.nodes:
        geom.node_degree.append(0)
    
    for vess in geom.vessels:
        geom.node_degree[vess.n0] += 1
        geom.node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(geom.node_degree)):
        if (geom.node_degree[i] == 3):
            geom.bifur_nodes.append(i)
    
    # Initialize the bifurcation status array
    for i in range(len(geom.bifur_nodes)):
        geom.bifur.append(0)
        
    # Output network info
#    Output.print_nodes(geom.nodes)
#    Output.print_geom.vessels(geom.geom.vessels)
    
    # Prescribe boundary conditions
    #assert (len(BC_nodes) == len(PBC)), "Number of boundary nodes and number of prescribed pressures is not the same"
    
    x_min, y_min = min(geom.nodes, key=lambda item:item)
    x_max, y_max = max(geom.nodes, key=lambda item:item)   
    
    for i in range(len(geom.nodes)):
        if (geom.nodes[i][0] == x_min):
            geom.BCs.append((i, i_Part))
        elif (geom.nodes[i][0] == x_max):
            geom.BCs.append((i, i_Pvein))
        

    if (i_yes_dirichlet_BCs):
        for vess in geom.vessels:
            if (geom.nodes[vess.n1][0] == x_max):
                vess.dirichlet = vess.num_cells


#------------------------------------------------------------------------------
# Give vessels information on their upstream and downstream neighbours
def find_vessel_neighbours(geom, BC_nodes):

    for vess in geom.vessels:
        for vess2 in geom.vessels:
            if (vess.ID != vess2.ID):
                if ((vess.n0 == vess2.n0) or (vess.n0 == vess2.n1)):
                    vess.neigh0.append(vess2.ID)
                    
                if ((vess.n1 == vess2.n0) or (vess.n1 == vess2.n1)):
                    vess.neigh1.append(vess2.ID)
    
    # Apply periodic conditions at the boundary
    free_n0 = []
    free_n1 = []
    
    for vess in geom.vessels:
        if len(vess.neigh0) == 0:
            free_n0.append(vess.ID)
            
        if len(vess.neigh1) == 0:
            free_n1.append(vess.ID)
    
    if ((i_sim_type != "ideal sprout front") and (i_sim_type != "Y branch")):
        assert (len(free_n0) == len(free_n1)), "Periodic Boundary Error: Number of free ends is not equal"
        
        if (i_sim_type != "ideal cap bed"):
            for i in range(len(free_n0)):
                geom.vessels[free_n0[i]-1].neigh0.append(free_n1[i])
                geom.vessels[free_n1[i]-1].neigh1.append(free_n0[i])
        else:
            if (i_yes_flip_periodic_BCs == True):
                for i in range(len(free_n0)):
                    geom.vessels[free_n0[i]-1].neigh0.append(free_n1[len(free_n0)-1-i])
                    geom.vessels[free_n1[i]-1].neigh1.append(free_n0[len(free_n0)-1-i])
            else:
                for i in range(len(free_n0)):
                    geom.vessels[free_n0[i]-1].neigh0.append(free_n1[i])
                    geom.vessels[free_n1[i]-1].neigh1.append(free_n0[i])     
     
    elif (i_sim_type == "ideal sprout front"):
        geom.vessels[free_n0[0]-1].neigh0.append(free_n0[1])
        geom.vessels[free_n0[1]-1].neigh0.append(free_n0[0])
    elif (i_sim_type == "Y branch"):
        geom.vessels[free_n0[0]-1].neigh0.append(free_n1[0])
        geom.vessels[free_n0[0]-1].neigh0.append(free_n1[1])
        geom.vessels[free_n1[0]-1].neigh1.append(free_n0[0])
        geom.vessels[free_n1[1]-1].neigh1.append(free_n0[0])
        
    for n0 in free_n0:
        BC_nodes.append(geom.vessels[n0-1].n0)
        
    for n1 in free_n1:
        BC_nodes.append(geom.vessels[n1-1].n1)
    
    BC_nodes.sort()