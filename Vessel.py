#------------------------------------------------------------------------------
# The Vessel class represents a segment of blood vessel, componsed of endothelial
# cells, and perfused with blood flow (Newtonian fluid, Poiseuille flow)
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

from math import sqrt
from math import pi
from math import exp
import numpy as np
#import random as rand

from Input import *
from LTETools import *
from ECell import *

#------------------------------------------------------------------------------
class Vessel:
    # Vessel constructor
    def __init__(self, ID=0, node1=0, node2=0, cell_num =0, mu = 1.):
         
        # Vessel number identifier
        ID : int
        self.ID = ID
        
        # Vessel type identifier
        type : str
        self.type = "capillary"
        
        # Node number identifiers (from master node list)
        n0 : int
        n1 : int
        
        self.n0 = node1                                                         
        self.n1 = node2                                                         
        
        # Vessel length
        L  : float        
        self.L = 1.
        
        # Vessel unit vector
        unit : Vect
        self.unit = Vect()
        
        # Number of cells residing in the vessel
        num_cells : int
        self.num_cells = cell_num
        self.old_num_cells = cell_num
        
        # Vessel diameter
        D  : float
        self.D = cell_num*i_cell_size/pi
        
        self.Dold = 0.
        
        # Dynamic viscosity of the fluid
        mu : float
        self.mu = mu
        
        # Vessel conductance
        G  : float
        self.G = 0.
        
        # Pressure at first node
        P0 : float
        self.P0 = 0.
        
        # Pressure at the second node
        P1 : float
        self.P1 = 0.
        
        # Flow through the vessel
        Q  : float 
        self.Q = 0.
        
        # List of resident endothelial cells
        cells : list[ECell]
        self.cells = []
        
        for i in range(self.num_cells):
            # Seed initial agent positions
            if (i_rand_init == True):
                xi_init = rand.uniform(0,1)
                zeta_init = rand.uniform(0,1)
            else:
                xi_init = 0.5
                zeta_init = (1/self.num_cells/2) + i*(1/self.num_cells)

            pol_vect = Vect(1., 0., 0.)
            
            self.cells.append(ECell(0, self.ID, xi_init, zeta_init, pol_vect))
            
#        for i in range(1,self.num_cells):
#            self.cells[i].zeta = 0.5 + 1e-3 + (i-1)*(0.5/(self.num_cells-1))
            
        # Vessel shearing parameter (equivalant of conductance but for shear stress)
        H : float
        self. H = 0.
        
        # Wall Shear Stress within the vessel segment
        tau : float
        self.tau = 0.
        
        # Information on upstream neighbours
        neigh0 : list[int]
        self.neigh0 = []
        
        # Information on downstream neighbours
        neigh1 : list[int]
        self.neigh1 = []
        
        # Indicator of a Dirichlet condition
        dirichlet : int
        self.dirichlet = 0
        
    
    # Calculate the length of the vessel based on node positions
    def calc_length(self, node_pos1, node_pos2):
        x0 = node_pos1[0]
        y0 = node_pos1[1]
        x1 = node_pos2[0]
        y1 = node_pos2[1]
        
        self.L = sqrt((x1 - x0)**2 + (y1 - y0)**2)
        
        self.unit.x = x1 - x0
        self.unit.y = y1 - y0
        self.unit.unit()
        
    
    # Update diameter
    def update_diameter(self):
        self.Dold = self.D
        
        D0 = self.num_cells*i_cell_size/pi
        
        mean_rho = 0
        
        for cell in self.cells:
            mean_rho += cell.rho
        
        try:
            mean_rho /= self.num_cells
        except:
            mean_rho = 0.        

        mean_rho = 1.0
        self.D = D0*mean_rho
        
        
    # Update Vessel conductance based on diameter, length, and dynamic viscosity 
    def update_conductance(self):
        if (self.D != 0.):
            self.G = pi*(self.D**4.)/(128.*self.mu*self.L)
            self.H = (32*self.mu)/(pi*self.D**3.)
        else:
            self.G = 1e-3
            self.H = 0.


    # Calculate flow through the vessel based off of the pressure difference across the vessel
    def calc_flow(self):
        self.Q = -self.G*(self.P1 - self.P0)
               
        self.tau = self.H*abs(self.Q)
      
        
    # String conversion function for printing    
    def __str__(self):
        return "Vessel {0} is a {1}: node1 = {2}; node2 = {3}; L = {4:.2f} um; num_cells = {5}; D = {6:.2f} um; P0 = {7:.2f} Pa; P1 = {8:.2f} Pa; Q = {9:.2f} nL/hr; tau = {10:.2f} Pa".format(self.ID, self.type, self.n0, self.n1, round(self.L, 2), self.num_cells, round(self.D, 2), round(self.P0/12.96, 2), round(self.P1/12.96, 2), round(self.Q/1e6, 2), round(self.tau/12.96, 2))


#------------------------------------------------------------------------------       
# main program
if __name__ == "__main__":
    pass

