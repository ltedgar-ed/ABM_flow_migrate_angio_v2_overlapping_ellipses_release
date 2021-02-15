#------------------------------------------------------------------------------
# The ECell class represents a migratory endothelial cell within the blood vessel wall
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import random as rand

from LTETools import *

#------------------------------------------------------------------------------
class ECell:
    # ECell constructor
    def __init__(self, ID =0, vessID =0, xi =0.5, zeta =0., polvect =Vect(1., 0., 0.,)):
        
        # Cell number identifier
        ID : int 
        self.ID = ID
        
        # Identifier for which Vessel the cell resides in
        vessID : int
        self.vessID = vessID
        
        # Natural coordinates indicating position within the Vessel segment
        xi : float
        self.xi = xi
        
        # Natural coordinates indicating circumfrential position as a percentage of the circumfrence
        zeta : float
        self.zeta = zeta
        
        # Natural coordinates indicating radial stretch
        rho : float
        self.rho = 1.
        
        # Remaining distance - distance cell has left to travel once it moves to a new segment
        rem_dist : Vect
        self.rem_dist = Vect()
        
        # Cell polarity vector
        polarity : Vect
        self.polarity = polvect
    
        # Migration idicator (+1 for upstream/against vessel unit vector, -1 for downstream/along vessel unit vector, 0 for not migrating)
        migrate : int
        self.migrate = 0
        
        # Net force
        net_force : Vect
        self.net_force = Vect()
        
        # Migration velocity
        vel : Vect
        self.vel = Vect()
        
        # Migration net force for recording
        mig_net_force : Vect
        self.mig_net_force = Vect()
        
        # Migration veloctiy for recording
        mig_vel : Vect
        self.mig_vel = Vect()
        
        # Migration velocity (old)
        old_vel : Vect
        self.old_vel = Vect()
        
        # Migration acceleration
        mig_acc : Vect
        self.acc = Vect()
        
        # POE overlap
        pos_overlap : float
        self.pos_overlap = 0.
        
        neg_overlap : float
        self.neg_overlap = 0.
        
        
    def __str__(self):
        #return "Cell {0} is in Vessel {1}; polarity = {2}; migrate = {3}".format(self.ID, self.vessID, round(self.polarity, 2), self.migrate)
        return "Cell {0} is in Vessel {1} with xi {2} and zeta {3}".format(self.ID, self.vessID, round(self.xi, 3),  round(self.zeta, 3))

    
#------------------------------------------------------------------------------       
# main program
if __name__ == "__main__":
    pass
