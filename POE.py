#------------------------------------------------------------------------------
# The POE class contains the model of force transmission between ECs
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2020

from math import sqrt
from math import pi
from math import sin
from math import cos

from Vessel import *
from LTETools import *
from Input import *
from ECell import *

#------------------------------------------------------------------------------    
# Draw distance vector from cell i to cell j residing in the same vessel
def distance_between_cells_same_vess(cell_i, cell_j, Z, R):
    z_i = Z*cell_i.xi
    theta_i = 2*pi*cell_i.zeta

    z_j = Z*cell_j.xi
    theta_j = 2*pi*cell_j.zeta

    dz = z_j - z_i

    dtheta = theta_j - theta_i

    if (dtheta > pi):
        dtheta = dtheta - 2*pi
        
    if (dtheta < -pi):
        dtheta = dtheta + 2*pi
        
    dc = R*dtheta
    
    r_vect = Vect(dz, dc, 0.)
    
    return r_vect


#------------------------------------------------------------------------------ 
# Draw distance vector from cell i to cell j residing in the upstream neighbour vessel
def distance_between_cell_and_upstream_neigh(cell_i, Z1, R1, cell_j, Z2, R2):
    z_i = Z1*cell_i.xi
    theta_i = 2*pi*cell_i.zeta
    
    z_j = Z2*cell_j.xi
    theta_j = 2*pi*cell_j.zeta
    
    dz = -z_i - (Z2 - z_j)
    
    dtheta = theta_j - theta_i

    if (dtheta > pi):
        dtheta = dtheta - 2*pi
        
    if (dtheta < -pi):
        dtheta = dtheta + 2*pi
        
    dc = (R2 - R1)*theta_i + R2*dtheta
    
    r_vect = Vect(dz, dc, 0.)
    
    return r_vect


#------------------------------------------------------------------------------
# Draw distance vector from cell i to cell j residing in the downstream neighbour vessel
def distance_between_cell_and_downstream_neigh(cell_i, Z1, R1, cell_j, Z2, R2):
    z_i = Z1*cell_i.xi
    theta_i = 2*pi*cell_i.zeta
    
    z_j = Z2*cell_j.xi
    theta_j = 2*pi*cell_j.zeta
    
    dz = (Z1 - z_i) + z_j
    
    dtheta = theta_j - theta_i

    if (dtheta > pi):
        dtheta = dtheta - 2*pi
        
    if (dtheta < -pi):
        dtheta = dtheta + 2*pi
        
    dc = (R2 - R1)*theta_i + R2*dtheta
    
    r_vect = Vect(dz, dc, 0.)
    
    return r_vect


#------------------------------------------------------------------------------
# Calculate the force due to overlap
def calc_overlap_force(cell_i, cell_j, r_vect):
    
    L_ij = r_vect.length()
    r_vect_hat = r_vect.copy()
    r_vect_hat.unit()
    
    if (L_ij == 0.):
        return Vect(0., 0., 0.)
    
    x_vect = Vect(1., 0., 0.)
    
    alpha_i = findangle2D(x_vect, cell_i.polarity)
    alpha_j = findangle2D(x_vect, cell_j.polarity)
    
    phi_ij = findangle2D(x_vect, r_vect) - alpha_i
    phi_ji = findangle2D(x_vect, -r_vect) - alpha_j
    
    Ain = i_poe_Ain
    Bin = i_poe_Bin
    Aout = i_poe_Aout
    Bout = i_poe_Bout
    
    Ri_in = (Ain*Bin)/sqrt((Bin*cos(phi_ij))**2 + (Ain*sin(phi_ij))**2)
    Ri_out = (Aout*Bout)/sqrt((Bout*cos(phi_ij))**2 + (Aout*sin(phi_ij))**2)
    Rj_in = (Ain*Bin)/sqrt((Bin*cos(phi_ji))**2 + (Ain*sin(phi_ji))**2)
    Rj_out = (Aout*Bout)/sqrt((Bout*cos(phi_ji))**2 + (Aout*sin(phi_ji))**2)

    d_ij = (Ri_in + Rj_in) - L_ij
    
#    if (Ain >= Bin):
#        dmax_in = 2*Ain
#    else:
#        dmax_in = 2*Bin
#
    dmax_in = 2*Ain
    #dmax_in = Ri_in + Rj_in
    
    d_ij = d_ij/dmax_in
    
    if (d_ij > 0):
        cell_i.pos_overlap += d_ij
    else:
        cell_i.neg_overlap += d_ij
        
    
    dout = ((Ri_in + Rj_in) - (Ri_out + Rj_out))/dmax_in
    
    katt = i_poe_katt
    krep = i_poe_krep
    
#    if (d_ij >= 0.9*dmax_in):
#        krep *= 1000*krep
        
    if (d_ij >= 0.):
        f_ij = -krep*(d_ij**2)*r_vect_hat
    else:
        if (d_ij >= dout):
            f_ij = katt*(d_ij**2)*r_vect_hat
        else:
            f_ij = Vect(0., 0., 0.)
            
    return f_ij
    