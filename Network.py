#------------------------------------------------------------------------------
# The Network class contains the methods that require knowledge of the vascular network as a whole
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

from math import pi
from math import acos
from math import cos
from math import sin
from math import sqrt
from numpy import sign

from Vessel import *
from LTETools import *
from Input import *
from ECell import *

import POE
import Bifurcation

#------------------------------------------------------------------------------    
# Realign polarity vectors in cells based on weights
def realign_polarity(geom):
    # Loop through all cells in all vessels
    for vess in geom.vessels:
        for cell in vess.cells:
            
            # Original polarity vector
            pol_old = cell.polarity.copy()
            
            # Phi1 - Realignment angle due to persistence
            phi1 = 0.
            
            # Phi2 - Realignment angle due to flow
            # Flow vector
            flow_vect = -float(sign(vess.Q))*vess.unit.copy()
            
            if (cell.mig_direction == 'with'):
                flow_vect = float(sign(vess.Q))*vess.unit.copy()
            
            #flow_vect = -float(sign(vess.Q))*Vect(1., 0., 0.)
            
#            if (abs(vess.Q) < 1e-10):
#                flow_vect = pol_old.copy()
                
            phi2 = findangle2D(pol_old, flow_vect)          
            
#            if (i_polarity_shift == True):
#                phi2 += cell.pol_shift
            
            # Phi3 - Realignment angle due to random walk component
            # Random walk component
            rand_vect = Vect(rand.uniform(-1, 1), rand.uniform(-1, 1), 0.)
            rand_vect.unit()
            phi3 = findangle2D(pol_old, rand_vect)
            
            # Calculate new realignment angle from weighted average
            theta = i_w1*phi1 + i_w2*phi2 + i_w3*phi3
            
            # Calculate the new polarity vector by rotating the old vector by theta
            Q = Tensor2O([cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 0])
            pol_new = Q*pol_old
            pol_new.unit()
            cell.polarity = pol_new.copy()
            
            assert cell.polarity.length() != 0., "polarity vector length dropped to zero"


def cell_randomise_force(geom):
    
    for vess in geom.vessels:
        for cell in vess.cells:
            cell.net_force.x = rand.uniform(-i_vess_length/2, i_vess_length/2)
            cell.net_force.y = rand.uniform(-i_vess_length/2, i_vess_length/2)
            

#------------------------------------------------------------------------------            
# Calculate force balance on cells
def cell_force_balance(geom):
    ## Zero out the net force for all cells to begin
    for vess in geom.vessels:
        for cell in vess.cells:
            cell.net_force.zero()
            cell.pos_overlap = 0.
            cell.neg_overlap = 0.
            
    
#    ## Calculate the FSI force (r component of net force)
#    for vess in geom.vessels:
#        Z = vess.L
#        R = vess.D/2
#        
#        Pint = (vess.P0 + vess.P1)/2
#        
#        f_fsi = (Pint - i_Pext)*(2*pi*R)*Z
#        
#        for cell in vess.cells:
#            cell.net_force.z = f_fsi
            
    ## Calculate the overlap forces (z and c components of net force)
    ## Add contribution to net force from other cells in the same vessel        
    for vess in geom.vessels:
        for cell_i in vess.cells:
            for cell_j in vess.cells:
                Z = vess.L
                R = vess.D/2

                # If two cells are on top of each other, nudge cell j to create an overlap reaction force
                if (cell_i.xi == cell_j.xi) and (cell_i.zeta == cell_j.zeta) and (cell_i.ID != cell_j.ID):
                    cell_j.zeta *= 0.90
                
                if (cell_i.ID != cell_j.ID):
                    r_vect = POE.distance_between_cells_same_vess(cell_i, cell_j, Z, R)
                    f_ij = POE.calc_overlap_force(cell_i, cell_j, r_vect)
                    cell_i.net_force += f_ij
                    
#                    # Add circumfrential contribution to radial net force
#                    # If an attractive force...
#                    if (r_vect.y >= 0.) and (f_ij.y >= 0.):
#                        cell_i.net_force.z -= f_ij.y
#                        
#                    if (r_vect.y < 0.) and (f_ij.y < 0.):
#                        cell_i.net_force.z -= -f_ij.y
#                    
#                    # If a repulsive force...
#                    if (r_vect.y >= 0.) and (f_ij.y < 0.):
#                        cell_i.net_force.z -= f_ij.y
#                        
#                    if (r_vect.y < 0.) and (f_ij.y > 0.):
#                        cell_i.net_force.z -= -f_ij.y
                
                
    ## Add contribution to net force from cells in upstream neighbour(s)
    for vess1 in geom.vessels:
        for k in range(len(vess1.neigh0)):
            vess2 = geom.vessels[vess1.neigh0[k]-1]
            
            for cell_i in vess1.cells:
                for cell_j in vess2.cells:
                    Z1 = vess1.L
                    R1 = vess1.D/2
                    Z2 = vess2.L
                    R2 = vess2.D/2
                    
                    r_vect = POE.distance_between_cell_and_upstream_neigh(cell_i, Z1, R1, cell_j, Z2, R2)
                    f_ij = POE.calc_overlap_force(cell_i, cell_j, r_vect)
                    cell_i.net_force += f_ij
    
    
    ## Add contribution to net force from cells in downstream neighbour(s)
    for vess1 in geom.vessels:
        for k in range(len(vess1.neigh1)):
            vess2 = geom.vessels[vess1.neigh1[k]-1]
            
            for cell_i in vess1.cells:
                for cell_j in vess2.cells:
                    Z1 = vess1.L
                    R1 = vess1.D/2
                    Z2 = vess2.L
                    R2 = vess2.D/2
                    
                    r_vect = POE.distance_between_cell_and_downstream_neigh(cell_i, Z1, R1, cell_j, Z2, R2)
                    f_ij = POE.calc_overlap_force(cell_i, cell_j, r_vect)
                    cell_i.net_force += f_ij


#    for vess in geom.vessels:
#        R = vess.D/2
#        
#        for cell in vess.cells:
#            cell.net_force.y *= 0.9
            
    
#------------------------------------------------------------------------------            
# Perform cell migration
def cell_migration(geom, mig_on, fsi_on, iters):
    ## Set the migration force magnitude
    kmig = i_poe_kmig
    dt = i_poe_dt
    
    if (mig_on == False):
        kmig = 0.
        
        
#    if (i_fsi_on == False):
#        for vess in geom.vessels:
#            for cell in vess.cells:
#                cell.net_force.z = 0.
#    
#    if (i_fsi_on == True) and (fsi_on == False):
#        for vess in geom.vessels:
#            for cell in vess.cells:
#                cell.vel.net_force = 0.
                
                
    ## Calculate the migration velocity for each cell from force balance
    for vess in geom.vessels:
        for cell in vess.cells:
            # Determine the components of the polarity vector in z
            cell_pol_z = cell.polarity*vess.unit
            
            # Determine the components of the polarity vector in c 
            angle = findangle2D(Vect(1.,0.,0.), vess.unit)
            ROT = Tensor2O([cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 0])
                        
            theta = 2*pi*cell.zeta
            
            # Circumfrence unit vector
            e_c = Vect(0., cos(theta), -sin(theta))
            
            cell_pol_c = cell.polarity*(ROT*e_c)

            # Polarity vector in z and c
            cell_pol = Vect(cell_pol_z, cell_pol_c, 0.)
            cell_pol.unit()
            
            # Calculate cell migration velocity through force balance
            cell.old_vel = cell.vel.copy()
            
            cell.vel = (cell.net_force + kmig*cell_pol)/i_poe_eta
            
            v_c = cell.vel.y
            R = vess.D/2.
            v_theta = v_c/R
            cell.vel.y = v_theta
            
            cell.acc = (cell.vel - cell.old_vel)/dt
            
    
    ## Zero out migration fields
    for vess in geom.vessels:
        for cell in vess.cells:
            cell.migrate = 0
            cell.rem_dist.zero()
           
            
    ## Displace cells and determine they are migrating into a neighbouring segment
    for vess in geom.vessels:
        #vess.cell_move = False
        
        for cell in vess.cells:
            # Calculate cell displacement
            u = cell.vel*dt
            u_z = u.x
            u_theta = u.y
            u_r = u.z
            
            cell.disp = Vect(u_z, u_theta, u_r)
               
            Z = vess.L
            R = vess.D/2
            
            #u_theta = u_c/R
            
            # Update z, c, and theta position
            old_z = Z*cell.xi
            old_c = 2*pi*R*cell.zeta
            old_theta = 2*pi*cell.zeta
            
            new_z = old_z + u_z
            new_theta = old_theta + u_theta
            
            # Apply periodicity to theta
            while (new_theta > 2*pi):
                new_theta = new_theta - 2*pi
                
            while (new_theta < 0.):
                new_theta = new_theta + 2*pi
            
            new_r = R + u_r
            
            D0 = vess.num_cells*i_cell_size/pi
            R0 = D0/2
            
            # If migrating out the +1 terminus...
            if (new_z > Z - 1e-5):
                if (u_z != 0.):
                    alpha = (Z - old_z)/u_z
                else:
                    alpha = 1.
                
                cell.migrate = 1
                new_z = Z
                new_theta = old_theta + alpha*u_theta
                
                cell.rem_dist.x = (1 - alpha)*u_z
                cell.rem_dist.y = (1 - alpha)*u_theta
                
                #vess.cell_move = True
            
            # If migrating out the 0 terminus...
            if (new_z < 1e-5):
                if (u_z != 0.):
                    alpha = -old_z/u_z
                else:
                    alpha = 1.
                    
                cell.migrate = -1
                new_z = 0.
                new_theta = old_theta + alpha*u_theta
                
                cell.rem_dist.x = (1 - alpha)*u_z
                cell.rem_dist.y = (1 - alpha)*u_theta
                
                #vess.cell_move = True
            
            # Update the cell
            cell.xi = new_z/Z
            cell.zeta = new_theta/(2*pi)
            
            if (i_fsi_linear_PD == False):
                cell.rho = new_r/R0
            else:
                cell.rho = new_r/R
                
        
    ## Move the cells that are migrating to their new vessel        
    for vess in geom.vessels:
        cells_that_migrated = []
        
        for i in range(len(vess.cells)):
            
            # If cell is migrating out the 0 terminus...
            if (vess.cells[i].migrate == -1):
                migrating_cell = vess.cells[i]
                
                # If only one neighbour...
                if (len(vess.neigh0) == 1):
                    migrating_cell.vessID = vess.neigh0[0]
                    geom.vessels[vess.neigh0[0]-1].cells.append(migrating_cell)
                
                # If more than one neighbour...
                elif (len(vess.neigh0) == 2):
                    Bifurcation.handle_bifurcation(migrating_cell, vess, vess.n0, geom.vessels[vess.neigh0[0]-1], geom.vessels[vess.neigh0[1]-1])
                
                # Migrate the cell
                if (migrating_cell.migrate != 0):
                    cells_that_migrated.append(i)
                    migrating_cell.migrate = 0
                    vess2 = geom.vessels[migrating_cell.vessID-1]
                    Z2 = vess2.L
                    R2 = vess2.D/2
                    
                    # Define the rotation transform between the two vessels
                    angle = findangle2D(vess.unit, vess2.unit)
                    ROT = Tensor2O([cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 0])
            
#                    if (migrating_cell.rem_dist.x > 0):
#                        assert migrating_cell.rem_dist.x <= 0, "rem_dist wrong sign"
#                    
                    # Update z position assuming cell is entering neighbour at +1 terminus
                    new_z = Z2 + migrating_cell.rem_dist.x
                    
                    # Correct if we were wrong in this assumption
                    for g in range(len(vess2.neigh0)):
                        if (vess2.neigh0[g] == vess.ID):
                            new_z = -migrating_cell.rem_dist.x                        
                    
                    # Update xi position
                    migrating_cell.xi = new_z/Z2
                    
                    old_theta = 2*pi*migrating_cell.zeta

                    new_theta = old_theta + migrating_cell.rem_dist.y

                    while (new_theta > 2*pi):
                        new_theta = new_theta - 2*pi
                        
                    while (new_theta < 0.):
                        new_theta = new_theta + 2*pi
                   
                    # Finish updating cell
                    migrating_cell.zeta = new_theta/(2*pi)
                    migrating_cell.polarity = ROT*migrating_cell.polarity
            
            
            # If cell is migrating out the +1 terminus...
            if (vess.cells[i].migrate == 1):
                migrating_cell = vess.cells[i]
                
                # If only one neighbour...
                if (len(vess.neigh1) == 1):
                    migrating_cell.vessID = vess.neigh1[0]
                    geom.vessels[vess.neigh1[0]-1].cells.append(migrating_cell)
                               
                # If more than one neighbour...
                elif (len(vess.neigh1) == 2):
                    Bifurcation.handle_bifurcation(migrating_cell, vess, vess.n1, geom.vessels[vess.neigh1[0]-1], geom.vessels[vess.neigh1[1]-1])

                # Migrate the cell 
                if (migrating_cell.migrate != 0):
                    cells_that_migrated.append(i)
                    migrating_cell.migrate = 0
                    vess2 = geom.vessels[migrating_cell.vessID-1]
                    Z2 = vess2.L
                    R2 = vess2.D/2
                    
                    # Define the rotation transform between the two vessels
                    angle = findangle2D(vess.unit, vess2.unit)
                    ROT = Tensor2O([cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 0])
        
#                    if (migrating_cell.rem_dist.x < 0.):
#                        assert migrating_cell.rem_dist.x >= 0, "rem_dist wrong sign"
#                    
                    # Update z position assuming cell is entering neighbour at 0 terminus
                    new_z = migrating_cell.rem_dist.x
                    
                    # Correct if we were wrong in this assumption
                    for g in range(len(vess2.neigh1)):
                        if (vess2.neigh1[g] == vess.ID):
                            new_z = Z2 - migrating_cell.rem_dist.x  
                    
                    # Update xi position
                    migrating_cell.xi = new_z/Z2
                    
                    old_theta = 2*pi*migrating_cell.zeta

                    new_theta = old_theta + migrating_cell.rem_dist.y

                    while (new_theta > 2*pi):
                        new_theta = new_theta - 2*pi
                        
                    while (new_theta < 0.):
                        new_theta = new_theta + 2*pi
                   
                    # Finish updating cell
                    migrating_cell.zeta = new_theta/(2*pi)
                    migrating_cell.polarity = ROT*migrating_cell.polarity
                    
        # Remove the cells that left the vessel (or rather, keep only the cells that didn't leave)
        new_cells = []
        
        for i in range(len(vess.cells)):
            if i not in cells_that_migrated:
                new_cells.append(vess.cells[i])
                
        vess.cells = new_cells.copy()


    ## Update cell number in vessels after migration
    for vess in geom.vessels:
        vess.num_cells = len(vess.cells)
        
    ## Apply Dirichlet boundary conditions
    max_cell_num = 0
    
    # Determine the max cell ID number for adding new cells
    for vess in geom.vessels:
        for cell in vess.cells:
            if (cell.ID > max_cell_num):
                max_cell_num = cell.ID
    
    # Find vessels that have Dirichlet conditions and apply
    for vess in geom.vessels:
        if (vess.dirichlet != 0 and vess.num_cells != vess.dirichlet):
            # If cell number is less than dirichlet condition, add cells
            if (vess.num_cells < vess.dirichlet):
                while (vess.num_cells < vess.dirichlet):
                    max_cell_num += 1
                    xi_init = 0.5
                    zeta_init = 1e-3
                    vess.cells.append(ECell(max_cell_num, vess.ID, xi_init, zeta_init, vess.unit.copy()))
                    vess.num_cells = len(vess.cells)

            # If cell number is more than dirichlet condition, remove cells randomly
            elif (vess.num_cells > vess.dirichlet):
                while (vess.num_cells > vess.dirichlet):
                    vess.cells.pop(rand.randrange(len(vess.cells)))
                    vess.num_cells = len(vess.cells)
  
    
    ## Update vessels after migration
    for vess in geom.vessels:
        vess.old_num_cells = vess.num_cells
        vess.num_cells = len(vess.cells)
        vess.update_diameter()
        vess.update_conductance()

    geom.count_cells()
    
    
    ## If its a migration step, record the migration velocity and net force
    if (mig_on == True):
        for vess in geom.vessels:
            for cell in vess.cells:
                cell.mig_vel = cell.vel.copy()
                cell.mig_net_force = cell.net_force.copy()


#------------------------------------------------------------------------------
# Intercalation
def intercalate(geom, init_flag):
    z_c_max_iter = 1000
    acc_tol = i_del_vel_min/i_poe_dt
    
#    if (init_flag == True):
#        #z_c_tol /= 10
        

    if (i_find_acc_free == True):
        vel_norm = geom.calc_vel_norm()    
        acc_norm = geom.calc_acc_norm()
        max_vel = geom.calc_max_vel()
        max_acc = geom.calc_max_acc()
        z_c_iters = 0
    
        converged = False
        
        while (converged == False):
            z_c_iters += 1
            cell_force_balance(geom)
            cell_migration(geom, False, False, z_c_iters)
            vel_norm = geom.calc_vel_norm()
            acc_norm = geom.calc_acc_norm()
            max_vel = geom.calc_max_vel()
            max_acc = geom.calc_max_acc()
            #print(str(z_c_iters) + " " + str(acc_norm))
            

            if (acc_norm < acc_tol):
                converged = True

            if  (z_c_iters == z_c_max_iter):
                converged = True
                
            
        print("          z,c: {0} ({1}, {2})".format(z_c_iters, acc_norm, max_acc))
        geom.total_interc_iters += z_c_iters






#    fsi_max_iter = 10000
#    fsi_tol = 1e-4
#    
#    if (i_find_acc_free == True):
#        vel_norm = geom.calc_vel_norm()    
#        acc_norm = geom.calc_acc_norm()
#        max_acc = geom.calc_max_acc()
#        z_c_iters = 0
#        
#        while (abs(acc_norm) > z_c_tol) and (z_c_iters < z_c_max_iter):
#            z_c_iters += 1
#            cell_force_balance(geom)
#            cell_migration(geom, False, False, z_c_iters)
#            vel_norm = geom.calc_vel_norm()
#            acc_norm = geom.calc_acc_norm()
#            max_acc = geom.calc_max_acc()
#            
#            #print(vel_norm)
#        
#        print("          z,c: {0} ({1}, {2})".format(z_c_iters, acc_norm, max_acc))
#        
#        diam_norm = 100*fsi_tol
#        fsi_iters = 0
#        
#        if (i_fsi_on == True):
#            while (diam_norm > fsi_tol) and (fsi_iters < fsi_max_iter):
#                cell_force_balance(geom)
#                cell_migration(geom, False, True) # Need to update flow with diameter updates?
#                diam_norm = geom.calc_diam_norm()
#                fsi_iters += 1
#                print(diam_norm)
#        
#            print("          fsi: {0} ({1})".format(fsi_iters, diam_norm))