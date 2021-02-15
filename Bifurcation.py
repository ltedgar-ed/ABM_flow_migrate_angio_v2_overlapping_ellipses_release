#------------------------------------------------------------------------------
# The Bifurcation class contains the bifurcation rules for cell migration

# Lowell Taylor Edgar
# University of Edinburgh
# 2020

from LTETools import *
from Input import *

#------------------------------------------------------------------------------
# Handle a bifurcation point while migrating
def handle_bifurcation(migrating_cell, parent_vess, bif_node, branch1, branch2):
    # Determine if flow within branches are incoming or outgoing of the bifurcation
    if (branch1.n1 == bif_node):
        if (branch1.Q > 0.):
            branch1_status = "incoming"
        else:
            branch1_status = "outgoing"
            
    if (branch1.n0 == bif_node):
        if (branch1.Q > 0.):
            branch1_status = "outgoing"
        else:
            branch1_status = "incoming" 
            
    if (branch2.n1 == bif_node):
        if (branch2.Q > 0.):
            branch2_status = "incoming"
        else:
            branch2_status = "outgoing"
            
    if (branch2.n0 == bif_node):
        if (branch2.Q > 0.):
            branch2_status = "outgoing"
        else:
            branch2_status = "incoming" 
    
#    if (abs(branch1.Q) < 10):
#        branch1_status = "zero"
#    
#    if (abs(branch2.Q) < 10):
#        branch2_status = "zero"
        
            
    # If flow is only incoming in one branch, always choose that branch (ie, against flow)
    if (branch1_status == "incoming" and branch2_status == "outgoing"):
        migrating_cell.vessID = branch1.ID
        branch1.cells.append(migrating_cell)
        return
    
    if (branch2_status == "incoming" and branch1_status == "outgoing"):
        migrating_cell.vessID = branch2.ID
        branch2.cells.append(migrating_cell)
        return
    
#    # If flow is zero in one branch and not the other, choose the branch with flow
#    if (branch1_status != "zero" and branch2_status == "zero"):
#        migrating_cell.vessID = branch1.ID
#        branch1.cells.append(migrating_cell)
#        return
#    
#    if (branch1_status == "zero" and branch2_status != "zero"):
#        migrating_cell.vessID = branch2.ID
#        branch2.cells.append(migrating_cell)
#        return
    
    # If a flow-converging bifurcation, enact bifurcation rule
    if (branch1_status == "incoming" and branch2_status == "incoming"):        
        bifurcation_rule(migrating_cell, parent_vess, bif_node, branch1, branch2)
        return

    # If a flow-diverging bifurcation, enact bifurcation rule
    if (branch1_status == "outgoing" and branch2_status == "outgoing"):        
        bifurcation_rule(migrating_cell, parent_vess, bif_node, branch1, branch2)
        return
    
#        # If a flow-converging bifurcation, enact bifurcation rule
#    if (branch1_status == "zero" and branch2_status == "zero"):        
#        bifurcation_rule(migrating_cell, parent_vess, bif_node, branch1, branch2)
#        return
    
    print('handle_bifurcation made it this far without returning!')
    migrating_cell.migrate = 0
    return
    

#------------------------------------------------------------------------------
# Determine cell behavior using the specified bifurcation rule
def bifurcation_rule(migrating_cell, parent_vess, bif_node, branch1, branch2):
    # Find the angle between parent vessel and both branches
    angle1 = findangle2D(parent_vess.unit, branch1.unit)
    angle2 = findangle2D(parent_vess.unit, branch2.unit)
        
    
    # Set the left branch to be the branch with the greatest angle
    if (angle1 > angle2):
        left_branch = branch1
        right_branch = branch2
    else:
        left_branch = branch2
        right_branch = branch1
       
        
    # Bifurcation Rule 1 - Simple geometric 
    if (i_bifurcation_rule == 1):
        if (parent_vess.n1 == bif_node):
            if (migrating_cell.zeta <= 0.5):
                migrating_cell.vessID = left_branch.ID
                left_branch.cells.append(migrating_cell)
                return
            
            if (migrating_cell.zeta > 0.5):
                migrating_cell.vessID = right_branch.ID
                right_branch.cells.append(migrating_cell)
                return
        else:
            if (migrating_cell.zeta <= 0.5):
                migrating_cell.vessID = right_branch.ID
                right_branch.cells.append(migrating_cell)
                return
            
            if (migrating_cell.zeta > 0.5):
                migrating_cell.vessID = left_branch.ID
                left_branch.cells.append(migrating_cell)
                return
     
        
    # Bifurcation Rule 2 - Geometric with diameter control
    if (i_bifurcation_rule == 2):
        if ((left_branch.num_cells + right_branch.num_cells) == 0.):
            print('bifurcation_rule zero cells in both branches')
            migrating_cell.migrate = 0
            return
           
        zeta_left = left_branch.num_cells/(left_branch.num_cells + right_branch.num_cells)
        zeta_right = right_branch.num_cells/(left_branch.num_cells + right_branch.num_cells)
        
        if (parent_vess.n1 == bif_node):        
            left_pt = migrating_cell.zeta
            right_pt = migrating_cell.zeta
    
            if ((left_pt >= 0.25) and (left_pt <= 0.25 + zeta_left/2.)) or ((left_pt - 1. < 0.25) and (left_pt - 1. > 0.25 - zeta_left/2.)):
                migrating_cell.vessID = left_branch.ID
                left_branch.cells.append(migrating_cell)
                return
            
            if ((left_pt >= 0.25) and (left_pt <= 0.25 + zeta_left/2.)) or ((left_pt < 0.25) and (left_pt > 0.25 - zeta_left/2.)):
                migrating_cell.vessID = left_branch.ID
                left_branch.cells.append(migrating_cell)
                return
            
            if ((right_pt >= 0.75) and (right_pt < 0.75 + zeta_right/2.)) or ((right_pt < 0.75) and (right_pt >= 0.75 - zeta_right/2.)):
                migrating_cell.vessID = right_branch.ID
                right_branch.cells.append(migrating_cell)
                return 
            
            if ((right_pt + 1. >= 0.75) and (right_pt + 1. < 0.75 + zeta_right/2.)) or ((right_pt < 0.75) and (right_pt >= 0.75 - zeta_right/2.)):
                migrating_cell.vessID = right_branch.ID
                right_branch.cells.append(migrating_cell)
                return
        else:
            right_pt = migrating_cell.zeta
            left_pt = migrating_cell.zeta
    
            if ((right_pt >= 0.25) and (right_pt <= 0.25 + zeta_right/2.)) or ((right_pt < 0.25) and (right_pt > 0.25 - zeta_right/2.)):
                migrating_cell.vessID = right_branch.ID
                right_branch.cells.append(migrating_cell)
                return
            
            if ((right_pt >= 0.25) and (right_pt <= 0.25 + zeta_right/2.)) or ((right_pt - 1. < 0.25) and (right_pt - 1. > 0.25 - zeta_right/2.)):
                migrating_cell.vessID = right_branch.ID
                right_branch.cells.append(migrating_cell)
                return
            
            if ((left_pt + 1. >= 0.75) and (left_pt + 1. < 0.75 + zeta_left/2.)) or ((left_pt < 0.75) and (left_pt >= 0.75 - zeta_left/2.)):
                migrating_cell.vessID = left_branch.ID
                left_branch.cells.append(migrating_cell)
                return
    
            if ((left_pt >= 0.75) and (left_pt < 0.75 + zeta_left/2.)) or ((left_pt < 0.75) and (left_pt >= 0.75 - zeta_left/2.)):
                migrating_cell.vessID = left_branch.ID
                left_branch.cells.append(migrating_cell)
                return
            
    # If we somehow didn't trigger a return
    print('bifrucation_rule made it this far without returning!')
    return
    

#------------------------------------------------------------------------------
# The simplest bifurcation rule, turn left or right based just on what side you're on
def enhanced_geometric(migrating_cell, parent_vess, bif_node, branch1, branch2):
    
    # Find the angle between parent vessel and both branches
    angle1 = findangle2D(parent_vess.unit, branch1.unit)
    angle2 = findangle2D(parent_vess.unit, branch2.unit)
        
    # Set the left branch to be the branch with the greatest angle
    if (angle1 > angle2):
        left_branch = branch1
        right_branch = branch2
    else:
        left_branch = branch2
        right_branch = branch1
    
    zeta_left = left_branch.num_cells/(left_branch.num_cells + right_branch.num_cells)
    zeta_right = right_branch.num_cells/(left_branch.num_cells + right_branch.num_cells)
        
    # If exiting parent vessel from +1 terminus
    if (parent_vess.n1 == bif_node):        
        cell_pt = migrating_cell.zeta

        if ((cell_pt >= 0.25) and (cell_pt <= 0.25 + zeta_left/2.)) or ((cell_pt - 1. < 0.25) and (cell_pt - 1. > 0.25 - zeta_left/2.)):
            assign_branch(left_branch, migrating_cell)
            return
        
        if ((cell_pt >= 0.25) and (cell_pt <= 0.25 + zeta_left/2.)) or ((cell_pt < 0.25) and (cell_pt > 0.25 - zeta_left/2.)):
            assign_branch(left_branch, migrating_cell)
            return
        
        if ((cell_pt >= 0.75) and (cell_pt < 0.75 + zeta_right/2.)) or ((cell_pt < 0.75) and (cell_pt >= 0.75 - zeta_right/2.)):
            assign_branch(right_branch, migrating_cell)
            return 
        
        if ((cell_pt + 1. >= 0.75) and (cell_pt + 1. < 0.75 + zeta_right/2.)) or ((cell_pt < 0.75) and (cell_pt >= 0.75 - zeta_right/2.)):
            assign_branch(right_branch, migrating_cell)
            return
    # If exiting from 0 terminus
    else:
        cell_pt = migrating_cell.zeta

        if ((cell_pt >= 0.25) and (cell_pt <= 0.25 + zeta_right/2.)) or ((cell_pt < 0.25) and (cell_pt > 0.25 - zeta_right/2.)):
            assign_branch(right_branch, migrating_cell)
            return
        
        if ((cell_pt >= 0.25) and (cell_pt <= 0.25 + zeta_right/2.)) or ((cell_pt - 1. < 0.25) and (cell_pt - 1. > 0.25 - zeta_right/2.)):
            assign_branch(right_branch, migrating_cell)
            return
        
        if ((cell_pt + 1. >= 0.75) and (cell_pt + 1. < 0.75 + zeta_left/2.)) or ((cell_pt < 0.75) and (cell_pt >= 0.75 - zeta_left/2.)):
            assign_branch(left_branch, migrating_cell)
            return

        if ((cell_pt >= 0.75) and (cell_pt < 0.75 + zeta_left/2.)) or ((cell_pt < 0.75) and (cell_pt >= 0.75 - zeta_left/2.)):
            assign_branch(left_branch, migrating_cell)
            return
     

#------------------------------------------------------------------------------
# The simplest bifurcation rule, turn left or right based just on what side you're on
def basic_geometric(migrating_cell, parent_vess, bif_node, branch1, branch2):
    
    # Find the angle between parent vessel and both branches
    angle1 = findangle2D(parent_vess.unit, branch1.unit)
    angle2 = findangle2D(parent_vess.unit, branch2.unit)
        
    # Set the left branch to be the branch with the greatest angle
    if (angle1 > angle2):
        left_branch = branch1
        right_branch = branch2
    else:
        left_branch = branch2
        right_branch = branch1
       
    # If exiting parent vessel from +1 terminus
    if (parent_vess.n1 == bif_node):
        if (migrating_cell.zeta <= 0.5):
            assign_branch(left_branch, migrating_cell)
            return
        
        if (migrating_cell.zeta > 0.5):
            assign_branch(right_branch, migrating_cell)
            return
    # If exiting from 0 terminus
    else:
        if (migrating_cell.zeta <= 0.5):
            assign_branch(right_branch, migrating_cell)
            return
        
        if (migrating_cell.zeta > 0.5):
            assign_branch(left_branch, migrating_cell)
            return


#------------------------------------------------------------------------------
# Assign a cell to a chosen target branch
def assign_branch(target_branch, migrating_cell):
    
    migrating_cell.vessID = target_branch.ID
    target_branch.cells.append(migrating_cell)
