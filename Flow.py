#------------------------------------------------------------------------------
# The Flow class contains the methods needed to solve for flow within a generalized
# vascular network consisting of Vessels
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

from Vessel import *
import numpy

#------------------------------------------------------------------------------
def solve_for_flow(geom):

    n = len(geom.nodes)                                                              # Number of nodes       
    v = len(geom.vessels)                                                            # Number of Vessels
    nBCs = len(geom.BCs)                                                             # Number of boundary conditions (known pressures)
    
    p = numpy.zeros(n)                                                          # Nodal pressure array
    A = numpy.zeros((n,n))                                                      # Conductance matrix
    b = numpy.zeros(n)                                                          # Solution array
    

    # Construct the conductance matrix and solution array
    for i in range(n):
          for vess in geom.vessels:
    
              if (vess.n1 == i):
                   if (vess.n0 in geom.BCs):
                       b[i] -= vess.G*p[vess.n0]
                       A[i][i] -= vess.G
                   else:
                       A[i][i] -= vess.G
                       A[i][vess.n0] += vess.G
                     
              if (vess.n0 == i):
                   if (vess.n1 in geom.BCs):
                       b[i] -= vess.G*p[vess.n1]
                       A[i][i] -= vess.G
                   else:
                       A[i][i] -= vess.G
                       A[i][vess.n1] += vess.G
    
    mean_G = 0.
    
    for vess in geom.vessels:
        mean_G += vess.G
        
    mean_G /= len(geom.vessels)
    
    if (i_inlet_flow_on == False):
        # Apply the pressure boundary conditions
        for BC in geom.BCs:
            for i in range(0,n):
                A[BC[0]][i] = 0.
            
            A[BC[0]][BC[0]] = mean_G
            b[BC[0]] = mean_G*BC[1]
            p[BC[0]] = BC[1]
    
    # Solve the system of equations for unknown pressures
    p = numpy.linalg.solve(A, b)
    #A_pinv = numpy.linalg.pinv(A)
    
    #p = numpy.dot(A_pinv, b)
    
    # Calculate flow within each Vessel
    for vess in geom.vessels:
        vess.P0 = p[vess.n0]
        vess.P1 = p[vess.n1]
        
        vess.calc_flow()
    
    
    # Determine if bifurcations are flow-converging or not
         
    # Find all vessels involved with each bifurcation
    for i in range(len(geom.bifur_nodes)):
        geom.bifur[i] = 0
        bif_node = geom.bifur_nodes[i]
        vess_n0 = []
        vess_n1 = []
        
        for vess in geom.vessels:
            if (vess.n0 == bif_node):
                vess_n0.append(vess.ID)
                
            if (vess.n1 == bif_node):
                vess_n1.append(vess.ID)
        
        # Determine if the bifurcation if flow-converging or flow-diverging
        if (len(vess_n0) == 2 and len(vess_n1) == 1):
            v01 = vess_n0[0]-1
            v02 = vess_n0[1]-1
            v11 = vess_n1[0]-1
            
            if (geom.vessels[v01].Q > 0 and geom.vessels[v02].Q > 0 and geom.vessels[v11].Q > 0):
                geom.bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
            if ((geom.vessels[v01].Q > 0 and geom.vessels[v02].Q < 0) or (geom.vessels[v01].Q < 0 and geom.vessels[v02].Q > 0)) and geom.vessels[v11].Q > 0:
                geom.bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))
            
            if (geom.vessels[v01].Q < 0 and geom.vessels[v02].Q < 0 and geom.vessels[v11].Q < 0):
                geom.bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))
            
            if ((geom.vessels[v01].Q > 0 and geom.vessels[v02].Q < 0) or (geom.vessels[v01].Q < 0 and geom.vessels[v02].Q > 0)) and geom.vessels[v11].Q < 0:
                geom.bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
        
        if (len(vess_n0) == 1 and len(vess_n1) == 2):
            v01 = vess_n0[0]-1
            v11 = vess_n1[0]-1
            v12 = vess_n1[1]-1
            
            if (geom.vessels[v01].Q > 0 and geom.vessels[v11].Q > 0 and geom.vessels[v12].Q > 0):
                geom.bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))
                
            if (geom.vessels[v01].Q < 0 and geom.vessels[v11].Q < 0 and geom.vessels[v12].Q < 0):
                geom.bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
            if geom.vessels[v01].Q > 0 and ((geom.vessels[v11].Q < 0 and geom.vessels[v12].Q > 0) or (geom.vessels[v11].Q > 0 and geom.vessels[v12].Q < 0)):
                geom.bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
            if geom.vessels[v01].Q < 0 and ((geom.vessels[v11].Q < 0 and geom.vessels[v12].Q > 0) or (geom.vessels[v11].Q > 0 and geom.vessels[v12].Q < 0)):
                geom.bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))