#call important libraries
import numpy as np
import igraph as ig
import scipy                              
import sys
from scipy.integrate import simpson




########################### new function #########################################
def cell_lineage_model(num_gene, noise_amp, network_mat, branch, g_root, tau, pow ):
    '''
    Description
    ----
    It gives a cell lineage graph contains gene expression of all genes throughout the lineage for 
    each cell. 

    Parameters
    ----
    num_gene : int, number of genes.
    noise_amp : float, noise amplitude.
    network_mat: In case independent gene expression dynamics, it's an 1D numpy array with shape 
    (num_gene, ).  
    branch : int, number of final branch. During taking a snapshot along a branch, it is better not to 
    include final branch.
    g_root : Initial node of ipython graph. It contains steady state gene expression.
    tau : float, cell division time.
    pow : int, time step size for cell division procedure is 10**-pow. 
    
    Returns
    ----
    g : igraph-python graph, with 2**branch nodes. Every node has three attributes: 
        name: list, [node number of parent cell, own node number], 
        gene: list, contains each time update of each gene
        time : list, contains each time step
    '''
    max_cell = int(2**branch) #maximum size of the graph

    # initialization of root cell
    g = ig.Graph()
    g.add_vertices(int(max_cell + 1))

    g.vs[0]['name'] = [0, 0]
    g.vs[0]['gene'] = []
    g.vs[0]['gene'].append(g_root)
    g.vs[0]['time'] = []
    g.vs[0]['time'].append(0)

    # parameters for cell lineage 
    dt = np.float32(10**-pow)
    steps_1 = int(np.round(tau / dt))


    # various initialization for cell lineage
    update_list = [0]  #Keeping track of cells
    counter = 0
    t = 0
    node = 0
    #to ensure positive gene expression state. Might have to change ...
    #... depending on noise amplitude.
    con_arr = np.ones(num_gene, dtype = np.float32)
    

    #cell lineage dynamics
    while counter < max_cell: #keep track of number of total cells
        remove_list = []     #removing parent cells
        adding_list = []     #adding daughter cells
        for i in range (0, len(update_list), 1):
            node = update_list[i]
            
            arr = np.copy(g.vs[node]['gene'][-1])
            #When we work with independent gene expression dynamics, network_mat...
            #... is an array of size (num_gene, ).
            du = (network_mat * arr) * dt + con_arr * dt + np.random.normal(0, noise_amp *dt**0.5, num_gene)
            # du = np.matmul(network_mat, arr) * dt + con_arr * dt + np.random.normal(0, noise_amp *dt**0.5, num_gene)

            # update for next time point
            arr += du
            
            # storing data
            g.vs[node]['gene'].append(np.copy(arr))
            g.vs[node]['time'].append(t)

            # cell division procedure
            if len(g.vs[node]['time']) > steps_1:
                new_node_1 = counter + 1
                new_node_2 = counter + 2

                if new_node_1 > int(max_cell) or new_node_2 > int(max_cell):
                    break
                
                remove_list.append(node)
                adding_list.append(new_node_1), adding_list.append(new_node_2)
                g.add_edges([(node,new_node_1)]), g.add_edges([(node,new_node_2)])

                #new_node_1
                g.vs[new_node_1]['name'] = [node, counter + 1]
                g.vs[new_node_1]['gene'] = []
                g.vs[new_node_1]['gene'].append(np.copy(g.vs[node]['gene'][-1]))
                g.vs[new_node_1]['time'] = [0]

                #new_node_2
                g.vs[new_node_2]['name'] = [node,counter + 2]
                g.vs[new_node_2]['gene'] = []
                g.vs[new_node_2]['gene'].append(np.copy(g.vs[node]['gene'][-1]))
                g.vs[new_node_2]['time'] = [0]

                #update counter
                counter = counter + 2
        
        
        update_list = update_list + adding_list
        update_list = [x for x in update_list if x not in remove_list]
        t = np.round(t + dt, pow + 1)   

    return g