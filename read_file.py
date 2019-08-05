import numpy as np
import math
import matplotlib.pyplot as plt
import aem
#import main

def extract_file(filename):
    """ 
    Extract list of nodes and elements from file
    and returns Elements and Nodes object
    Eg file:
        24 35                                         # 24 nodes and 35 element definitions
        1 0.000000e+000 0.000000e+000 0.000000e+000
        2 2.200000e+000 0.000000e+000 0.000000e+000
        ...
        20 102 21 22                                  # boundary elements
        22 102 23 24 
        ...
        33 204 7 6 23 22                              # rectangular elements
        28 204 12 11 18 17 
        ...
    """
    # Open mesh file
    f = open(filename,'r')
    l = list(f)
    f.close()

    #    
    
    lsplit=[]
    
    # split and convert each row of the list to floats
    for i in l:
        each_row = i.split()
        for j in range(len(each_row)):
            each_row[j] = float(each_row[j])
        lsplit.append(each_row)
    
    # the first no is the number of nodes and second is no of elements
    node_num = int(lsplit[0][0])
    ele_num = int(lsplit[0][1])

    nodes = lsplit[1:node_num+1]
    nodes_dic = {}
    for i in range(len(nodes)):
        k = int(nodes[i][0])
        nodes_dic[k] = nodes[i][1:]
    
    remaining =lsplit[1+node_num:1+node_num+ele_num]
    
    # convert remaining elements into integer
    for i in range(len(remaining)):
        for j in range(len(remaining[i])):
            remaining[i][j] = int(remaining[i][j])
            
    # split remining as boundries and elements
    bc = []
    ele = []
    bc_dic = {}
    ele_dic ={}
    for i in remaining:
        if int(i[1]) == 102:
            bc.append(i)
            bc_dic[i[0]] = i[2:]
        else:
            ele.append(i)
            ele_dic[i[0]] = i[2:]
            
    bc = np.array(bc)
    ele = np.array(ele)
    return {"nodes":np.array(nodes),"bc":bc,"ele":ele,\
             "nodes_dic":nodes_dic, "bc_dic":bc_dic, "ele_dic":ele_dic}

# convert b

def distance(node_1,node_2):
    return math.sqrt((node_2[0]-node_1[0])**2 + (node_2[1]-node_1[1])**2)




def both_in(a,b,x):
    """true if both a and b are elements of x"""
    return a in x and b in x

#def body():

#############################################################

# Ectract from file
a = extract_file("prob2.dat")
nodes_array = a["nodes"]
ele_array = a["ele"]
bc = a["bc"]
nodes_dic = a["nodes_dic"]
ele_dic = a["ele_dic"]
bc_dic = a["bc_dic"]

# Create nodes object
nodes = aem.Nodes()
nodes.addNewNodes(nodes_array)

# Creating elements object
elements = aem.Elements()
elements.addNewElements(ele_array,nodes, 2.07e9, 79.3e9)

# Plot elements
fig1 = plt.figure(1)
fig1.clf()

for e in elements.list:
    node_1 = e.nodes[0].getCoord()
    node_2 = e.nodes[1].getCoord()
    node_3 = e.nodes[2].getCoord()
    node_4 = e.nodes[3].getCoord()
    cords = np.vstack([node_1,node_2,node_3,node_4])
    plt.plot(cords[:,0],cords[:,1])
    plt.text(e.x,e.y,e.ele_no)

plt.axis('equal')
plt.show()
plt.draw()


########
""" Spring no."""
no_spr = 20
ele_com = []
for i in range(len(elements.list)):
    for k in range(i+1,len(elements.list)):
        i_nodes = [elements.list[i].nodes[0].node_no,\
                    elements.list[i].nodes[1].node_no,\
                    elements.list[i].nodes[2].node_no,\
                    elements.list[i].nodes[3].node_no]
        k_nodes = [elements.list[k].nodes[0].node_no,\
                    elements.list[k].nodes[1].node_no,\
                    elements.list[k].nodes[2].node_no,\
                    elements.list[k].nodes[3].node_no]

        if both_in(i_nodes[0], i_nodes[1], k_nodes):
            ele_com.append([elements.list[i],elements.list[k],1,3,no_spr])
        elif both_in(i_nodes[1], i_nodes[2], k_nodes):
            ele_com.append([elements.list[i],elements.list[k],4,2,no_spr])
        elif both_in(i_nodes[2], i_nodes[3], k_nodes):
            ele_com.append([elements.list[i],elements.list[k],3,1,no_spr])
        elif both_in(i_nodes[3], i_nodes[0], k_nodes):
            ele_com.append([elements.list[i],elements.list[k],2,4,no_spr])
            
ele_com_mat = ele_com


ele_com_mat = np.array(ele_com_mat)

# create Springs object and then add new springs to it using ele_com_mat
springs = aem.Springs()
springs.addNewSprings_auto(ele_com_mat)

glob_mat = elements.global_stiffness_matrix(springs)

n = len(elements.list)

all_dis = range((3*n))



# Boundary Condision
fixed_ele = [int(x) for x in raw_input("Enter element numbers to be fixed: ").split()]
fixed_dis = []
for ele_no in fixed_ele:
    x = elements.ele2id[ele_no]
    fixed_dis.append(3*x)
    fixed_dis.append(3*x+1)
    fixed_dis.append(3*x+2)

unknown_dis_list = []
for x in all_dis:
    if not(x in fixed_dis):
        unknown_dis_list.append(x)
our_map = {}
for i in range(len(unknown_dis_list)):
    our_map[unknown_dis_list[i]] =  i


loads = np.zeros((len(unknown_dis_list),1))
load_ele = [int(x) for x in raw_input("Enter element no. on with to provide vertical load: ").split()]
for ele_no in load_ele:
    i = elements.ele2id[ele_no]
    loads[(our_map[i*3+1]),0] = -1000

unknown_dis = np.array([unknown_dis_list])
K=glob_mat[unknown_dis.T,unknown_dis]

dis = np.linalg.solve(K,loads)

for ele_no in load_ele:
    i = elements.ele2id[ele_no]
    print dis[(our_map[i*3+1]),0]

del_y = [dis[i,0] for i in range(len(dis)) if (i-1)%3==0]
ele_nos = [elements.list[i].ele_no for i in range(1,len(elements.list))]

fig2=plt.figure(2)
fig2.clf()
plt.plot(ele_nos,del_y)



#eprint (dis[(our_map[ elements.ele2id[46]*3+1]),0]+dis[(our_map[ elements.ele2id[45]*3+1]),0])/2

#if __name__ == '__main__':
#    body()