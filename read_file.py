import numpy
import math
import matplotlib.pyplot as plt
import aem
import main


f = open("prob3.dat",'r')
l = list(f)
f.close()
lsplit=[]
for i in l:
    lsplit.append(i.split())


for i in range(len(lsplit)):
    for j in range(len(lsplit[i])):
        lsplit[i][j]= float(lsplit[i][j])

node_num = int(lsplit[0][0])
ele_num = int(lsplit[0][1])
nodes_array = lsplit[1:node_num+1]
nodes = {}
for i in range(len(nodes_array)):
    k = int(nodes_array[i][0])
    nodes[k] = nodes_array[i][1:]

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
        

bc = numpy.array(bc)
ele = numpy.array(ele)

# convert b

def distance(node_1,node_2):
    return math.sqrt((node_2[0]-node_1[0])**2 + (node_2[1]-node_1[1])**2)

def midpoint(node_1,node_2):
    x1 = node_1[0]
    y1 = node_1[1]
    x2 = node_2[0]
    y2 = node_2[1]
    return [(x1 + x2)/2, (y1 + y2)/2]

def make_element_array(ele,BC,E,G):  
    """ make element array"""  
    elements = []
    for i in range(len(ele)):
        ele_no = ele[i,0]
        node_1 = nodes[ele[i,2]]
        node_2 = nodes[ele[i,3]]
        node_3 = nodes[ele[i,4]]
        node_4 = nodes[ele[i,5]]
        b = distance(node_1, node_2)
        a = distance(node_2, node_3)
        center = midpoint(node_1, node_3)
        x = center[0]
        y = center[1]
        elements.append(aem.Element(i, b, a, E, G, x, y,node_1,node_2,node_3,node_4))
        ele[i,0] = i
    return elements

# Creating element array
elements = make_element_array(ele,bc,2.07e9,79.3e9)

fig1 = plt.figure(1)
fig1.clf()

for e in elements:
    node_1 = e.node_1
    node_2 = e.node_2
    node_3 = e.node_3
    node_4 = e.node_4
    cords = numpy.vstack([node_1,node_2,node_3,node_4])
    plt.plot(cords[:,0],cords[:,1])
    plt.text(e.x,e.y,e.ele_no)

plt.axis('equal')
plt.show()

def both_in(a,b,x):
    """true if both a and b are elements of x"""
    return a in x and b in x

########

element_dic = {}
for i in elements:
    element_dic[i.ele_no] = i

no_spr = 40
ele_com = []
for i in range(ele.shape[0]):
    for k in range(i+1,ele.shape[0]):
        if both_in(ele[i,2], ele[i,3], ele[k,2:]):
            ele_com.append([element_dic[i],element_dic[k],1,3,no_spr])
        elif both_in(ele[i,3],ele[i,4], ele[k,2:]):
            ele_com.append([element_dic[i],element_dic[k],4,2,no_spr])
        elif both_in(ele[i,4],ele[i,5], ele[k,2:]):
            ele_com.append([element_dic[i],element_dic[k],3,1,no_spr])
        elif both_in(ele[i,5],ele[i,2], ele[k,2:]):
            ele_com.append([element_dic[i],element_dic[k],2,4,no_spr])
            
ele_com_mat = ele_com
############3

#f = open("prob2_con.dat")
#l = list(f)
#lsplit=[]
#for i in l:
#    lsplit.append(i.split())
#
#for i in range(len(lsplit)):
#    for j in range(len(lsplit[i])):
#        lsplit[i][j]= int(float(lsplit[i][j]))
#
#
#element_dic = {}
#for i in elements:
#    element_dic[i.ele_no] = i
#
#ele_com_mat = []
#
#for i in lsplit:
#    ele_1 = i[0]
#    ele_2 = i[1]
#    ele_com_mat.append([element_dic[ele_1],element_dic[ele_2],i[2],i[3],80])

ele_com_mat = numpy.array(ele_com_mat)

spr_com_mat = main.make_auto_spr_com_mat(ele_com_mat)

spr_array = main.make_spring_array(spr_com_mat)

glob_mat = main.make_global_stiff_mat(elements, spr_array)

n = len(elements)

all_dis = range((3*n))

# Boundary Condision
fixed_ele = [2,20,51]
fixed_dis = []
for x in fixed_ele:
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

loads = numpy.zeros((len(unknown_dis_list),1))
load_ele = [14]
for i in load_ele:
    loads[(our_map[i*3+1]),0] = -1000

unknown_dis = numpy.array([unknown_dis_list])
K=glob_mat[unknown_dis.T,unknown_dis]

dis = numpy.linalg.solve(K,loads)

for i in load_ele:
    print dis[(our_map[i*3+1]),0]