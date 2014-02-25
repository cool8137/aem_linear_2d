import numpy

f = open("smallMesh2d.dat",'r')
l = list(f)
lsplit=[]
for i in l:
    lsplit.append(i.split())


for i in range(len(lsplit)):
    for j in range(len(lsplit[i])):
        lsplit[i][j]= float(lsplit[i][j])

node_num = int(lsplit[0][0])
ele_num = int(lsplit[0][1])
nodes = numpy.array(lsplit[1:node_num+1])

remaining =lsplit[1+node_num:1+node_num+ele_num]
bc = []
ele =[]
for i in remaining:
    if int(i[1]) == 102:
        bc.append(i)
    else:
        ele.append(i)

bc = numpy.array(bc)
ele = numpy.array(ele)


