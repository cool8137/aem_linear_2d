import math
import numpy as np
import aem
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix


class AemError():
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def pos_atan(y,x):
    """
    returns invers tangent with a positive angle
    """
    if (x == 0):
        return math.pi/2.
    elif (y/x > 0):
        return math.atan(y/x)
    else:
        return math.atan(y/x) + math.pi

def make_spring_array(spr_comb_mat):
    """
    Input: spring combination matrix.
    Output: an array of spring objects.
    
    spr_comb_mat = [{spr_no} {ele_1} {ele_2} {edg_1} {edg_2} {d} {x_1} {y_1} {x_2} {y_2}]
    {spr_no} = spring no.
    {ele_1} = col. vector of first element no.
    {ele_2} = col. vector of second element no.
    {edg_1} = col. vector of first element's edge no.
    {edg_2} = col. vector of second element's edge no.
    {d} = thickness of spring
    {x_1} = x coordinate of first point of spring
    {y_1} = y coordinate of first point of spring
    {x_2} = x coordinate of second point of spring
    {y_2} = y coordinate of second point of spring
    """
    spring_array = []
    for i in range(len(spr_comb_mat)):
        spr_no = spr_comb_mat[i,0]
        ele_1 = spr_comb_mat[i,1]
        ele_2 = spr_comb_mat[i,2]
        edg_1 = spr_comb_mat[i,3]
        edg_2 = spr_comb_mat[i,4]
        d = spr_comb_mat[i,5]
        x_1 = spr_comb_mat[i,6]
        y_1 = spr_comb_mat[i,7]
        x_2 = spr_comb_mat[i,8]
        y_2 = spr_comb_mat[i,9]
        curr_spring = aem.Spring(spr_no,ele_1,ele_2,edg_1,edg_2,d,x_1,y_1,x_2,y_2)
        spring_array.append(curr_spring)

    return spring_array


def make_auto_spr_com_mat(ele_com_mat):
    """
    Make auto generated springs when 
    Input: element combination matrix
    Output: spring combination matrix

    ele_comb_mat = [{ele_1} {ele_2} {edg_1} {edg_2} {n}]
    {ele_1} = col. vector of first element no.
    {ele_2} = col. vector of second element no.
    {edg_1} = col. vector of first element's edge no.
    {edg_2} = col. vector of second element's edge no.
    {n} = no. of springs per edge
    """

    # first row is defined so as to define the shape of matrix.
    # don't forget to delete the first row of spr_com_mat
    spr_com_mat = np.matrix([[0,0,0,0,0,0,0,0,0,0]])
    spr_no = 0

    for i in range(len(ele_com_mat)):
        ele_1 = ele_com_mat[i,0]
        ele_2 = ele_com_mat[i,1]
        edg_1 = ele_com_mat[i,2]
        edg_2 = ele_com_mat[i,3]
        n = ele_com_mat[i,4]
        
        ## check common edge of elements. if ele_2 has smaller edge than second ele will be ele_1
        #if (ele_com_mat[i,0].edg_len[ele_com_mat[i,2]] > ele_com_mat[i,1].edg_len[ele_com_mat[i,3]]):
        #    ele_2 = ele_com_mat[i,0]
        #    ele_1 = ele_com_mat[i,1]
        #    edg_2 = ele_com_mat[i,2]
        #    edg_1 = ele_com_mat[i,3]
        #    n = ele_com_mat[i,4]        
        #    
        #else:
        #    ele_1 = ele_com_mat[i,0]
        #    ele_2 = ele_com_mat[i,1]
        #    edg_1 = ele_com_mat[i,2]
        #    edg_2 = ele_com_mat[i,3]
        #    n = ele_com_mat[i,4]

        gov_height = ele_1.edg_len[edg_1] # it is the height of ele_1
        
        # height of each spring
        d = 1.0 * gov_height / n

        # each edge will have n springs
        for j in range(n):
            pos_from_center_1 = gov_height/2. - ((1.0*d/2) + (j*d)) #second term is pos_from_top_1
            
            # Calc angle (beta) from horizontal to the line connecting 
            # center of first element with the starting point of spring
            beta =   {2:ele_1.r,\
                     3:(ele_1.r+math.pi/2),\
                     4:(ele_1.r+math.pi),\
                     1:(ele_1.r+3.*math.pi/2.)}[edg_1] + math.pi/2
        
            #Calc coordinate for first point of spring
            x_1 = ele_1.x + pos_from_center_1 * math.cos(beta)
            y_1 = ele_1.y + pos_from_center_1 * math.sin(beta)

            # Calc dis from first point of spring to second (ie length of spring)
            spr_len = ele_1.edg_len[edg_1%2+1]/2. + ele_2.edg_len[edg_2%2+1]/2.

            # Calc coordinate of second point of spring
            x_2 = x_1 + spr_len * math.cos(beta-math.pi/2)
            y_2 = y_1 + spr_len * math.sin(beta-math.pi/2)
            
            # Stack all info into the spr_com_mat matrix
            spr_com_mat = np.vstack([spr_com_mat, [spr_no,ele_1,ele_2,edg_1,edg_2,d,x_1,y_1,x_2,y_2]])
            spr_no = spr_no + 1

    # delete first row that contains zeros only
    spr_com_mat = np.delete(spr_com_mat,(0),axis=0)
    return spr_com_mat

def make_BC_array(spr_comb_mat):
    """
    Input: BCspring combination matrix.
    Output: an array of BCspring objects.
    
    spr_comb_mat = [{spr_no} {ele_1} {edg_1} {BC} {d} {x_1} {y_1} {x_2} {y_2}]
    {spr_no} = spring no.
    {ele_1} = col. vector of first element no.
    {edg_1} = col. vector of first element's edge no.
    {BC} = BC.
    {d} = thickness of spring
    {x_1} = x coordinate of first point of spring
    {y_1} = y coordinate of first point of spring
    {x_2} = x coordinate of second point of spring
    {y_2} = y coordinate of second point of spring
    """
    spring_array = []
    for i in range(len(spr_comb_mat)):
        spr_no = spr_comb_mat[i,0]
        ele_1 = spr_comb_mat[i,1]
        edg_1 = spr_comb_mat[i,2]
        BC = spr_comb_mat[i,3]
        d = spr_comb_mat[i,4]
        x_1 = spr_comb_mat[i,5]
        y_1 = spr_comb_mat[i,6]
        x_2 = spr_comb_mat[i,7]
        y_2 = spr_comb_mat[i,8]
        curr_spring = aem.BCspring(spr_no,ele_1,edg_1,BC,d,x_1,y_1,x_2,y_2)
        spring_array.append(curr_spring)

    return spring_array

def make_auto_BC_com_mat(ele_com_mat):
    """
    Make auto generated BCsprings when 
    Input: element combination matrix
    Output: spring combination matrix

    ele_comb_mat = [{ele_1} {edg_1} {BC} {n}]
    {ele_1} = col. vector of first element no.
    {edg_1} = col. vector of first element's edge no.
    {BC} = BC
    {n} = no. of springs per edge
    """

    # first row is defined so as to define the shape of matrix.
    # don't forget to delete the first row of spr_com_mat
    spr_com_mat = np.matrix([[0,0,0,0,0,0,0,0,0]])
    spr_no = 0

    for i in range(len(ele_com_mat)):
        ele_1 = ele_com_mat[i,0]
        edg_1 = ele_com_mat[i,1]
        BC = ele_com_mat[i,2]
        n = ele_com_mat[i,3]

        gov_height = ele_1.edg_len[edg_1] # it is the height of ele_1
        
        # height of each spring
        d = 1.0 * gov_height / n

        # each edge will have n springs
        for j in range(n):
            pos_from_center_1 = gov_height/2. - ((1.0*d/2) + (j*d)) #second term is pos_from_top_1
            
            # Calc angle (beta) from horizontal to the line connecting 
            # center of first element with the starting point of spring
            beta =   {2:ele_1.r,\
                     3:(ele_1.r+math.pi/2),\
                     4:(ele_1.r+math.pi),\
                     1:(ele_1.r+3.*math.pi/2.)}[edg_1] + math.pi/2
        
            #Calc coordinate for first point of spring
            x_1 = ele_1.x + pos_from_center_1 * math.cos(beta)
            y_1 = ele_1.y + pos_from_center_1 * math.sin(beta)

            # Calc dis from first point of spring to edge (ie length of spring)
            spr_len = ele_1.edg_len[edg_1%2+1]/2.

            # Calc coordinate of second point of spring (at the center of edge)
            x_2 = x_1 + spr_len * math.cos(beta-math.pi/2)
            y_2 = y_1 + spr_len * math.sin(beta-math.pi/2)
            
            # Stack all info into the spr_com_mat matrix
            spr_com_mat = np.vstack([spr_com_mat, [spr_no,ele_1,edg_1,BC,d,x_1,y_1,x_2,y_2]])
            spr_no = spr_no + 1

    # delete first row that contains zeros only
    spr_com_mat = np.delete(spr_com_mat,(0),axis=0)
    return spr_com_mat


def make_global_stiff_mat(ele_array,spr_array,BC_array=[]):
    """
    numering: 
    0 = element one horizontal
    1 = element one vertical
    2 = element one rotation
    3 = element two horizontal
    4 = element two vertical
    5 = element two rotation
    6 = element three horizaontal
    ...
    12 = three element one horizaontal
    ...
    """
    N = len(ele_array)
    #global_mat = lil_matrix((3*N,3*N))
    global_mat = np.zeros((3*N,3*N))

    for i in range(len(spr_array)):
        spr = spr_array[i]

        spr_mat = spr.stiffness_matrix()
        #print spr_mat.shape
        
        ele_i = spr.ele_1.ele_no
        ele_j = spr.ele_2.ele_no
        
        # Location vector
        lv = np.array([[(ele_i*3), (ele_i*3+1), (ele_i*3+2),\
                        (ele_j*3), (ele_j*3+1), (ele_j*3+2)]])
        
        global_mat[lv.T,lv] = global_mat[lv.T,lv] + spr_mat

        # Note while extracting array from array
        # A[[[r1],[r2],[r3]],[c1,c2,c3,c4]]

    for i in range(len(BC_array)):
        BC = BC_array[i]

        BC_mat = BC.stiffness_matrix()
        #print spr_mat.shape
        
        ele_i = spr.ele_1.ele_no
        
        # Location vector
        lv = np.array([[(ele_i*3), (ele_i*3+1), (ele_i*3+2)]])
        
        global_mat[lv.T,lv] = global_mat[lv.T,lv] + BC_mat

        # Note while extracting array from array
        # A[[[r1],[r2],[r3]],[c1,c2,c3,c4]]
    
    return global_mat
        
def main():
    element = [aem.Element(0, 0.2, 0.1, 2.07e9, 79.3e9, 0.000000e+000, 0.000000e+000),\
            aem.Element(1, 0.2, 0.1, 2.07e9, 79.3e9, 1.000000e-001, 0.000000e+000),\
            aem.Element(2, 0.2, 0.1, 2.07e9, 79.3e9, 2.000000e-001, 0.000000e+000),\
            aem.Element(3, 0.2, 0.1, 2.07e9, 79.3e9, 3.000000e-001, 0.000000e+000),\
            aem.Element(4, 0.2, 0.1, 2.07e9, 79.3e9, 4.000000e-001, 0.000000e+000),\
            aem.Element(5, 0.2, 0.1, 2.07e9, 79.3e9, 5.000000e-001, 0.000000e+000),\
            aem.Element(6, 0.2, 0.1, 2.07e9, 79.3e9, 6.000000e-001, 0.000000e+000),\
            aem.Element(7, 0.2, 0.1, 2.07e9, 79.3e9, 7.000000e-001, 0.000000e+000),\
            aem.Element(8, 0.2, 0.1, 2.07e9, 79.3e9, 8.000000e-001, 0.000000e+000),\
            aem.Element(9, 0.2, 0.1, 2.07e9, 79.3e9, 9.000000e-001, 0.000000e+000),\
            aem.Element(10, 0.2, 0.1, 2.07e9, 79.3e9, 1.000000e+000, 0.000000e+000),\
            aem.Element(11, 0.2, 0.1, 2.07e9, 79.3e9, 1.100000e+000, 0.000000e+000),\
            aem.Element(12, 0.2, 0.1, 2.07e9, 79.3e9, 1.200000e+000, 0.000000e+000),\
            aem.Element(13, 0.2, 0.1, 2.07e9, 79.3e9, 1.300000e+000, 0.000000e+000),\
            aem.Element(14, 0.2, 0.1, 2.07e9, 79.3e9, 1.400000e+000, 0.000000e+000),\
            aem.Element(15, 0.2, 0.1, 2.07e9, 79.3e9, 1.500000e+000, 0.000000e+000),\
            aem.Element(16, 0.2, 0.1, 2.07e9, 79.3e9, 1.600000e+000, 0.000000e+000),\
            aem.Element(17, 0.2, 0.1, 2.07e9, 79.3e9, 1.700000e+000, 0.000000e+000),\
            aem.Element(18, 0.2, 0.1, 2.07e9, 79.3e9, 1.800000e+000, 0.000000e+000),\
            aem.Element(19, 0.2, 0.1, 2.07e9, 79.3e9, 1.900000e+000, 0.000000e+000),\
            aem.Element(20, 0.2, 0.1, 2.07e9, 79.3e9, 2.00000e+000, 0.000000e+000)]
            
    ele_com_mat = np.array([[element[0],element[1],2,4,20],\
                            [element[1],element[2],2,4,20],\
                            [element[2],element[3],2,4,20],\
                            [element[3],element[4],2,4,20],\
                            [element[4],element[5],2,4,20],\
                            [element[5],element[6],2,4,20],\
                            [element[6],element[7],2,4,20],\
                            [element[7],element[8],2,4,20],\
                            [element[8],element[9],2,4,20],\
                            [element[9],element[10],2,4,20],\
                            [element[10],element[11],2,4,20],\
                            [element[11],element[12],2,4,20],\
                            [element[12],element[13],2,4,20],\
                            [element[13],element[14],2,4,20],\
                            [element[14],element[15],2,4,20],\
                            [element[15],element[16],2,4,20],\
                            [element[16],element[17],2,4,20],\
                            [element[17],element[18],2,4,20],\
                            [element[18],element[19],2,4,20],\
                            [element[19],element[20],2,4,20]])
    
    spr_com_mat = make_auto_spr_com_mat(ele_com_mat)
    
    spr_array = make_spring_array(spr_com_mat)
    
    BC_com_mat = make_auto_BC_com_mat(np.array([[element[0],4,3,20]]))
    
    BC_array = make_BC_array(BC_com_mat)
    
    glob_mat = make_global_stiff_mat(element, spr_array,BC_array)
    
    n = len(element)
    
    #unknown_dis = np.array([range(3,(3*n))])
    
    loads = np.hstack([np.zeros((1,(3*n-3))),[[0,-1000,0]]]).T
    
    K=glob_mat
    
    dis43 = np.linalg.solve(K,loads)
    
    #####
    
    element2 = [aem.Element(0, 0.3, 300, 2.07e9, 79.3e9, 0.000000e+000, 0.000000e+000),\
            aem.Element(1, 0.3, 300, 2.07e9, 79.3e9, 3.0000e-001, 0.000000e+000)]
    ele_com_mat2 = np.array([[element2[0],element2[1],2,4,20]])
    
    spr_com_mat2 = make_auto_spr_com_mat(ele_com_mat2)
    
    spr_array2 = make_spring_array(spr_com_mat2)
    
    glob_mat2 = make_global_stiff_mat(element2, spr_array2)
    
    n2 = len(element2)
    
    unknown_dis2 = np.array([[3,4,5]])
    
    loads2 = np.array([[10000,0,0]]).T
    
    K2=glob_mat2[unknown_dis2.T,unknown_dis2]
    
    dis2 = np.linalg.solve(K2,loads2)
    
    #####
    
    element3 = [aem.Element(0, 0.2, 0.2, 2.07e9, 79.3e9, 0.000000e+000, 0.000000e+000),\
            aem.Element(1, 0.2, 0.2, 2.07e9, 79.3e9, 2.000000e-001, 0.000000e+000)]
    ele_com_mat3 = np.array([[element3[0],element3[1],2,4,20]])
    
    spr_com_mat3 = make_auto_spr_com_mat(ele_com_mat3)
    
    spr_array3 = make_spring_array(spr_com_mat3)
    
    BC_com_mat3 = make_auto_BC_com_mat(np.array([[element3[0],4,3,2]]))
    
    BC_array3 = make_BC_array(BC_com_mat3)
    
    glob_mat3 = make_global_stiff_mat(element3, spr_array3, BC_array3)
    
    n3 = len(element2)
    
    loads3 = np.array([[0,0,0,1000,0,0]]).T
    
    K3=glob_mat3
    
    dis3 = np.linalg.solve(K3,loads3)
    
    
    element = [aem.Element(0, 0.2, 0.1, 2.07e9, 79.3e9, 0.000000e+000, 0.000000e+000),\
            aem.Element(1, 0.2, 0.1, 2.07e9, 79.3e9, 1.000000e-001, 0.000000e+000),\
            aem.Element(2, 0.2, 0.1, 2.07e9, 79.3e9, 2.000000e-001, 0.000000e+000),\
            aem.Element(3, 0.2, 0.1, 2.07e9, 79.3e9, 3.000000e-001, 0.000000e+000),\
            aem.Element(4, 0.2, 0.1, 2.07e9, 79.3e9, 4.000000e-001, 0.000000e+000),\
            aem.Element(5, 0.2, 0.1, 2.07e9, 79.3e9, 5.000000e-001, 0.000000e+000),\
            aem.Element(6, 0.2, 0.1, 2.07e9, 79.3e9, 6.000000e-001, 0.000000e+000),\
            aem.Element(7, 0.2, 0.1, 2.07e9, 79.3e9, 7.000000e-001, 0.000000e+000),\
            aem.Element(8, 0.2, 0.1, 2.07e9, 79.3e9, 8.000000e-001, 0.000000e+000),\
            aem.Element(9, 0.2, 0.1, 2.07e9, 79.3e9, 9.000000e-001, 0.000000e+000),\
            aem.Element(10, 0.2, 0.1, 2.07e9, 79.3e9, 1.000000e+000, 0.000000e+000),\
            aem.Element(11, 0.2, 0.1, 2.07e9, 79.3e9, 1.100000e+000, 0.000000e+000),\
            aem.Element(12, 0.2, 0.1, 2.07e9, 79.3e9, 1.200000e+000, 0.000000e+000),\
            aem.Element(13, 0.2, 0.1, 2.07e9, 79.3e9, 1.300000e+000, 0.000000e+000),\
            aem.Element(14, 0.2, 0.1, 2.07e9, 79.3e9, 1.400000e+000, 0.000000e+000),\
            aem.Element(15, 0.2, 0.1, 2.07e9, 79.3e9, 1.500000e+000, 0.000000e+000),\
            aem.Element(16, 0.2, 0.1, 2.07e9, 79.3e9, 1.600000e+000, 0.000000e+000),\
            aem.Element(17, 0.2, 0.1, 2.07e9, 79.3e9, 1.700000e+000, 0.000000e+000),\
            aem.Element(18, 0.2, 0.1, 2.07e9, 79.3e9, 1.800000e+000, 0.000000e+000),\
            aem.Element(19, 0.2, 0.1, 2.07e9, 79.3e9, 1.900000e+000, 0.000000e+000),\
            aem.Element(20, 0.2, 0.1, 2.07e9, 79.3e9, 2.00000e+000, 0.000000e+000)]
            
    ele_com_mat = np.array([[element[0],element[1],2,4,80],\
                            [element[1],element[2],2,4,80],\
                            [element[2],element[3],2,4,80],\
                            [element[3],element[4],2,4,80],\
                            [element[4],element[5],2,4,80],\
                            [element[5],element[6],2,4,80],\
                            [element[6],element[7],2,4,80],\
                            [element[7],element[8],2,4,80],\
                            [element[8],element[9],2,4,80],\
                            [element[9],element[10],2,4,80],\
                            [element[10],element[11],2,4,80],\
                            [element[11],element[12],2,4,80],\
                            [element[12],element[13],2,4,80],\
                            [element[13],element[14],2,4,80],\
                            [element[14],element[15],2,4,80],\
                            [element[15],element[16],2,4,80],\
                            [element[16],element[17],2,4,80],\
                            [element[17],element[18],2,4,80],\
                            [element[18],element[19],2,4,80],\
                            [element[19],element[20],2,4,80]])
    
    spr_com_mat = make_auto_spr_com_mat(ele_com_mat)
    
    spr_array = make_spring_array(spr_com_mat)
    
    glob_mat = make_global_stiff_mat(element, spr_array)
    
    n = len(element)
    
    unknown_dis = np.array([range(3,(3*n))])
    
    loads = np.hstack([np.zeros((1,(3*n-6))),[[0,-1000,0]]]).T
    
    K=glob_mat[unknown_dis.T,unknown_dis]
    
    dis = np.linalg.solve(K,loads)
    
    #if __name__ == '__main__':
    #    main()
    
    element = [aem.Element(3, 0.2, 0.2, 2.07e9, 79.3e9, 0.000000e+000, 0.000000e+000),\
            aem.Element(6, 0.2, 0.2, 2.07e9, 79.3e9, 2.000000e-001, 0.000000e+000),\
            aem.Element(0, 0.2, 0.2, 2.07e9, 79.3e9, 4.000000e-001, 0.000000e+000),\
            aem.Element(5, 0.2, 0.2, 2.07e9, 79.3e9, 6.000000e-001, 0.000000e+000),\
            aem.Element(9, 0.2, 0.2, 2.07e9, 79.3e9, 8.000000e-001, 0.000000e+000),\
            aem.Element(8, 0.2, 0.2, 2.07e9, 79.3e9, 1.000000e+000, 0.000000e+000),\
            aem.Element(2, 0.2, 0.2, 2.07e9, 79.3e9, 1.200000e+000, 0.000000e+000),\
            aem.Element(1, 0.2, 0.2, 2.07e9, 79.3e9, 1.400000e+000, 0.000000e+000),\
            aem.Element(10, 0.2, 0.2, 2.07e9, 79.3e9, 1.600000e+000, 0.000000e+000),\
            aem.Element(4, 0.2, 0.2, 2.07e9, 79.3e9, 1.800000e+000, 0.000000e+000),\
            aem.Element(7, 0.2, 0.2, 2.07e9, 79.3e9, 2.00000e+000, 0.000000e+000)]
            
    ele_com_mat = np.array([[element[3],element[6],2,4,80],\
                            [element[6],element[0],2,4,80],\
                            [element[0],element[5],2,4,80],\
                            [element[5],element[9],2,4,80],\
                            [element[9],element[8],2,4,80],\
                            [element[8],element[2],2,4,80],\
                            [element[2],element[1],2,4,80],\
                            [element[1],element[10],2,4,80],\
                            [element[10],element[4],2,4,80],\
                            [element[4],element[7],2,4,80]])
    
    spr_com_mat = make_auto_spr_com_mat(ele_com_mat)
    
    spr_array = make_spring_array(spr_com_mat)
    
    glob_mat = make_global_stiff_mat(element, spr_array)
    
    n = len(element)
    
    unknown_dis = np.array([range(3,(3*n))])
    
    loads = np.hstack([np.zeros((1,(3*n-6))),[[0,-1000,0]]]).T
    
    # Boundary Condision
    fixed_ele = [3]
    fixed_dis = []
    for x in fixed_ele:
        fixed_dis.append(3*x)
        fixed_dis.append(3*x+1)
        fixed_dis.append(3*x+2)
    
    all_dis = range((3*n))
    
    unknown_dis_list = []
    for x in all_dis:
        if not(x in fixed_dis):
            unknown_dis_list.append(x)
    
    loads = np.zeros((len(unknown_dis_list),1))
    load_ele = [7]
    for i in load_ele:
        loads[i*3+1,0] = -1000
    
    K=glob_mat[unknown_dis.T,unknown_dis]
    
    dis = np.linalg.solve(K,loads)
    
    for i in load_ele:
        print dis[i*3+1,0]

if __name__ == '__main__':
    main()