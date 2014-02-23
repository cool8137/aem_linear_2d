import math
import numpy
import aem

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
    
    spr_comb_mat = [{spr_no} {ele_1} {ele_2} {edg_1} {edg_2} {d} {alpha_1} {alpha_2} {L}]
    {spr_no} = spring no.
    {ele_1} = col. vector of first element no.
    {ele_2} = col. vector of second element no.
    {edg_1} = col. vector of first element's edge no.
    {edg_2} = col. vector of second element's edge no.
    {d} = thickness of spring
    {alpha_1} = angle towards spring contact wrt element one's vertical
    {alpha_2} = angle towards spring contact wrt element two's vertical
    {L_1} = length from center of element one to the contact point
    {L_2} = length from center of element two to the contact point
    """
    spring_array = []
    for i in range(len(spr_comb_mat)):
        spr_no = spr_comb_mat[i,0]
        ele_1 = spr_comb_mat[i,1]
        ele_2 = spr_comb_mat[i,2]
        edg_1 = spr_comb_mat[i,3]
        edg_2 = spr_comb_mat[i,4]
        d = spr_comb_mat[i,5]
        alpha_1 = spr_comb_mat[i,6]
        alpha_2 = spr_comb_mat[i,7]
        L_1 = spr_comb_mat[i,8]
        L_2 = spr_comb_mat[i,9]
        curr_spring = aem.Spring(spr_no,ele_1,ele_2,edg_1,edg_2,d,alpha_1,alpha_2,L_1,L_2)
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
    spr_com_mat = numpy.matrix([[0,0,0,0,0,0,0,0,0,0]])
    spr_no = 0

    for i in range(len(ele_com_mat)):
        ele_1 = ele_com_mat[i,0]
        ele_2 = ele_com_mat[i,1]
        edg_1 = ele_com_mat[i,2]
        edg_2 = ele_com_mat[i,3]
        n = ele_com_mat[i,4]
        
        # calc base and height for element one
        if (edg_1 == 1 or edg_1 == 3):
            height_1 = ele_1.a  # length of base (edge perp to contact edge)
            base_1 = ele_1.b    # length of contact edge
            
        elif (edg_1 == 2 or edg_1 == 4):
            height_1 = ele_1.b
            base_1 = ele_1.a
            
        else:
            raise AemError('edg_1 at row '+str(i)+' of ele_com_mat is invalid.')

        # calc base and height for element two
        if (edg_2 == 1 or edg_2 == 3):
            height_2 = ele_2.a  # length of base (edge perp to contact edge)
            base_2 = ele_2.b    # length of contact edge
            
        elif (edg_2 == 2 or edg_2 == 4):
            height_2 = ele_2.b
            base_2 = ele_2.a
            
        else:
             raise AemError('edg_2 at row '+str(i)+' of ele_com_mat is invalid.')

        # height of each spring
        d = 1.0 * height_1 / n

        # each edge will have n springs
        for j in range(n):
            pos_from_top_1 = (1.0*d/2) + (j*d)  #position of spring from top of element
            pos_from_top_2 = height_2 - pos_from_top_1 # postiton for ele 2 is compliment because ele two is mirro

            alpha_1 = pos_atan(base_1,(height_1/2.-pos_from_top_1))
            alpha_2 = pos_atan((base_2),(height_2/2.-pos_from_top_2)) # negative because ele two is a mirror to ele one

            L_1 = (base_1/2.) / math.sin(alpha_1)
            L_2 = (base_2/2.) / math.sin(alpha_2)

            spr_com_mat = numpy.vstack([spr_com_mat, [spr_no,ele_1,ele_2,edg_1,edg_2,d,alpha_1,alpha_2,L_1,L_2]])
            spr_no = spr_no + 1

    # delete first row that contains zeros only
    spr_com_mat = numpy.delete(spr_com_mat,(0),axis=0)
    return spr_com_mat


def make_global_stiff_mat(spr_array):
    pass
    
         
element = [aem.Element(1,4,5,2e9,2e9,0,0),aem.Element(2,4,5,2e9,2e9,5,0)]

ele_com_mat = numpy.matrix([[element[0],element[1],2,4,5]])

spr_com_mat = make_auto_spr_com_mat(ele_com_mat)

spr_array = make_spring_array(spr_com_mat)