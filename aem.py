# list of classes used by the Applied Element Method (AEM) program
import math
import numpy

class Element:
    """Element class
    ele_no = element number
    edg_len = {1:a,2:b,3:a,4:b}
    b = height of element (edge 2 and 4)
    a = length of element (edge 1 and 3)
    T = thickness of element
    E = Young's modulus
    G = Shear modulus
    x = x coordinate of center of element
    y = y coordinate of center of element
    r = rotation angle in radian of element
    [note: when r is zero, edg no. 2 is towards the positive x axis]
    """
    def __init__(self,ele_no,b,a,E,G,x,y,r=0,T=0.15):
        self.ele_no = ele_no
        self.edg_len = {1:a,2:b,3:a,4:b,5:a,6:b}
        self.E = E
        self.G = G
        self.T = T
        self.x = x
        self.y = y
        self.r = r
        

class Spring:
    """Spring class
    ele_1 = element one
    ele_2 = element two
    edg_1 = edge no. of element one
    edg_2 = edge no. of element two
    d = height of spring
    T = thickness of spring
    K_n = Normal spring constant
    K_s = Shear spring constant
    x_1 = x coordinate of first point of spring
    y_1 = y coordinate of first point of spring
    x_2 = x coordinate of second point of spring
    y_2 = y coordinate of second point of spring
    """

    def __init__(self,spr_no, ele_1,ele_2,edg_1,edg_2,d,x_1,y_1,x_2,y_2):

        self.spr_no = spr_no
        
        self.ele_1 = ele_1
        self.ele_2 = ele_2
        
        self.edg_1 = edg_1
        self.edg_2 = edg_2

        self.d = d
        T = min(ele_1.T,ele_2.T)
        self.T = T

        K_n_1 = 1.0*(ele_1.E * d * T)/(self.ele_1.edg_len[edg_1+1]/2.)
        K_n_2 = 1.0*(ele_2.E * d * T)/(self.ele_2.edg_len[edg_2+1]/2.)
        K_n = K_n_1*K_n_2/(K_n_1 + K_n_2)
        self.K_n = K_n

        K_s_1 = 1.0*(ele_1.G * d * T)/(self.ele_1.edg_len[edg_1+1]/2.)
        K_s_2 = 1.0*(ele_2.G * d * T)/(self.ele_2.edg_len[edg_2+1]/2.)
        K_s = K_s_1*K_s_2/(K_s_1 + K_s_2)
        self.K_s = K_s

        self.x_1 = x_1
        self.y_1 = y_1        
        self.x_2 = x_2
        self.y_2 = y_2        


    def stiffness_matrix(self):
        """
        Returns spring stiffness matrix (6X6)
        """
        mat = numpy.matrix(numpy.zeros((6,6)))
        
        L = math.sqrt((self.x_2 - self.x_1)**2 + (self.y_2 - self.y_1)**2)
        l = (self.x_2 - self.x_1)/L
        m = (self.y_2 - self.y_1)/L
        
        # Distance from center of ele_1 to first point of spring (d11)
        dx11 = (self.x_1 - self.ele_1.x)
        dy11 = (self.y_1 - self.ele_1.y)
        d11 = math.sqrt(dx11**2 + dy11**2)
        
        # Distance from center of ele_2 to second point of spring (d22)
        dx22 = (self.x_2 - self.ele_2.x)
        dy22 = (self.y_2 - self.ele_2.y)
        d22 = math.sqrt(dx22**2 + dy22**2)
        
        # Perp. distance from center of ele_1 to the common edge (d1e)
        d1e = self.ele_1.edg_len[self.edg_1 + 1]/2 # using mod 2 +1 to get the base length
        
        # Perp. distance from center of ele_2 to the common edge (d2e)
        d2e = self.ele_2.edg_len[self.edg_2%2 + 1]/2

        # top left (3X3)
        # unit horizontal disp on ele 1
        mat[0,0] = (self.K_n * l * l) + (self.K_s * m * m) # horizontal reaction at ele 1
        mat[0,1] = (self.K_n * l * m) - (self.K_s * m * l) # vertical reaction at ele 1
        mat[0,2] = (-self.K_n * l * l * dy11 - self.K_n * l * m * dx11) - (self.K_s * m * d1e) # moment reaction at ele 1
        # unit vertical disp on ele 1
        mat[1,0] = mat[0,1] # hor reaction at ele 1
        mat[1,1] = (self.K_n * m * m) + (self.K_s * l * l) # ver reaction at ele 1
        mat[1,2] = (-self.K_n * m * l * dy11 - self.K_n * m * m * dx11) + (self.K_s * l * d1e) # moment reaction at ele 1
        # unit rotation on ele 1
        mat[2,0] = mat[0,2] # hor reaction at ele 1
        mat[2,1] = mat[1,2] # ver reaction at ele 1
        mat[2,2] = (self.K_n * d11 * d11) + (self.K_s * d1e * d1e) # moment reaction at ele 1

        # top right (3X3)
        # unit hor disp on ele 1
        mat[0,3] = -((self.K_n * l * l) + (self.K_s * m * m)) # hor reaction at ele 2
        mat[0,4] = -((self.K_n * l * m) - (self.K_s * m * l)) # ver reaction at ele 2
        mat[0,5] = -((-self.K_n * l * l * dy22 - self.K_n * l * m * dx22) + (self.K_s * m * d2e)) #1 rotation reaction on ele_2
        # unit ver disp on ele 1
        mat[1,3] = mat[0,4] # hor reaction at ele 2
        mat[1,4] = -((self.K_n * m * m) + (self.K_s * l * l)) # ver reaction at ele 2
        mat[1,5] = -((-self.K_n * m * l * dy22 - self.K_n * m * m * dx22) - (self.K_s * l * d2e)) #2 rotation reaction on ele_2
        # unit rotation on ele 1
        mat[2,3] = -((-self.K_n * l * l * dy11 - self.K_n * l * m * dx11) - (self.K_s * m * d1e)) #1 rotation on ele_1 causing  horizontal reaction at two
        mat[2,4] = -((-self.K_n * m * l * dy11 - self.K_n * m * m * dx11) + (self.K_s * l * d1e)) #2 rotation on ele_1 causing vertical reaction at two
        mat[2,5] = -((self.K_n * d11 * d22) - (self.K_s * d1e * d2e)) # roation on ele_1 casing reation on ele_2

        # bottom left (3X3)
        mat[3,0] = mat[0,3]
        mat[3,1] = mat[1,3]
        mat[3,2] = mat[2,3]
        mat[4,0] = mat[0,4]
        mat[4,1] = mat[1,4]
        mat[4,2] = mat[2,4]
        mat[5,0] = mat[0,5]
        mat[5,1] = mat[1,5]
        mat[5,2] = mat[2,5]

        # bottom right (3X3)
        mat[3,3] = mat[0,0]
        mat[3,4] = mat[0,1]
        mat[3,5] = (-self.K_n * l * l * dy22 - self.K_n * l * m * dx22) - (self.K_s * m * d2e)
        mat[4,3] = mat[1,0]
        mat[4,4] = mat[1,1]
        mat[4,5] = (-self.K_n * m * l * dy22 - self.K_n * m * m * dx22) + (self.K_s * l * d2e)
        mat[5,3] = mat[3,5]
        mat[5,4] = mat[4,5]
        mat[5,5] = (self.K_n * d22 * d22) + (self.K_s * d2e * d2e)

##################################

#        # top left (3X3)
#        mat[0,0] = (self.K_n * math.sin(angle_sum)**2)\
#                    + (self.K_s * math.cos(angle_sum)**2)
#        mat[0,1] = (-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
#                    + (self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
#        mat[0,2] = (-self.K_n * self.L_1 * math.sin(angle_sum) * math.cos(self.alpha_1))\
#                    + (self.K_s * self.L_1 * math.cos(angle_sum) * math.sin(self.alpha_1))
#        mat[1,0] = mat[0,1]
#        mat[1,1] = (self.K_n * math.cos(angle_sum)**2) \
#                    + (self.K_s * math.sin(angle_sum)**2)
#        mat[1,2] = (self.K_n * self.L_1 * math.cos(angle_sum) * math.cos(self.alpha_1))\
#                    + (self.K_s * self.L_1 * math.sin(angle_sum) * math.sin(self.alpha_1))
#        mat[2,0] = mat[0,2]
#        mat[2,1] = mat[1,2]
#        mat[2,2] = (self.K_n * self.L_1**2 * math.cos(self.alpha_1)**2)\
#                    + (self.K_s * self.L_1**2 * math.sin(self.alpha_1)**2)
#
#        # top right (3X3)
#        mat[0,3] = -(self.K_n * math.sin(angle_sum)**2)\
#                    - (self.K_s * math.cos(angle_sum)**2)
#        mat[0,4] = -(-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
#                    + -(self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
#        mat[0,5] = -(-self.K_n * self.L_2 * math.sin(angle_sum) * math.cos(self.alpha_2))\
#                    + -(self.K_s * self.L_2 * math.cos(angle_sum) * math.sin(self.alpha_2))
#        mat[1,3] = mat[0,4]
#        mat[1,4] = -(self.K_n * math.cos(angle_sum)**2) \
#                    + -(self.K_s * math.sin(angle_sum)**2)
#        mat[1,5] = -(self.K_n * self.L_2 * math.cos(angle_sum) * math.cos(self.alpha_1))\
#                    + -(self.K_s * self.L_2 * math.sin(angle_sum) * math.sin(self.alpha_1))
#        mat[2,3] = mat[0,5]
#        mat[2,4] = mat[1,5]
#        mat[2,5] = -(self.K_n * self.L_2**2 * math.cos(self.alpha_2)**2)\
#                    + -(self.K_s * self.L_2**2 * math.sin(self.alpha_2)**2)
#
#        # bottom left (3X3)
#        mat[3,0] = -(self.K_n * math.sin(angle_sum)**2)\
#                    - (self.K_s * math.cos(angle_sum)**2)
#        mat[3,1] = -(-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
#                    + -(self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
#        mat[3,2] = -(-self.K_n * self.L_2 * math.sin(angle_sum) * math.cos(self.alpha_2))\
#                    + -(self.K_s * self.L_2 * math.cos(angle_sum) * math.sin(self.alpha_2))
#        mat[4,0] = mat[0,4]
#        mat[4,1] = -(self.K_n * math.cos(angle_sum)**2) \
#                    + -(self.K_s * math.sin(angle_sum)**2)
#        mat[4,2] = -(self.K_n * self.L_2 * math.cos(angle_sum) * math.cos(self.alpha_1))\
#                    + -(self.K_s * self.L_2 * math.sin(angle_sum) * math.sin(self.alpha_1))
#        mat[5,0] = mat[0,5]
#        mat[5,1] = mat[1,5]
#        mat[5,2] = -(self.K_n * self.L_2**2 * math.cos(self.alpha_2)**2)\
#                    + -(self.K_s * self.L_2**2 * math.sin(self.alpha_2)**2)
#
#        # bottom right (3X3)
#        mat[3,3] = (self.K_n * math.sin(angle_sum)**2)\
#                    + (self.K_s * math.cos(angle_sum)**2)
#        mat[3,4] = (-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
#                    + (self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
#        mat[3,5] = (-self.K_n * self.L_1 * math.sin(angle_sum) * math.cos(self.alpha_1))\
#                    + (self.K_s * self.L_1 * math.cos(angle_sum) * math.sin(self.alpha_1))
#        mat[4,3] = mat[0,1]
#        mat[4,4] = (self.K_n * math.cos(angle_sum)**2) \
#                    + (self.K_s * math.sin(angle_sum)**2)
#        mat[4,5] = (self.K_n * self.L_1 * math.cos(angle_sum) * math.cos(self.alpha_1))\
#                    + (self.K_s * self.L_1 * math.sin(angle_sum) * math.sin(self.alpha_1))
#        mat[5,3] = mat[0,2]
#        mat[5,4] = mat[1,2]
#        mat[5,5] = (self.K_n * self.L_1**2 * math.cos(self.alpha_1)**2)\
#                    + (self.K_s * self.L_1**2 * math.sin(self.alpha_1)**2)

        return mat

        
print "aem classes imported"