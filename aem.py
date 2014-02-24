# list of classes used by the Applied Element Method (AEM) program
import math
import numpy

class Element:
    """Element class
    ele_no = element number
    edg_len = edg length
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
    def __init__(self,ele_no,b,a,E,G,x,y,r=0,T=1):
        self.ele_no = ele_no
        self.edg_len = {1:a,2:b,3:a,4:b}
        self.b = b
        self.a = a
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
    theta_1 = angle to spring contact wrt horizontal
    theta_2 = angle to spring contact wrt horizontal
    alpha_1 = angle to spring contact wrt element one's vertical
    alpha_2 = angle to spring contact wrt element two's vertical
    L_1 = length from center of element one to the contact point
    L_2 = length from center to element two to the contact point
    """

    def __init__(self,spr_no, ele_1,ele_2,edg_1,edg_2,d,alpha_1,alpha_2,L_1,L_2):

        self.spr_no = spr_no
        
        self.ele_1 = ele_1
        self.ele_2 = ele_2
        
        self.edg_1 = edg_1
        self.edg_2 = edg_2

        self.d = d
        T = min(ele_1.T,ele_2.T)
        self.T = T

        K_n_1 = 1.0*(ele_1.E * d * T)/(self.ele_1.a/2.)
        K_n_2 = 1.0*(ele_2.E * d * T)/(self.ele_2.a/2.)
        K_n = K_n_1*K_n_2/(K_n_1 + K_n_2)
        self.K_n = K_n

        K_s_1 = 1.0*(ele_1.G * d * T)/(self.ele_1.a/2.)
        K_s_2 = 1.0*(ele_2.G * d * T)/(self.ele_2.a/2.)
        K_s = K_s_1*K_s_2/(K_s_1 + K_s_2)
        self.K_s = K_s

        def calc_theta(edg, alpha, r):
            if (edg == 2):
                return (alpha - math.pi/2.) + ele_1.r
            elif(edg == 3):
                return (alpha - math.pi/2.) + ele_1.r + math.pi / 2.
            elif(edg == 4):
                return (alpha - math.pi/2.) + ele_1.r + math.pi
            elif(edg == 1):
                return (alpha - math.pi/2.) + ele_1.r + 3. * math.pi / 2.
            else:
                 raise AemError('edg_2 at row '+str(i)+' of ele_com_mat is invalid.')
            
        self.theta_1 = calc_theta(edg_1, alpha_1, ele_1.r)
        self.theta_2 = calc_theta(edg_2, alpha_2, ele_2.r)

        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2

        self.L_1 = L_1
        self.L_2 = L_2


    def stiffness_matrix(self):
        """
        Returns spring stiffness matrix (6X6)
        """
        mat = numpy.matrix(numpy.zeros((6,6)))
        
        angle_sum = self.theta_1 + self.alpha_1

        # top left (3X3)
        mat[0,0] = (self.K_n * math.sin(angle_sum)**2)\
                    + (self.K_s * math.cos(angle_sum)**2)
        mat[0,1] = (-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
                    + (self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
        mat[0,2] = (-self.K_n * self.L_1 * math.sin(angle_sum) * math.cos(self.alpha_1))\
                    + (self.K_s * self.L_1 * math.cos(angle_sum) * math.sin(self.alpha_1))
        mat[1,0] = mat[0,1]
        mat[1,1] = (self.K_n * math.cos(angle_sum)**2) \
                    + (self.K_s * math.sin(angle_sum)**2)
        mat[1,2] = (self.K_n * self.L_1 * math.cos(angle_sum) * math.cos(self.alpha_1))\
                    + (self.K_s * self.L_1 * math.sin(angle_sum) * math.sin(self.alpha_1))
        mat[2,0] = mat[0,2]
        mat[2,1] = mat[1,2]
        mat[2,2] = (self.K_n * self.L_1**2 * math.cos(self.alpha_1)**2)\
                    + (self.K_s * self.L_1**2 * math.sin(self.alpha_1)**2)

        # top right (3X3)
        mat[0,3] = -(self.K_n * math.sin(angle_sum)**2)\
                    - (self.K_s * math.cos(angle_sum)**2)
        mat[0,4] = -(-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
                    + -(self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
        mat[0,5] = -(-self.K_n * self.L_2 * math.sin(angle_sum) * math.cos(self.alpha_2))\
                    + -(self.K_s * self.L_2 * math.cos(angle_sum) * math.sin(self.alpha_2))
        mat[1,3] = mat[0,4]
        mat[1,4] = -(self.K_n * math.cos(angle_sum)**2) \
                    + -(self.K_s * math.sin(angle_sum)**2)
        mat[1,5] = -(self.K_n * self.L_2 * math.cos(angle_sum) * math.cos(self.alpha_1))\
                    + -(self.K_s * self.L_2 * math.sin(angle_sum) * math.sin(self.alpha_1))
        mat[2,3] = mat[0,5]
        mat[2,4] = mat[1,5]
        mat[2,5] = -(self.K_n * self.L_2**2 * math.cos(self.alpha_2)**2)\
                    + -(self.K_s * self.L_2**2 * math.sin(self.alpha_2)**2)

        # bottom left (3X3)
        mat[3,0] = -(self.K_n * math.sin(angle_sum)**2)\
                    - (self.K_s * math.cos(angle_sum)**2)
        mat[3,1] = -(-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
                    + -(self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
        mat[3,2] = -(-self.K_n * self.L_2 * math.sin(angle_sum) * math.cos(self.alpha_2))\
                    + -(self.K_s * self.L_2 * math.cos(angle_sum) * math.sin(self.alpha_2))
        mat[4,0] = mat[0,4]
        mat[4,1] = -(self.K_n * math.cos(angle_sum)**2) \
                    + -(self.K_s * math.sin(angle_sum)**2)
        mat[4,2] = -(self.K_n * self.L_2 * math.cos(angle_sum) * math.cos(self.alpha_1))\
                    + -(self.K_s * self.L_2 * math.sin(angle_sum) * math.sin(self.alpha_1))
        mat[5,0] = mat[0,5]
        mat[5,1] = mat[1,5]
        mat[5,2] = -(self.K_n * self.L_2**2 * math.cos(self.alpha_2)**2)\
                    + -(self.K_s * self.L_2**2 * math.sin(self.alpha_2)**2)

        # bottom right (3X3)
        mat[3,3] = (self.K_n * math.sin(angle_sum)**2)\
                    + (self.K_s * math.cos(angle_sum)**2)
        mat[3,4] = (-self.K_n * math.sin(angle_sum) * math.cos(angle_sum))\
                    + (self.K_s * math.sin(angle_sum) * math.cos(angle_sum))
        mat[3,5] = (-self.K_n * self.L_1 * math.sin(angle_sum) * math.cos(self.alpha_1))\
                    + (self.K_s * self.L_1 * math.cos(angle_sum) * math.sin(self.alpha_1))
        mat[4,3] = mat[0,1]
        mat[4,4] = (self.K_n * math.cos(angle_sum)**2) \
                    + (self.K_s * math.sin(angle_sum)**2)
        mat[4,5] = (self.K_n * self.L_1 * math.cos(angle_sum) * math.cos(self.alpha_1))\
                    + (self.K_s * self.L_1 * math.sin(angle_sum) * math.sin(self.alpha_1))
        mat[5,3] = mat[0,2]
        mat[5,4] = mat[1,2]
        mat[5,5] = (self.K_n * self.L_1**2 * math.cos(self.alpha_1)**2)\
                    + (self.K_s * self.L_1**2 * math.sin(self.alpha_1)**2)

        return mat

        
print "aem classes imported"
