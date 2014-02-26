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
    def __init__(self,ele_no,b,a,E,G,x,y,node_1=None,node_2=None,node_3=None,node_4=None,r=0,T=0.15):
        self.ele_no = ele_no
        self.edg_len = {1:a,2:b,3:a,4:b,5:a,6:b}
        self.E = E
        self.G = G
        self.T = T
        self.x = x
        self.y = y
        self.r = r
        self.node_1 = node_1
        self.node_2 = node_2
        self.node_3 = node_3
        self.node_4 = node_4
        

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

        print "self.x_1 = ",self.x_1
        print "self.y_1 = ",self.y_1
        print "self.x_2 = ",self.x_2
        print "self.y_2 = ",self.y_2
        print "self.ele_1.x = ",self.ele_1.x
        print "self.ele_1.y = ",self.ele_1.y
        print "self.ele_2.x = ",self.ele_2.x
        print "self.ele_2.y = ",self.ele_2.y        
        
        L = math.sqrt((self.x_2 - self.x_1)**2 + (self.y_2 - self.y_1)**2)
        l = (self.x_2 - self.x_1)/L
        m = (self.y_2 - self.y_1)/L
        
        # Distance from center of ele_1 to first point of spring (d11)
        dx11 = -(self.x_1 - self.ele_1.x)
        dy11 = -(self.y_1 - self.ele_1.y)
        d11 = math.sqrt(dx11**2 + dy11**2)
        
        # Distance from center of ele_2 to second point of spring (d22)
        dx22 = -(self.x_2 - self.ele_2.x)
        dy22 = -(self.y_2 - self.ele_2.y)
        d22 = math.sqrt(dx22**2 + dy22**2)
        
        # Perp. distance from center of ele_1 to the common edge (d1e)
        d1e = self.ele_1.edg_len[self.edg_1 + 1]/2 # using mod 2 +1 to get the base length
        
        # Perp. distance from center of ele_2 to the common edge (d2e)
        d2e = self.ele_2.edg_len[self.edg_2%2 + 1]/2

        K_n = self.K_n
        K_s = self.K_s
        
        print "m = ",m
        print "l = ",l
        print "dx11 = ",dx11
        print "dx22 = ",dx22
        print "dy11 = ",dy11
        print "dy22 = ",dy22
        
        d11 = (m*dx11 - l*dy11)
        d22 = (m*dx22 - l*dy22)
        print d11
        print d22

        L = numpy.array([[l, m, 0, 0, 0, 0],\
                          [-m, l, 0, 0, 0, 0],\
                          [0, 0, 1, 0, 0, 0],\
                          [0, 0, 0, l, m, 0],\
                          [0, 0, 0, -m, l, 0],\
                          [0, 0, 0, 0, 0, 1]])
        #L = numpy.array([[l, m, 0, 0, 0, -m/d1e],\
        #                  [-m, l, 0, 0, 0, l/d1e],\
        #                  [0, 0, 1, 0, 0, -d2e/d1e],\
        #                  [0, 0, 0, l, m, 0],\
        #                  [0, 0, 0, -m, l, 0],\
        #                  [m/d2e, -l/d2e, -d1e/d2e, 0, 0, 1]])

        k = numpy.array([[K_n, 0, -K_n*d11, -K_n, 0, K_n*d22],\
                          [0, K_s, K_s*d1e, 0, -K_s, K_s*d2e],\
                          [-K_n*d11, K_s*d1e, K_n*(d11)**2+K_s*(d1e)**2, K_n*d11, -K_s*d1e, -K_n*d11*d22+K_s*d1e*d2e],\
                          [-K_n, 0, K_n*d11, K_n, 0, -K_n*d22],\
                          [0, -K_s, -K_s*d1e, 0, K_s, -K_s*d2e],\
                          [K_n*d22, K_s*d2e, -K_n*d11*d22+K_s*d1e*d2e, -K_n*d22, -K_s*d2e, K_n*(d22)**2+K_s*(d2e)**2]])
             

        K = numpy.dot(L.T,numpy.dot(k,L))

        return K

class BCspring:
    """BCspring class: spring connecting to the boundary condition
    ele_1 = element one
    edg_1 = edge no. of element one
    BC = 1 or 2 or 3 
        where     1: only normal ristriction on the edge
                  2: only tangential ristriction on the edge
                  3: both normal and tangential ristriction on the edge  
    d = height of spring
    T = thickness of spring
    K_n = Normal spring constant
    K_s = Shear spring constant
    x_1 = x coordinate of first point of spring
    y_1 = y coordinate of first point of spring
    x_2 = x coordinate of center of edge of the BC
    y_2 = y coordinate of center of edge of the BC
    """

    def __init__(self,spr_no, ele_1,edg_1,BC,d,x_1,y_1,x_2,y_2):
        self.spr_no = spr_no
        
        self.ele_1 = ele_1
        self.edg_1 = edg_1

        self.BC = BC

        self.d = d
        T = ele_1.T
        self.T = T

        K_n_1 = 1.0*(ele_1.E * d * T)/(self.ele_1.edg_len[edg_1+1]/2.)
        K_n = K_n_1
        self.K_n = K_n
        
        K_s_1 = 1.0*(ele_1.G * d * T)/(self.ele_1.edg_len[edg_1+1]/2.)
        K_s = K_s_1
        self.K_s = K_s

        self.x_1 = x_1
        self.y_1 = y_1
        
        self.x_2 = x_2
        self.y_2 = y_2

    def stiffness_matrix(self):
        """
        Returns spring stiffness matrix (6X6)
        """

        print "self.x_1 = ",self.x_1
        print "self.y_1 = ",self.y_1
        print "self.x_2 = ",self.x_2
        print "self.y_2 = ",self.y_2
        print "self.ele_1.x = ",self.ele_1.x
        print "self.ele_1.y = ",self.ele_1.y
             
        
        L = math.sqrt((self.x_2 - self.x_1)**2 + (self.y_2 - self.y_1)**2)
        l = (self.x_2 - self.x_1)/L
        m = (self.y_2 - self.y_1)/L
        
        # Distance from center of ele_1 to first point of spring (d11)
        dx11 = -(self.x_1 - self.ele_1.x)
        dy11 = -(self.y_1 - self.ele_1.y)
        d11 = math.sqrt(dx11**2 + dy11**2)
        
        # Perp. distance from center of ele_1 to the common edge (d1e)
        d1e = self.ele_1.edg_len[self.edg_1 + 1]/2 # using mod 2 +1 to get the base length

        K_n = self.K_n
        K_s = self.K_s
        
        print "m = ",m
        print "l = ",l
        print "dx11 = ",dx11
        print "dy11 = ",dy11
        
        d11 = (m*dx11 - l*dy11)
        print d11

        L = numpy.array([[l, m, 0],\
                          [-m, l, 0],\
                          [0, 0, 1]])

        #L = numpy.array([[l, m, 0, 0, 0, -m/d1e],\
        #                  [-m, l, 0, 0, 0, l/d1e],\
        #                  [0, 0, 1, 0, 0, -d2e/d1e],\
        #                  [0, 0, 0, l, m, 0],\
        #                  [0, 0, 0, -m, l, 0],\
        #                  [m/d2e, -l/d2e, -d1e/d2e, 0, 0, 1]])

        k = numpy.array([[K_n, 0, -K_n*d11],\
                          [0, K_s, K_s*d1e],\
                          [-K_n*d11, K_s*d1e, K_n*(d11)**2+K_s*(d1e)**2]])

        K = numpy.dot(L.T,numpy.dot(k,L))

        return K


        
    

print "aem classes imported"