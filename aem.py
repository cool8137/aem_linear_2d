# list of classes used by the Applied Element Method (AEM) program
import math
import numpy as np

def pos_atan(y,x):
    """
    returns invers tangent with a positive angle
    """
    if (x == 0):
        return math.pi/2.
    #elif (y/x > 0):
    #    return math.atan(y/x)
    else:
        return math.atan(y/x) # + math.pi

class AemError():
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

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
    nodes = [node_1, node_2, node_3, node_4]
    r = rotation angle in radian of element
    u = displacements
    [note: when r is zero, edg no. 2 is towards the positive x axis]
    """
    def __init__(self,ele_no,b,a,E,G,x,y,nodes=None,r=0,T=0.2):
        self.ele_no = ele_no
        self.edg_len = {1:a,2:b,3:a,4:b,5:a,6:b}
        self.E = E
        self.G = G
        self.T = T
        self.x = x
        self.y = y
        self.r = r
        self.nodes = nodes
        self.u = np.array([[0,0,0]])
        if nodes == None:
            node_1 = Node(ele_no*10+1, x+a/2*math.cos(r) - b/2*math.sin(r),\
                                y+a/2*math.sin(r) + b/2*math.cos(r))
            node_2 = Node(ele_no*10+1, x+-a/2*math.cos(r) - b/2*math.sin(r),\
                                y-a/2*math.sin(r) + b/2*math.cos(r))
            node_3 = Node(ele_no*10+1, x-a/2*math.cos(r) + b/2*math.sin(r),\
                                y-a/2*math.sin(r) - b/2*math.cos(r))
            node_4 = Node(ele_no*10+1, x+a/2*math.cos(r) + b/2*math.sin(r),\
                                y+a/2*math.sin(r) - b/2*math.cos(r))
            self.nodes = [node_1, node_2, node_3, node_4]
    
    def __str__(self):
        return str(self.__dict__)
        #return ' - '.join(str(getattr(self,key)) for key in self.__dict__ if not key.startswith('_'))

    def addDis(self,u):
        self.u = np.vstack([self.u,u])
        
    def getPlotMat(self,state=0,mag=100):
        u_x = mag * self.u[state,0]
        u_y = mag * self.u[state,1]
        u_w = mag * self.u[state,2]
        
        coord_mat = np.array([[0,0]])
        # loop and add coordinates of each nodes to coord_mat
        # but don't forget to delete first col and also to add node_1 again on the last col
        # note coord_mat = [x1,x2,x3;y1,y2,y3]'
        for node in self.nodes:
            coord = np.array([[node.x,node.y]]).T
            # Transform center to the origin
            coord_from_origin = coord - np.array([[self.x,self.y]]).T
            # Rotate about u_w using rotation matrix w_mat
            w_mat = np.array([[math.cos(u_w), -math.sin(u_w)],[math.sin(u_w), math.cos(u_w)]])
            coord_after_rotation = np.dot(w_mat, coord_from_origin)
            # translate back to center + u_x and u_y
            coord_final = coord_after_rotation + np.array([[self.x + u_x, self.y + u_y ]]).T
            # add the point to coord_mat
            coord_mat = np.vstack([coord_final.T,coord_mat])
        
        # replace last coord [0,0] with the first column
        coord_mat[-1,:] = coord_mat[0,:]
        
        return coord_mat
        
class Elements:
    """ Elements class
    Element objects are listed
    """
    def __init__(self):
        self.list = []
        self.dic = {}
        self.ele2id = {}
        
    def addElement(self, element):
        self.list.append(element)
        new_id = len(self.dic.keys())
        if not element.ele_no in self.dic.keys():
            self.dic[element.ele_no] = element
            self.ele2id[element.ele_no] = new_id
        else:
            raise AemError("Element no. "+str(element.ele_no)+" already exists!!!")



    def addNewElements(self,ele_array, nodes, E, G, T = 0.1):
        """ 
        creates new Element objects and adds it to its list
        input:
            ele_array = [{ele_no} {code} {node_1} {node_2} {node_3} {node_4}]
        """  

        def midpoint(node_1,node_2):
            x1 = node_1.x
            y1 = node_1.y
            x2 = node_2.x
            y2 = node_2.y
            return [(x1 + x2)/2, (y1 + y2)/2]
        
        for i in ele_array:
            ele_no = i[0]
            node_1 = nodes.dic[i[2]]
            node_2 = nodes.dic[i[3]]
            node_3 = nodes.dic[i[4]]
            node_4 = nodes.dic[i[5]]
                        
            a = node_1.distanceTo(node_2)
            b = node_2.distanceTo(node_3)
            center = midpoint(node_1, node_3)
            x = center[0]
            y = center[1]
            
            r = pos_atan(node_1.y-node_2.y,node_1.x-node_2.x)
            self.addElement(Element(ele_no, b, a, E, G, x, y,[node_1,node_2,node_3,node_4],r=r,T=T))

    def global_stiffness_matrix(self,springs=None,bcsprings=None):
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
        N = len(self.list)
        #global_mat = lil_matrix((3*N,3*N))
        global_mat = np.zeros((3*N,3*N))
        
        if springs != None:
            for i in range(len(springs.list)):
                spr = springs.list[i]
        
                spr_mat = spr.stiffness_matrix()
                #print spr_mat.shape
                
                ele_i = self.ele2id[spr.ele_1.ele_no]
                ele_j = self.ele2id[spr.ele_2.ele_no]
                
                # Location vector
                lv = np.array([[(ele_i*3), (ele_i*3+1), (ele_i*3+2),\
                                (ele_j*3), (ele_j*3+1), (ele_j*3+2)]])
                
                global_mat[lv.T,lv] = global_mat[lv.T,lv] + spr_mat
        
                # Note while extracting array from array
                # A[[[r1],[r2],[r3]],[c1,c2,c3,c4]]
        if bcsprings != None:
            for i in range(len(bcsprings.list)):
                BC = bcsprings.list[i]
        
                BC_mat = BC.stiffness_matrix()
                #print spr_mat.shape
                
                ele_i = self.ele2id[BC.ele_1.ele_no]
                
                # Location vector
                lv = np.array([[(ele_i*3), (ele_i*3+1), (ele_i*3+2)]])
                
                global_mat[lv.T,lv] = global_mat[lv.T,lv] + BC_mat
        
                # Note while extracting array from array
                # A[[[r1],[r2],[r3]],[c1,c2,c3,c4]]
        
        return global_mat
        
    def load_vector(self, loads):
        """ 
        returns load_vector
        input: loads object
        """
        load_vec = np.zeros((3*len(self.list),1))

        for load in loads:
            ele_id = self.ele2id[load.ele.ele_no]
            load_vec[3*ele_id,0] += load.F_x
            load_vec[3*ele_id + 1, 0] += load.F_y
            load_vec[3*ele_id + 2, 0] += load.M

        return load_vec


class Node:
    """ 
    Node class    
    """
    
    def __init__(self, node_no, x, y):
        self.node_no = node_no
        self.x = x
        self.y = y

    def getCoord(self):
        return [self.x,self.y]
        
    def distanceTo(self,node_2):
        return math.sqrt((node_2.x - self.x)**2 + (node_2.y - self.y)**2)

    

class Nodes:
    """ Nodes class
    It consists of Node objects listed
    """
    def __init__(self):
        self.list = []
        self.dic = {}
        
    def addNode(self, node):
        self.list.append(node)
        self.dic[node.node_no] = node

    def addNewNodes(self, nodes_array):
        """
        Creates new Node objects and adds it to the Nodes objects list
        input: nodes_array = [{node_no} {x} {y}]
        """
        for i in nodes_array:
            node_no = i[0]
            x = i[1]
            y = i[2]
            self.addNode(Node(node_no,x,y))

    
    
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

        #print "self.x_1 = ",self.x_1
        #print "self.y_1 = ",self.y_1
        #print "self.x_2 = ",self.x_2
        #print "self.y_2 = ",self.y_2
        #print "self.ele_1.x = ",self.ele_1.x
        #print "self.ele_1.y = ",self.ele_1.y
        #print "self.ele_2.x = ",self.ele_2.x
        #print "self.ele_2.y = ",self.ele_2.y        
        
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
        d1e = self.ele_1.edg_len[self.edg_1 + 1]/2. # using mod 2 +1 to get the base length
        
        # Perp. distance from center of ele_2 to the common edge (d2e)
        d2e = self.ele_2.edg_len[self.edg_2%2 + 1]/2.

        K_n = self.K_n
        K_s = self.K_s
        
        #print "m = ",m
        #print "l = ",l
        #print "dx11 = ",dx11
        #print "dx22 = ",dx22
        #print "dy11 = ",dy11
        #print "dy22 = ",dy22
        
        d11 = (m*dx11 - l*dy11)
        d22 = (m*dx22 - l*dy22)
        #print d11
        #print d22

        L = np.array([[l, m, 0, 0, 0, 0],\
                          [-m, l, 0, 0, 0, 0],\
                          [0, 0, 1, 0, 0, 0],\
                          [0, 0, 0, l, m, 0],\
                          [0, 0, 0, -m, l, 0],\
                          [0, 0, 0, 0, 0, 1]])

        k = np.array([[K_n, 0, -K_n*d11, -K_n, 0, K_n*d22],\
                          [0, K_s, K_s*d1e, 0, -K_s, K_s*d2e],\
                          [-K_n*d11, K_s*d1e, K_n*(d11)**2+K_s*(d1e)**2, K_n*d11, -K_s*d1e, -K_n*d11*d22+K_s*d1e*d2e],\
                          [-K_n, 0, K_n*d11, K_n, 0, -K_n*d22],\
                          [0, -K_s, -K_s*d1e, 0, K_s, -K_s*d2e],\
                          [K_n*d22, K_s*d2e, -K_n*d11*d22+K_s*d1e*d2e, -K_n*d22, -K_s*d2e, K_n*(d22)**2+K_s*(d2e)**2]])
             

        K = np.dot(L.T,np.dot(k,L))
        #print K
        return K


class Springs:
    """ Springss class
    It consists of Spring objects listed
    """
    def __init__(self):
        self.list = []
        self.dic = {}
        
    def addSpring(self, spring):
        self.list.append(spring)
        if not spring.spr_no in self.dic.keys():
            self.dic[spring.spr_no] = spring
        else:
            raise AemError("Spring no. "+str(spring.spr_no)+" already exists!!!")

    
    def addNewSprings_auto(self, ele_com_mat):
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
        spr_no = len(self.dic.keys())
    
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
                beta =   {2:ele_1.r+ math.pi/2,\
                        3:(ele_1.r+math.pi),\
                        4:(ele_1.r+3.*math.pi/2),\
                        1:(ele_1.r)}[edg_1] 
            
                #Calc coordinate for first point of spring
                x_1 = ele_1.x + pos_from_center_1 * math.cos(beta)
                y_1 = ele_1.y + pos_from_center_1 * math.sin(beta)
    
                # Calc dis from first point of spring to second (ie length of spring)
                spr_len = ele_1.edg_len[edg_1%2+1]/2. + ele_2.edg_len[edg_2%2+1]/2.
    
                # Calc coordinate of second point of spring
                x_2 = x_1 + spr_len * math.cos(beta-math.pi/2)
                y_2 = y_1 + spr_len * math.sin(beta-math.pi/2)
                
                # Creat new spring and add it to the list of springs (ie the Springs object)
                self.addSpring(Spring(spr_no,ele_1,ele_2,edg_1,edg_2,d,x_1,y_1,x_2,y_2))
                spr_no = spr_no + 1
    

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

        K_n_1 = 1.0*(ele_1.E * d * T)/(self.ele_1.edg_len[edg_1%2+1]/2.)
        K_n = K_n_1
        self.K_n = K_n
        
        K_s_1 = 1.0*(ele_1.G * d * T)/(self.ele_1.edg_len[edg_1%2+1]/2.)
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
        #print "self.x_1 = ",self.x_1
        #print "self.y_1 = ",self.y_1
        #print "self.x_2 = ",self.x_2
        #print "self.y_2 = ",self.y_2
        #print "self.ele_1.x = ",self.ele_1.x
        #print "self.ele_1.y = ",self.ele_1.y     
        
        L = math.sqrt((self.x_2 - self.x_1)**2 + (self.y_2 - self.y_1)**2)
        l = (self.x_2 - self.x_1)/L
        m = (self.y_2 - self.y_1)/L
        
        # Distance from center of ele_1 to first point of spring (d11)
        dx11 = -(self.x_1 - self.ele_1.x)
        dy11 = -(self.y_1 - self.ele_1.y)
        #d11 = math.sqrt(dx11**2 + dy11**2)
        
        
        # Perp. distance from center of ele_1 to the common edge (d1e)
        d1e = self.ele_1.edg_len[self.edg_1 + 1]/2. # using mod 2 +1 to get the base length
        
        

        K_n = self.K_n
        K_s = self.K_s
        
        #print "m = ",m
        #print "l = ",l
        #print "dx11 = ",dx11
        #
        #print "dy11 = ",dy11

        
        d11 = (m*dx11 - l*dy11)
        #print d11


        L = np.array([[l, m, 0],\
                          [-m, l, 0],\
                          [0, 0, 1]])
        #L = np.array([[l, m, 0, 0, 0, -m/d1e],\
        #                  [-m, l, 0, 0, 0, l/d1e],\
        #                  [0, 0, 1, 0, 0, -d2e/d1e],\
        #                  [0, 0, 0, l, m, 0],\
        #                  [0, 0, 0, -m, l, 0],\
        #                  [m/d2e, -l/d2e, -d1e/d2e, 0, 0, 1]])

        k = np.array([[K_n, 0, -K_n*d11],\
                          [0, K_s, K_s*d1e],\
                          [-K_n*d11, K_s*d1e, K_n*(d11)**2+K_s*(d1e)**2]])
             

        K = np.dot(L.T,np.dot(k,L))
        #print K
        return K


class BCsprings:
    """ Springss class
    It consists of BCspring objects listed
    """
    def __init__(self, springs=[]):
        self.list = []
        self.dic = {}

        
    def addSpring(self, spring):
        self.list.append(spring)
        if not spring.spr_no in self.dic.keys():
            self.dic[spring.spr_no] = spring
        else:
            raise AemError("BCspring no. "+str(spring.spr_no)+" already exists!!!")

    
    def addNewSprings_auto(self,ele_com_mat):
        """
        Make auto generated BCsprings when 
        Input: element combination matrix
        Output: spring combination matrix
    
        ele_comb_mat = [{ele_1} {edg_1} {BC} {n}]
        {ele_1} = col. vector of first element no.
        {edg_1} = col. vector of first element's edge no.
        {BC} = BC type
        {n} = no. of springs per edge
        """
    
        # first row is defined so as to define the shape of matrix.
        # don't forget to delete the first row of spr_com_mat
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
                beta =   {2:(ele_1.r + math.pi/2.),\
                        3:(ele_1.r+math.pi),\
                        4:(ele_1.r+3*math.pi/2),\
                        1:(ele_1.r)}[edg_1]
            
                #Calc coordinate for first point of spring
                x_1 = ele_1.x + pos_from_center_1 * math.cos(beta)
                y_1 = ele_1.y + pos_from_center_1 * math.sin(beta)
    
                # Calc dis from first point of spring to edge (ie length of spring)
                spr_len = ele_1.edg_len[edg_1%2+1]/2.
    
                # Calc coordinate of second point of spring (at the center of edge)
                x_2 = x_1 + spr_len * math.cos(beta-math.pi/2)
                y_2 = y_1 + spr_len * math.sin(beta-math.pi/2)
                #print y_1
                #print x_1
                
                # Stack all info into the spr_com_mat matrix
                self.addSpring(BCspring(spr_no,ele_1,edg_1,BC,d,x_1,y_1,x_2,y_2))
                spr_no = spr_no + 1
                
class Load:
    """
    Load objects
    """
    
    def __init__(self, ele, edg_no=0, F=0, M=0):
        # input load
        # ele = element object
        # edg_no = edg on which the load will be applied
        # F is positive than the load causes positive (anti-clockwise) rotation
        self.ele = ele
        if edg_no == 1:
            self.F_x = F
            self.F_y = 0
            self.M = M + F * ele.edg_len[edg_no%2+1]/2.
        elif edg_no == 2:
            self.F_x = 0
            self.F_y = F
            self.M = M + F * ele.edg_len[edg_no%2+1]/2.
        elif edg_no == 3:
            self.F_x = -F
            self.F_y = 0
            self.M = M + F * ele.edg_len[edg_no%2+1]/2.
        elif edg_no == 4:
            self.F_x = 0
            self.F_y = -F
            self.M = M + F * ele.edg_len[edg_no%2+1]/2.
        else:
            self.F_x = 0
            self.F_y = 0
            self.M = M

        



print "aem classes imported"