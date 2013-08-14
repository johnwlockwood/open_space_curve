# -*- coding: utf-8 -*-
#transformations generator
#builds a transformation table
#zero order is 1 point
#first order is 4 points
#initial state is bottom-left,top-left, top-right, bottom-right
#[[0,0],[0,1],[1,1],[1,0]]

#Degrees	Radians
#	90°	pi/2
#	60°	pi/3
#	45°	pi/4
#	30°	pi/6
import numpy as np
import math
import transformations as tr
import logging
import copy
from time import clock, time

#firstState = np.matrix([[0,0],[0,1],[1,1],[1,0]])
#plus90 = np.matrix([[0,-1],[1,0]])
#firstState2 = np.matrix([[0,0,1],[0,1,1],[1,1,1],[1,0,1]])
#plus902 = np.matrix([[0,-1,0],[1,0,0],[0,0,1]])
#plus90roundcenter = np.matrix([[0,-1,0],[1,0,0],[0,0,1]])
#transmatrix = np.matrix([[1,0,1],[0,1,1],[0,0,1]])

#T = tr.translation_matrix([-0.5,-0.5,0])
#T.dot([1,1,0,1])


#matrix for rotating a normalized 3D point around the center on a z axis
#you can use this to rotate a normalized 2D point around it's center
# by adding a 0 for the z access, plus another 1
# [1,1] => [1,1,0,1]
#R = tr.rotation_matrix(-math.pi/2,[0,0,1],[0.5,0.5,0])
# then to tranform it, you get the dot product of the two
# R.dot(numpy.array([1,1,0,1])
"""
rotate an array of points:
>>> state0 = np.matrix([[0, 0, 0, 1],
...        [0, 1, 0, 1],
...        [1, 1, 0, 1],
...        [1, 0, 0, 1]])

first transpose the point matrix for alignment with the transform matrix
then get the dot product, then transpose the resulting point matrix back

>>> R.dot(state0.getT()).getT()

then the first two values of each array of the resulting matrix will
be the x and y of the 2D point

OR one point at a time:

>>> for point in state0.tolist():
...  print R.dot(point)[0:2]

"""



#R2 = tr.rotation_matrix(math.pi/2,[0,0,1],[0,0,0])
"""
if you want to transform a point by a matrix you do matrix.dot(point)

"""
MAX_ORDER = 29



STATE0M = np.matrix([[0, 0, 0, 1],
    [0, 1, 0, 1],
    [1, 1, 0, 1],
    [1, 0, 0, 1]])

def faceRepositionMatrixes():
    R40CCW = tr.rotation_matrix((math.pi/2)*0.4444444,[1,0,0],[0.5,0.5,0.5])
    R40CW = tr.rotation_matrix(-(math.pi/2)*0.4444444,[1,0,0],[0.5,0.5,0.5])
    return R40CW, R40CCW

def get_standard_matrixes():
    #90 degree clockwise around center of normalized square
    R90CW = tr.rotation_matrix(-math.pi/2,[0,0,1],[0.5,0.5,0.5])
    #R90CCW = tr.rotation_matrix(math.pi/2,[0,0,1],[0.5,0.5,0.5])

    #horizontal axis flip
    FH = tr.reflection_matrix([0.5,0.5,0],[0,1,0])

    #vertical axis flip
    FV = tr.reflection_matrix([0.5,0.5,0],[1,0,0])

    R90CW_FH = tr.concatenate_matrices(FH,R90CW)

    R90CW_FV = tr.concatenate_matrices(FV,R90CW)

    R180CW = tr.rotation_matrix(-math.pi,[0,0,1],[0.5,0.5,0])
    #R180CCW = tr.rotation_matrix(math.pi,[0,0,1],[0.5,0.5,0])
    #R270CCW = tr.rotation_matrix(2*math.pi,[0,0,1],[0.5,0.5,0])

    return R90CW,FH,FV,R90CW_FH,R90CW_FV,R180CW

def create_curve_states_2D():
    R90CW,FH,FV,R90CW_FH,R90CW_FV,R180CW = get_standard_matrixes()
    state0m = np.matrix([[0, 0, 0, 1],
        [0, 1, 0, 1],
        [1, 1, 0, 1],
        [1, 0, 0, 1]])


    state0 = np.round_(state0m.getT()[:2:].getT().astype(int)).astype(int)
    state1 = np.round_(R90CW_FH.dot(state0m.getT())[:2:].getT()).astype(int)
    state2 = np.round_(R90CW_FV.dot(state0m.getT())[:2:].getT()).astype(int)
    state3 = np.round_(R180CW.dot(state0m.getT())[:2:].getT()).astype(int)

    return np.array([state0,state1, \
    state2, state3])


#CURVESTATES2D = create_curve_states_2D()
#caching results to reduce module load time
CURVESTATES2D = np.array([[[0, 0],
        [0, 1],
        [1, 1],
        [1, 0]],

       [[0, 0],
        [1, 0],
        [1, 1],
        [0, 1]],

       [[1, 1],
        [0, 1],
        [0, 0],
        [1, 0]],

       [[1, 1],
        [1, 0],
        [0, 0],
        [0, 1]]])






def pad_state2d(curve_state):
    return np.matrix(np.insert(curve_state,[2,2], [0,1], axis=1),dtype=int)


def permutate_curve_state_matrix(curve_state):
    R90CW,FH,FV,R90CW_FH,R90CW_FV,R180CW = get_standard_matrixes()
    PERMUTATION_TRANSFORM_PATTERN = [
        R90CW_FH,
        None,
        None,
        R90CW_FV,
    ]
    state0m = pad_state2d(curve_state)
    permutated_states = []
    for trans in PERMUTATION_TRANSFORM_PATTERN:
        if not trans is None:
            permutated_states.append(
                np.round_(trans.dot(state0m.getT())[:2:].getT()).astype(int))
        else:
            permutated_states.append(curve_state)

    return np.array(permutated_states)


def index_of_state(curve_state):
    for index,state in enumerate(CURVESTATES2D):
        if np.allclose(curve_state,state):
            return index
    return -1


def create_permutation_keys():
    d2permutation_keys = []

    for curve_state in CURVESTATES2D:
        state_permed =  permutate_curve_state_matrix(curve_state)
        permed_keys = []
        for statep in state_permed:
            permed_keys.append(index_of_state(statep))

        d2permutation_keys.append(permed_keys)

    return d2permutation_keys


#PERMUTATION_KEYS_2D = create_permutation_keys()
#caching results to reduce module load time
PERMUTATION_KEYS_2D =  [[1, 0, 0, 2], [0, 1, 1, 3], [3, 2, 2, 0], [2, 3, 3, 1]]


twoFourArray = np.tile(np.array([2, 2]),(4,1))



def point_2D_from_geopt(geopt):
    return np.array([geopt.lon,geopt.lat])


def get_next_order_fold_index_and_point_2D(indexing_point,
    curve_point,state_key,order):
    increment = distance_increment_at_order_2D(order+1)
    next_states = CURVESTATES2D[state_key]
    next_bounding_boxes = []

    for fold_index, state in enumerate(next_states):
        if state[0] == 1:
            d1 = curve_point[0]+increment
        else:
            d1 = curve_point[0]-increment
        if state[1] == 1:
            d2 = curve_point[1]+increment
        else:
            d2 = curve_point[1]-increment

        next_curve_point = np.array([d1,d2])
        bbox = BoundingBox_2D.square_from_center_point_and_size(
            next_curve_point,increment)
        point_in = bbox.point_in(indexing_point)
        if point_in:
            return fold_index, next_curve_point


def get_next_order_point_2D(the_fold_index,curve_point,state_key,order):
    increment = distance_increment_at_order_2D(order+1)
    next_states = CURVESTATES2D[state_key]

    state = next_states[the_fold_index]
    if state[0] == 1:
        d1 = curve_point[0]+increment
    else:
        d1 = curve_point[0]-increment
    if state[1] == 1:
        d2 = curve_point[1]+increment
    else:
        d2 = curve_point[1]-increment

    return np.array([d1,d2])


def index_point_at_orders_2D(faceCoords, order):
    if faceCoords[0] < 0 or faceCoords[0] > 1 or faceCoords[1] < 0 or faceCoords[1] > 1:
        raise ValueError("please normalize your values")

    point = np.array([0.5,0.5])
    state_key = 0
    fold_index_at_order = []
    next_order = 0

    while next_order<=order:
        fold_index, point = get_next_order_fold_index_and_point_2D( \
            faceCoords, point,state_key,next_order)
        fold_index_at_order.append(fold_index)
        state_key = PERMUTATION_KEYS_2D[state_key][fold_index]
        next_order += 1

    return fold_index_at_order


def index_point_2D(faceCoords, order):
    if faceCoords[0] < 0 or faceCoords[0] > 1 or faceCoords[1] < 0 or faceCoords[1] > 1:
        raise ValueError("please normalize your values")

    cfTime = 0.0
    cfCount = 0

    point = np.array([0.5,0.5])
    state_key = 0
    index = 0
    next_order = 0
    curve_key = ""
    statesKey = ""
    fold_index = 0
    while next_order<=order:

        startCFtime = time()
        fold_index, point = get_next_order_fold_index_and_point_2D( \
            faceCoords, point,state_key,next_order)
        cfTime += time() - startCFtime
        cfCount += 1

        index = index*CHILDREN_COUNT_2D+fold_index
        statesKey += str(state_key)
        state_key = PERMUTATION_KEYS_2D[state_key][fold_index]
        curve_key += str(fold_index)
        next_order += 1


    curve_point = HilbertSpaceCurvePoint(index,point,order,fold_index=fold_index, \
        state_key=state_key,states_key=statesKey,curve_key=curve_key)

    timings = (cfTime,cfCount)

    return curve_point, timings


def point_from_fold_indexes_2D(fold_indexes):
    cfTime = 0.0
    cfCount = 0

    point = np.array([0.5,0.5])
    state_key = 0
    index = 0
    next_order = 0
    curve_key = ""
    statesKey = ""
    fold_index = 0
    for fold_index in fold_indexes:

        startCFtime = time()
        point = get_next_order_point_2D(
            fold_index, point,state_key,next_order)
        cfTime += time() - startCFtime
        cfCount += 1

        index = index*CHILDREN_COUNT_2D+fold_index
        statesKey += str(state_key)
        state_key = PERMUTATION_KEYS_2D[state_key][fold_index]
        curve_key += str(fold_index)
        next_order += 1

    #order = next_order-1
    #curve_point = HilbertSpaceCurvePoint(index,point,order,fold_index=fold_index,
    #    state_key=state_key,states_key=statesKey,curve_key=curve_key)

    timings = (cfTime,cfCount)

    return point, timings


def get_capkey(geopoint,order=MAX_ORDER):
    cfc = sphere_to_cube_face_2D(geopoint)
    #logging.debug("face %s"%cfc.face)
    #logging.debug("face %s"%cfc.position)
    return cfc.encode_capkey_to_order(order)

def distance_increment_at_order_2D(order):
    """ considering 0 as first order.
    so we add 2 because the first fold is
    a quarter of 1
    and 1.0/(1<<2) == 0.25
    """
    return 1.0/(1<<(order + 1))

#INCREMENT_DISTANCES = [distance_increment_at_order_2D(order) for order in range(0,34)]

def get_dimensional_position_of_children_2D(point,pemutation_state_key,order):
    increment = distance_increment_at_order_2D(order+1)#INCREMENT_DISTANCES[order+1]
    state = CURVESTATES2D[pemutation_state_key]
    positions = []

    for value in state:
        if value[0] == 1:
            d1 = point[0]+increment
        else:
            d1 = point[0]-increment
        if value[1] == 1:
            d2 = point[1]+increment
        else:
            d2 = point[1]-increment
        positions.append(np.array([d1,d2]))
    return np.array(positions)


def get_point_2D_at_index(index,order):
    index = int(index)
    point = np.array([0.5,0.5])
    dimensional_positions = get_dimensional_positions_at_orders_2D(index,order)
    key = 0
    while key <= order:
        increment = distance_increment_at_order_2D(key+1)#INCREMENT_DISTANCES[key+1]

        if dimensional_positions[key][0] == 1:
            d1 = point[0]+increment
        else:
            d1 = point[0]-increment

        if dimensional_positions[key][1] == 1:
            d2 = point[1]+increment
        else:
            d2 = point[1]-increment

        point = np.array([d1,d2])
        key += 1

    return point


def clampToGeoPt(lat,lon):
    try:
        from google.appengine.api.datastore_types import GeoPt
    except:
        from geopt_alt import GeoPt
    if lat > 90.0:
        lat = 90.0
    elif lat < -90.0:
        lat = -90.0

    if lon > 180.0:
        lon = 180.0
    elif lon < -180.0:
        lon = -180.0

    return GeoPt(lat,lon)

class BoundingBox_2D(object):
    def __init__(self,pointLB,pointRT):
        self.pointLB = pointLB
        self.pointRT = pointRT

    def point_in(self,point):
        return (self.pointLB[0] <= point[0] and point[0] <= self.pointRT[0]) \
        and (self.pointLB[1] <= point[1]and point[1] <= self.pointRT[1])

    def corner_points(self,res=1):
        if res > 1:
            pass
        return [self.pointLB,np.array([self.pointLB[0],self.pointRT[1] ]),
        self.pointRT,np.array([self.pointRT[0] ,self.pointLB[1]])]

    def corner_geopts(self):

        return [clampToGeoPt(self.pointLB[1],self.pointLB[0]),
        clampToGeoPt(self.pointRT[1],self.pointLB[0]),
        clampToGeoPt(self.pointRT[1],self.pointRT[0]),
        clampToGeoPt(self.pointLB[1],self.pointRT[0])]

    def get_center_as_geopt(self):
        latlon = self.get_center()
        return clampToGeoPt(latlon[0],latlon[1])

    def get_center(self):
        height = self.pointRT[1] - self.pointLB[1]
        width = self.pointRT[0] - self.pointLB[0]
        half_height = height * 0.5
        half_width = width * 0.5
        return np.array([self.pointLB[0]+half_height,self.pointLB[1]+half_width])


    def diagonal_size(self):
        return math.sqrt(distance_squared_of_2_2d_point(self.pointLB,
            self.pointRT))

    def is_inverted(self):
        if self.pointLB[0] > self.pointRT[0]:
            return True
        if self.pointLB[1] > self.pointRT[1]:
            return True


    def expand_corners_by(self,amount):
        amount_array = np.array([amount,amount])
        lb = np.subtract(self.pointLB,amount_array)
        rt = np.add(self.pointRT,amount_array)
        return lb,rt

    def contract_corners_by(self,amount):
        amount_array = np.array([amount,amount])
        lb = np.add(self.pointLB,amount_array)
        rt = np.subtract(self.pointRT,amount_array)
        return lb,rt

    @classmethod
    def square_from_center_point_and_size(cls,center_point,halfsize):
        #halfsize = size*0.5
        halfsize_array = np.array([halfsize,halfsize])
        pointLB = np.subtract(center_point,halfsize_array)
        pointRT = np.add(center_point,halfsize_array)
        return cls(pointLB,pointRT)

    @classmethod
    def from_points(cls,points):
        top = None
        bottom = None
        left = None
        right = None
        for point in points:
            if top is None or point[1] > top:
                top = point[1]
            if bottom is None or point[1] < bottom:
                bottom = point[1]
            if right is None or point[0] > right:
                right = point[0]
            if left is None or point[0] < left:
                left = point[0]

        return cls(np.array([left,bottom]),np.array([right,top]))




def distance_squared_of_2_2d_point(pointA,pointB):
    return math.pow((pointB[0]- pointA[0]),2)+math.pow((pointB[1]- pointA[1]),2)

def get_next_curve_key(curve_key):
    """
    expects a string like:
    21032
    expressing level and fold_index
>>> fold_index_key = "21033"
>>> next_fold_index_key = get_next_curve_key(fold_index_key)
>>> next_fold_index_key
'21100'

    """
    format = "%%0%dd"%len(curve_key)
    return format%int(np.base_repr(int(curve_key,4)+1,4))

def capkey_ranges_from_curve_point_ranges(face,curve_point_ranges):
    capkey_ranges = []
    face_prefix = "%s"%face
    for start,end in curve_point_ranges:
        capkey_ranges.append((face_prefix+start.get_curve_key(),
            face_prefix+get_next_curve_key(end.get_curve_key())))
    return capkey_ranges

class HilbertSpaceCurvePoint(object):
    def __init__(self,index,point,order,fold_index=None,state_key=None,
        curve_key=None,states_key=None):
        self.index = index
        self.point = point
        self.order = order
        self.curve_key = curve_key
        self.states_key = states_key
        self.fold_index = fold_index
        self.state_key = state_key

    def __str__(self):
        return "(%s, %s, %s, %s, %s)"%(self.index,self.point,
            self.order,self.curve_key,self.states_key)
    def __unicode__(self):
        return u'(%s, %s, %s, %s, %s)'%(self.index,self.point,
            self.order,self.curve_key,self.states_key)
    def __repr__(self):
        return u'openhilbert.HilbertSpaceCurvePoint(%s, %s, %s,curve_key=%s,states_key=%s)'%(self.index,self.point,self.order,self.curve_key,self.states_key)


    def get_curve_key(self):
        return self.curve_key

    def get_curve_key_from_index(self,children_count=4):
        format_string = "%0."+str(self.order + 1)+"d"
        return format_string%int(np.base_repr(self.index,children_count))

    def index_from_curve_key(self,children_count=4):
        key = self.get_curve_key()
        return long(key,children_count)

    def get_parent_curve_point(self):
        order = self.order
        if order > 0:
            order -= 1
            parent = copy.deepcopy(self)
            parent.order = order
            parent.fold_index = self.curve_key[-1]
            parent.state_key = self.states_key[-1]
            parent.index = parent.index_from_curve_key(4)
            parent.states_key = self.states_key[:-1]
            parent.curve_key = self.curve_key[:-1]
            return parent
        else:
            return None



def normalize_cube_position_value(value):
    return (value+1)/2

def denormalize_to_cube_position_value(value):
    return value*2-1

def get_surrounding_capkeys(geopt,order=MAX_ORDER):
    cfc = sphere_to_cube_face_2D(geopt)
    cfcs = cfc.surrounding_face_positions_at_order(order)
    #for cfc in cfcs:
    #    logging.debug(cfc.position)
    capkey_strings = [cfc.encode_capkey_to_order(order=order) for cfc in cfcs]
    return capkey_strings

faceCoordNormShiftP1 = np.array([1,1])
faceCoordNormShiftP2 = np.array([2,2])

faceCoordNormShiftAltP2 = np.array([0.5,0.5])

dim1 = 0
dim2 = 1
dims = [0,1]


class CubeFaceAndCoords(object):
    """
    """

    def __init__(self,face,faceCoords):
        if face > 5: face = 5
        if face < 0: face = 0
        self.face = face

        for index,value in enumerate(faceCoords):
            if value >= 1.0:
                faceCoords[index] = 0.999999999999999999999999999999
            elif value < 0.0:
                faceCoords[index] = 0.0

        self.position = np.array([faceCoords[0],faceCoords[1]],dtype=float)
        self.curve_points = []


    def get_point(self):
        return self.position

    def compute_curve_points(self,order):
        self.curve_points = index_point_at_orders_2D(self.position,order)
        return self.curve_points


    def compute_capkey_parts(self,order=MAX_ORDER):
        capkey_parts = [self.face]
        capkey_parts.extend(index_point_at_orders_2D(self.position,order))
        return [str(key) for key in capkey_parts]


    def compute_point_at_order(self,order):
        if self.curve_points and len(self.curve_points)>order:
            curve_points = self.curve_points
        else:
            curve_points = self.compute_curve_points(order)

        point = curve_points[order].point

        return copy.deepcopy(point)

    def encode_capkey_to_order(self,order=MAX_ORDER):
        """
        '31331320020312'
        """
        capkey_parts = self.compute_capkey_parts(order=order)


        return "".join(capkey_parts)


    def surrounding_face_positions_at_order(self,order=MAX_ORDER,
        scale_factor=1.0,res=1):
        increment = distance_increment_at_order_2D(order)*scale_factor #INCREMENT_DISTANCES[order]*scale_factor
        curve_point, timings = index_point_2D(self.position,order)

        bbox = BoundingBox_2D.square_from_center_point_and_size(
            curve_point.point,increment)

        corner_points = bbox.corner_points()
        #logging.debug(corner_points)
        cfcs = []
        for point in corner_points:
            cfcs.append(self.__class__(self.face,point))

        return cfcs




    @classmethod
    def center_CubeFaceAndCords_at_order(cfc,order=MAX_ORDER):
        curve_point, timings = index_point_2D(cfc.position,order)
        new_cfc = copy.deepcopy(cfc)
        new_cfc.position = curve_point.point
        return new_cfc

    @classmethod
    def from_capkey(cls,capkey):
        capkey_parts = [int(value) for value in capkey]
        face = int(capkey_parts[0])
        fold_indexes = capkey_parts[1:]
        position, timings = point_from_fold_indexes_2D(fold_indexes)

        return cls(face,position)



    def center_at_order_to_sphere_geopt(self,order):
        """
        """
        curve_point, timings = index_point_2D(self.position,order)
        position = self.to_position3(alt_position = curve_point.point)
        position = cubeCoords_to_sphereCoords(position)
        return cartesian_to_polar_coords_of_sphere(spherePositionFromFaceFilter(position))


    def boundary_face_positions_at_order(self,order=MAX_ORDER,
        scale_factor=1.0,res=1):
        increment = distance_increment_at_order_2D(order+1)*scale_factor #INCREMENT_DISTANCES[order+1]*scale_factor
        curve_point, timings = index_point_2D(self.position,order)
        bbox = BoundingBox_2D.square_from_center_point_and_size(
            curve_point.point,increment)

        corner_points = bbox.corner_points(res=res)
        cfcs = []
        for point in corner_points:
            cfcs.append(self.__class__(self.face,point))

        return cfcs


    def __repr__(self):
        return u'geospatial2d.CubeFaceAndCoords(%s,%s,%s)'%(self.face,self.position[0],self.position[1])





    def change_position(self,dim,sign,amount,from_position = None):
        self.curve_points = None #they are invalid with a position change
        ## so they will be recomputed the next time they are requested
        if from_position is None:
            from_position = self.position


        new_dimensional_position = from_position[dim]+sign*amount
        if 1.0 >= new_dimensional_position >= 0.0:
            self.position[dim]=new_dimensional_position
            self.position[(dim+1)%2] =from_position[(dim+1)%2]
            return self


        if new_dimensional_position <= 1.0:
            self.position[dim]=0.999999999999999
            self.position[(dim+1)%2] =from_position[(dim+1)%2]

        elif new_dimensional_position >= 0.0:
            self.position[dim]=0.000000000000001
            self.position[(dim+1)%2] =from_position[(dim+1)%2]

        return self


    def to_position3(self,alt_position=None):
        #shift from curvespace to cube face

        #logging.debug("step1 s1 %s  t1 %s"%(self.position[0],self.position[1]))
        if alt_position:
            orig_position = alt_position
        else:
            orig_position = self.position


        #change the coordinates from being centered on 0.5,0.5 to being centered on 0,0
        # muliply by 2 then subtract 1
        position = np.subtract(np.multiply(orig_position,faceCoordNormShiftP2),
        faceCoordNormShiftP1)
        #logging.debug("denormed pos %s  pos %s"%(position[0],position[1]))
        position3 = None
        if self.face == 0:
            position3 = np.array([1.0,position[0],position[1]])
        elif self.face == 1:
            position3 = np.array([-1.0,position[0]*-1,position[1]])
        elif self.face == 2:
            position3 = np.array([position[0]*-1,1.0,position[1]])
        elif self.face == 3:
            position3 = np.array([position[0],-1.0,position[1]])
        elif self.face == 4:
            position3 = np.array([position[1]*-1,position[0],1.0])
        elif self.face == 5:
            position3 = np.array([position[1],position[0],-1.0])

        return position3


    def to_sphere_geopt(self):
        """
        """
        cubeCoords = self.to_position3()
        sphereCoords = cubeCoords_to_sphereCoords(cubeCoords)
        return cartesian_to_polar_coords_of_sphere(spherePositionFromFaceFilter(sphereCoords))


oneArray = np.array([1,1,1])
twoArray = np.array([2,2,2])
threeArray = np.array([3,3,3])
halfArray = np.array([0.5,0.5,0.5])
oneThird = 1.0/3
thirdArray = np.array([oneThird,oneThird,oneThird])


def cubeCoords_to_sphereCoords(pA):

    ps1 = np.power(pA,twoArray)

    #ps2 = np.divide(ps1,twoArray)
    ps2 = np.multiply(ps1,halfArray)

    #ps3 = np.divide(ps1,threeArray)
    ps3 = np.multiply(ps1,thirdArray)

    ps4 = np.subtract(oneArray,ps2)

    #logging.debug("\nps1:%s\nps2%s\nps3:%s\nps4%s"%(ps1,ps2,ps3,ps4))
    sp1x = ps4[1]-ps2[2]+ps1[1]*ps3[2]
    sp1y = ps4[2]-ps2[0]+ps1[2]*ps3[0]
    sp1z = ps4[0]-ps2[1]+ps1[0]*ps3[1]

    sp1 = np.array([sp1x,sp1y,sp1z])

    sp2 = np.sqrt(sp1)

    sphereCoords = np.multiply(pA,sp2)

    return sphereCoords


def cartesian_to_polar_coords_of_sphere(position):
    """
    takes: sphere coordinates
    returns: a GeoPt
    """
    try:
        from google.appengine.api.datastore_types import GeoPt
    except ImportError:
        from geopt_alt import GeoPt
    ps1 = np.power(position,twoArray)
    psum = np.sum(ps1)
    psqrt = math.sqrt(psum)
    lat_radians = math.asin(position[2]/psqrt)

    lon_radians = math.atan2(position[1],position[0])
    geopt = GeoPt(math.degrees(lat_radians),math.degrees(lon_radians))
    return geopt

def polar_to_cartesian_coords_of_sphere(geopt):
    """
    """
    lat_radians = math.radians(geopt.lat) #degree to radians
    lon_radians = math.radians(geopt.lon)
    x = math.cos(lat_radians)*math.cos(lon_radians)
    y = math.cos(lat_radians)*math.sin(lon_radians)
    z = math.sin(lat_radians)
    return np.array([x,y,z])

R40CW, R40CCW = faceRepositionMatrixes()
def spherePositionToFaceFilter(pO):
    # run on position before cubizePoint3
    # transforms the point so that the latitudal center of
    # face 3 falls on 36 deg latitude making most of north america on one face
    return R40CCW.dot(np.array([pO[0],pO[1],pO[2],0.0]))[:3]

def spherePositionFromFaceFilter(pO):
    # run on a position after cubeCoords_to_sphereCoords
    return R40CW.dot(np.array([pO[0],pO[1],pO[2],0.0]))[:3]

def cubizePoint3(position):
    #thanks to:
    #http://stackoverflow.com/questions/2656899/mapping-a-sphere-to-a-cube
    import math
    pA = copy.deepcopy(position)
    pB = copy.deepcopy(pA)

    fP = np.absolute(pB)
    faceValue = None

    #inverseSqrt2 = 0.70710676908493042
    inverseSqrt2 = 0.7071067811865475

    #which of x, y or z has the largest absolute value
    # then is the value of that one positive or negative

    if fP[0] <= fP[1] >= fP[2]:
        a2 = pB[0] * pB[0] * 2.0
        b2 = pB[2] * pB[2] * 2.0
        inner = -a2 + b2 -3
        innersqrt = -math.sqrt((inner * inner) - 12.0 * a2)

        if pB[0] == 0.0 or pB[0] == -0.0:
            pA[0] = 0.0

        else:
            pA[0] = math.sqrt(innersqrt + a2 - b2 + 3.0) * inverseSqrt2

        if pB[2] == 0.0 or pB[2] == -0.0:
            pA[2] = 0.0

        else:
            pA[2] = math.sqrt(innersqrt - a2 + b2 + 3.0) * inverseSqrt2

        if pA[0] > 1.0: pA[0] = 1.0
        if pA[2] > 1.0: pA[2] = 1.0

        if pB[0] < 0: pA[0] = -pA[0]
        if pB[2] < 0: pA[2] = -pA[2]

        if pB[1] > 0:
            #face 2
            pA[1] = 1.0
            face = 2
            faceCoords = np.array([pA[0]*-1,pA[2]])
            s1 = pA[0]*-1
            t1 = pA[2]
            faceValue = pB[1]

        else:
            #face 3
            pA[1] = -1.0
            face = 3
            faceCoords = np.array([pA[0],pA[2]])
            s1 = pA[0]
            t1 = pA[2]
            faceValue = pB[1]

    elif fP[1] <= fP[0] >= fP[2]:
        a2 = pB[1] * pB[1] * 2.0
        b2 = pB[2] * pB[2] * 2.0
        inner = -a2 + b2 -3
        innersqrt = -math.sqrt((inner * inner) - 12.0 * a2)

        if pB[1] == 0.0 or pB[1] == -0.0:
            pA[1] = 0.0

        else:
            pA[1] = math.sqrt(innersqrt + a2 - b2 + 3.0) * inverseSqrt2


        if pB[2] == 0.0 or pB[2] == -0.0:
            pA[2] = 0.0

        else:
            pA[2] = math.sqrt(innersqrt - a2 + b2 + 3.0) * inverseSqrt2


        if pA[1] > 1.0: pA[1] = 1.0
        if pA[2] > 1.0: pA[2] = 1.0

        if pB[1] < 0: pA[1] = -pA[1]
        if pB[2] < 0: pA[2] = -pA[2]

        if pB[0] > 0:
            #face 0
            pA[0] = 1.0
            face = 0
            faceCoords = np.array([pA[1],pA[2]])
            s1 = pA[1]
            t1 = pA[2]
            faceValue = pB[0]
        else:
            #face 1
            pA[0] = -1.0
            face = 1
            faceCoords = np.array([pA[1]*-1,pA[2]])
            s1 = pA[1]*-1
            t1 = pA[2]
            faceValue = pB[0]

    else:
        a2 = pB[0] * pB[0] * 2.0
        b2 = pB[1] * pB[1] * 2.0
        inner = -a2 + b2 -3
        innersqrt = -math.sqrt((inner * inner) - 12.0 * a2)

        if pB[0] == 0.0 or pB[0] == -0.0:
            pA[0] = 0.0

        else:
            pA[0] = math.sqrt(innersqrt + a2 - b2 + 3.0) * inverseSqrt2

        if pB[1] == 0.0 or pB[1] == -0.0:
            pA[1] = 0.0

        else:
            pA[1] = math.sqrt(innersqrt - a2 + b2 + 3.0) * inverseSqrt2

        if pA[0] > 1.0: pA[0] = 1.0
        if pA[1] > 1.0: pA[1] = 1.0

        if pB[0] < 0: pA[0] = -pA[0]
        if pB[1] < 0: pA[1] = -pA[1]

        if pB[2] > 0:
            #face 4
            pA[2] = 1.0
            face = 4
            faceCoords = np.array([pA[1],pA[0]*-1])
            s1 = pA[1]
            t1 = pA[0]*-1
            faceValue = pB[2]

        else:
            #face 5
            pA[2] = -1.0
            face = 5
            faceCoords = np.array([pA[1],pA[0]])
            s1 = pA[1]
            t1 = pA[0]
            faceValue = pB[2]


    #shift to curve face coords
    #faceCoords = np.divide(np.add(faceCoords,faceCoordNormShiftP1), \
    #    faceCoordNormShiftP2)

    #logging.debug("faceCoords %s %s"%(faceCoords[0],faceCoords[1]))
    # this shifts the coordinates from centered on 0,0 to centered on 0.5,0.5
    # add 1 then multiply by 0.5
    faceCoords = np.multiply(np.add(faceCoords,faceCoordNormShiftP1),
        faceCoordNormShiftAltP2)

    #logging.debug("normed faceCoords %s %s"%(faceCoords[0],faceCoords[1]))

    cfc = CubeFaceAndCoords(face,faceCoords)
    cfc.faceValue = faceValue

    return cfc

def sphere_to_cube_face_2D(geopt):
    position = polar_to_cartesian_coords_of_sphere(geopt)
    return cubizePoint3(spherePositionToFaceFilter(position))


def boundingBox_from_points_2D(points):
    left = None
    right = None
    top = None
    bottom = None
    for point in points:
        if left is None or left > point[0]:
            left = point[0]
        if right is None or right < point[0]:
            right = point[0]
        if top is None or top < point[1]:
            top = point[1]
        if bottom is None or bottom > point[1]:
            bottom = point[1]

    return BoundingBox_2D(np.array([left,bottom]),np.array([right,top]))


def get_children_curves_2D(curve_point):
    """
    """
    children = []
    order = curve_point.order+1
    children_points = get_dimensional_position_of_children_2D(
        curve_point.point,
            curve_point.state_key,order)

    for fold_index in range(0,4,1):
        statesKey = curve_point.states_key
        curve_key = curve_point.curve_key

        state_key = PERMUTATION_KEYS_2D[curve_point.state_key][fold_index]
        statesKey += str(state_key)
        index = curve_point.index*CHILDREN_COUNT_2D + fold_index
        curve_key += str(fold_index)

        children.append(HilbertSpaceCurvePoint(index,children_points[fold_index],order,
            fold_index, state_key, states_key=statesKey,curve_key=curve_key))

    return children


##querying
# dimensions = 2
# 1<<dimensions
CHILDREN_COUNT_2D = 4

def get_index_at_previous_order_2D(index):
    return int(index/CHILDREN_COUNT_2D)

def gen_index_at_previous_orders_2D(index,order):
    for theOrder in sorted(range(0,order+1),reverse=True):
        yield index
        index = int(index/CHILDREN_COUNT_2D)

def get_fold_index_2D(index):
    return index%CHILDREN_COUNT_2D

def get_indexes_for_orders_2D(index,order):
    indexes = {order:index}
    while order > 0:
        order -= 1
        index = get_index_at_previous_order_2D(index)
        indexes[order]=index
    return indexes

def gen_indexes_for_orders_2D(index,order):
    return (order_index for order_index in gen_index_at_previous_orders_2D(index,order))

def get_fold_indexes_at_orders_2D(index,order):
    return dict(zip(sorted(range(0,order+1),reverse=True),
        (value%4 for value in gen_indexes_for_orders_2D(index,order))))






def get_states_at_orders_2D(index,order,fold_indexes_at_orders=None):
    starting_state_index = 0
    permutations = PERMUTATION_KEYS_2D[starting_state_index]
    if fold_indexes_at_orders is None:
        fold_indexes_at_orders = get_fold_indexes_at_orders_2D(index,order)

    next_state = starting_state_index
    states_at_orders = []
#    states_at_orders = [starting_state_index]
    for i in range(0,order+1):
        states_at_orders.append(next_state)
        fold_index = fold_indexes_at_orders[i]
        next_state = permutations[fold_index]
        permutations = PERMUTATION_KEYS_2D[next_state]

#    import gaepdb;gaepdb.set_trace()
        #experiment in list comprehension
#    states_at_orders = [next_state for next_state in [permutations[fold_index]] for fold_index in [fold_indexes_at_orders[i]] for i in range(0,order)]



    return states_at_orders

def get_dimensional_positions_at_orders_2D(index,order):
    index = int(index)
    fold_indexes_at_orders = get_fold_indexes_at_orders_2D(index,order)
    states_at_orders = get_states_at_orders_2D(index,order,fold_indexes_at_orders)
    positions_at_orders = {}
    for key,value in enumerate(states_at_orders):
        try:
            positions_at_orders[key]=CURVESTATES2D[value][fold_indexes_at_orders[key]]
        except TypeError, e:
            print "ERROR"
            print "value: ", value
            print "key: ",key
            raise e

    return positions_at_orders

def compare_range(a,b):
    if b[1] < a[1]:
        return 1
    elif b[1] == a[1]:
        return 0
    else:
        return -1

def compare_index(x, y):
    if y.index < x.index:
        return 1
    elif y.index == x.index:
        return 0
    return -1

def compare_range_length_desc(a,b):
    length_a = a[1].index-a[0].index
    length_b = b[1].index-b[0].index
    if length_b > length_a:
        return 1
    elif length_b == length_a:
        return 0
    else:
        return -1

def compress_ranges_2D(ranges):
    ranges.sort(compare_range)
    compact_ranges = []
    length = len(ranges)
    index = 0
    while index < length:
        subindex = index+1
        skip = 0
        while subindex<length:
            if ranges[index][1]+1 == ranges[subindex][0]:
                ranges[index]=(ranges[index][0],ranges[subindex][1])
                skip +=1
            else:
                compact_ranges.append(ranges[index])
                break
            if subindex==length-1:
                compact_ranges.append(ranges[index])
                break
            subindex+=1
        index += skip
        index += 1

    return compact_ranges

def curve_point_range_at_order_2D(curve_point,order):
    if curve_point.order >= order:
        return curve_point.index,curve_point.index
    else:
        total_count_at_order = 1<<(2*(order + 1))
        point_count_at_curve_order = 1<<(2*(curve_point.order + 1))
        count_at_order = total_count_at_order/point_count_at_curve_order
        starting_index = count_at_order*curve_point.index
        ending_index = starting_index+ count_at_order - 1
        return starting_index,ending_index

def long_distance_matching_lat_distance(lat_rads,distance_rads,
    miles_dist, radius = 3963.0):
    lat1 = lat_rads
    lon1 = 0.0
    lat2 = lat_rads
    lon2 = distance_rads
    distance = radius * math.acos(math.sin(lat1) *  math.sin(lat2) + \
    math.cos(lat1) * math.cos(lat2) * math.cos(lon2 - lon1))
    scale = distance_rads/distance*miles_dist
    return scale

def circumference_geopts(geopt,radius_in_miles=1.5,level_of_detail=4):
    try:
        from google.appengine.api.datastore_types import GeoPt
    except ImportError:
        from geopt_alt import GeoPt
    theta_lat = math.radians(geopt.lat)
    theta_lon = math.radians(geopt.lon)
    twopi = 2*math.pi
    if level_of_detail <= 0:
        level_of_detail = 4
    increment = (1.0/level_of_detail)*math.pi
    angle = twopi
    lat_distance = math.radians((1.0/69.1)*radius_in_miles)
    lon_distance = long_distance_matching_lat_distance(theta_lat,
        lat_distance,radius_in_miles)

    geopts = []
    while angle >= 0:
        delta1rads = math.cos(angle)*lat_distance+theta_lat
        delta1 = math.degrees(delta1rads)
        delta2 = math.degrees(math.sin(angle)*lon_distance+theta_lon)
        geopts.append(GeoPt(delta1,delta2))
        angle -= increment

    return geopts

CENTERPOINT_2D = np.array([0.5,0.5])
def first_order_curve_points():
    top_curve_points = []
    top_points = get_dimensional_position_of_children_2D(CENTERPOINT_2D,0,0)
    top_state_keys = PERMUTATION_KEYS_2D[0]
    for top_fold_index,state_key in enumerate(top_state_keys):
        statesKey = str(state_key)
        curve_key = str(top_fold_index)
        top_curve_points.append( HilbertSpaceCurvePoint(top_fold_index,
            top_points[top_fold_index],order=0,fold_index=top_fold_index,
            state_key=state_key,states_key=statesKey,curve_key=curve_key))

    return top_curve_points

FIRST_ORDER_CURVE_POINTS = first_order_curve_points()

def index_ranges_in_bounding_box_2D(bounding_box,to_order=None,
    added_precision=0,compress_ranges=True):
    """
"""

    checks = 0
    # find optimal order based on the size of box
    diagonal_size_of_box = bounding_box.diagonal_size()
    shiftstart = int(1.0/diagonal_size_of_box)
    starttime1 = time()
    if to_order is None:
        to_order = 0
        shift = shiftstart
        while shift > 0:
            shift = shiftstart>>to_order
            to_order += 1
    endtime1= time() - starttime1

    starttime2 = time()
    #logSomething = math.log(shiftstart,2) +2
    endtime2= time() - starttime2
    logging.debug("to_order time comp:%0.32f %0.32f"%(endtime1,endtime2,))
    #logging.debug("\nshiftstart:%s %s %s"%(shiftstart,to_order,int(logSomething)) )


    #adding or removing precision
    to_order += added_precision
    if to_order<0:
        to_order = 0

    open_list = copy.deepcopy(FIRST_ORDER_CURVE_POINTS)
    ranges = []
    logging.debug("to_order: %s"%to_order)
    startbboxtime = time()
    for order in range(0,to_order+1):
        if not len(open_list):
            break
        increment = distance_increment_at_order_2D(order+1)#INCREMENT_DISTANCES[order+1]
        next_order_open_list = []

        expanded_box = BoundingBox_2D(*bounding_box.expand_corners_by(increment))

        contracted_box = BoundingBox_2D(*bounding_box.contract_corners_by(increment))

        if contracted_box.is_inverted():
            for curve_point in open_list:
                checks+=1
                if expanded_box.point_in(curve_point.point):
                    if order == to_order:
                        ranges.append((curve_point.index,curve_point.index))
                    else:
                        next_order_open_list.extend(
                            get_children_curves_2D(curve_point))
        else:
            for curve_point in open_list:
                checks+=1
                if contracted_box.point_in(curve_point.point):
                    ranges.append(curve_point_range_at_order_2D(curve_point,to_order))
                elif expanded_box.point_in(curve_point.point):
                    checks+=1
                    if order == to_order:
                        ranges.append((curve_point.index,curve_point.index))
                    else:
                        next_order_open_list.extend(get_children_curves_2D(curve_point))

        open_list = next_order_open_list

    endbboxtime= time() - startbboxtime

    logging.debug("point-in-bbox checks %s, took %0.32f at %0.32f each"%(checks,
        endbboxtime,endbboxtime/checks))

    #compress ranges
    ranges.sort(compare_range)

    starttime3 = time()
    if compress_ranges:
        ranges = compress_ranges_2D(ranges)

    endtime3= time() - starttime3
    logging.debug("compress ranges: %0.32f"%(endtime3,))
    curve_ranges = []

    lenRanges = len(ranges)
    logging.debug("for %s ranges:"%(lenRanges,))
    starttime4 = time()
    get_point_time = 0.0
    get_point_count = 0

    index_point_time = 0.0
    index_point_count = 0

    cfTimeT = 0.0
    cfCountT = 0.0
    cfTime = 0.0

    for start,end in ranges:
        startpointtime = time()
        start_point = get_point_2D_at_index(start,to_order)
        get_point_time += time() - startpointtime
        get_point_count +=1

        startindextime = time()
        start_curve, timings = index_point_2D(start_point,to_order)
        cfTime,cfCount = timings
        cfTimeT += cfTime
        cfCountT += cfCount
        index_point_time += time() - startindextime
        index_point_count +=1

        if end == start:
            end_curve = copy.deepcopy(start_curve)
        else:
            startpointtime = time()
            end_point = get_point_2D_at_index(end,to_order)
            get_point_time += time() - startpointtime
            get_point_count +=1

            startindextime = time()
            end_curve, timings = index_point_2D(end_point,to_order)
            cfTime,cfCount = timings
            cfTimeT += cfTime
            cfCountT += cfCount
            index_point_time += time() - startindextime
            index_point_count +=1


        curve_ranges.append((start_curve,end_curve))

    endtime4 = time() - starttime4
    logging.debug("index ranges:%0.32f, that's %0.32f/range"%(endtime4,endtime4/lenRanges))
    logging.debug("get point time:%0.32f, that's %0.32f/curve of %s"%(get_point_time,
        get_point_time/get_point_count, get_point_count))
    logging.debug("index point time:%0.32f, that's %0.32f/point of %s"%(index_point_time,
        index_point_time/index_point_count,index_point_count))

    logging.debug("get_next_order_fold_index_and_point_2D time:%0.32f, that's %0.32f/point of %s\n"%(cfTimeT,
        cfTime/cfCountT,cfCountT))


    curve_ranges.sort(compare_range_length_desc)


    return curve_ranges

##

def get_common_curve_point(start,end):
    start_key = start.get_curve_key()
    end_key = end.get_curve_key()
    level = 0
    start_len = len(start_key)
    if start_len == len(end_key):
        for index in range(0,start_len):
            logging.info("? %s == %s"%(start_key[index],end_key[index]))
            if not start_key[index] == end_key[index]:
                logging.info("matched to level %s"%level)
                break
            level = index


    common_point = copy.deepcopy(start)

    common_point.order = order
    common_point.fold_index =common_point.curve_key[-1]
    common_point.state_key = common_point.states_key[-1]
    common_point.index = common_point.index_from_curve_key(4)
    common_key = common_point.get_curve_key()
    logging.info("common level of %s and %s is %s with key: %s"% \
        (start_key,end_key,common_point.level,common_key))
    #raise ValueError("blah")
    return common_point


def caps_from_curve_point_range_sets(curve_point_range_sets):
    """
    geocap is a curvepoint range reduced to the common curvepoint

    10230240321001,1023024031231 => 10230240

    curve_point_range_sets=[{'face':face, \
            'ranges':spec.index_ranges_in_bounding_box(query_bbox, \
            added_precision=added_precision)}]
    """
    geocaps = []

    for range_set in curve_point_range_sets:
        for curve_point_range in range_set['ranges']:
            geocaps.append(get_common_curve_point(curve_point_range[0],
                curve_point_range[1]))

    return geocaps

def endcaps_from_curve_point_range_sets(curve_point_range_sets):
    """
    geocap is a curvepoint range reduced to the common curvepoint

    10230240321001,1023024031231 => 10230240

    curve_point_range_sets=[{'face':face, \
            'ranges':spec.index_ranges_in_bounding_box(query_bbox, \
            added_precision=added_precision)}]
    """
    geocaps = []

    for range_set in curve_point_range_sets:
        for curve_point_range in range_set['ranges']:
            geocaps.append(curve_point_range[0])
            geocaps.append(curve_point_range[1])

    return geocaps

def capkey_boundary_multipolygon(geopt,order,scale_factor=0.99999999999,res=1):

    cfc = sphere_to_cube_face_2D(geopt)
    order -= 1
    cfcs = cfc.boundary_face_positions_at_order(order,
            scale_factor=scale_factor,res=res)
    #exclude center point
    points=[cfc.to_sphere_geopt() for cfc in cfcs]
    return points

def center_geopt_at_order(geopt,order):
    cfc = sphere_to_cube_face_2D(geopt)
    return center_at_order_to_sphere_geopt(order)

def geoPtToGeojsonable(geopt):
    #{ "type": "Point", "coordinates": [100.0, 0.0] }
    if geopt and hasattr(geopt,'lat') and hasattr(geopt,'lon'):
        return { "type": "Point", "coordinates": [geopt.lon, geopt.lat] }
    return {}

def geoPtFromGeoJsonalbe(geojsonalbe):
    try:
        from google.appengine.api.datastore_types import GeoPt
    except ImportError:
        from geopt_alt import GeoPt
    #{ "type": "Point", "coordinates": [ -93.491307000098629, 36.089674000689349, 0.0 ] }
    if geojsonalbe and geojsonalbe.get('type') == "Point":
        coords = geojsonalbe.get('coordinates')
        if coords:
            if len(coords) > 1:
                return GeoPt(coords[1],coords[0])
    else:
        logging.debug("couldn't convert %s"%geojsonalbe)

def capkey_boundary_from_encoded_capkey(encoded_capkey):
    cfcF = CubeFaceAndCoords.from_capkey(encoded_capkey)
    order = len(encoded_capkey)-2
    cfcs = cfcF.boundary_face_positions_at_order(int(order))
    #exclude center point
    points=[cfc.to_sphere_geopt() for cfc in cfcs]
    multipolygon = [[[ [point.lon, point.lat] for point in points]]]
    geojsonable ={ "type": "MultiPolygon", "coordinates": multipolygon}
    return geojsonable


def main():
    try:
        from google.appengine.api.datastore_types import GeoPt
    except ImportError:
        from geopt_alt import GeoPt
    state0m = STATE0M
    R90CW,FH,FV,R90CW_FH,R90CW_FV,R180CW = get_standard_matrixes()

    print "\nstate0"


    state0 = np.round_(state0m.getT()[:2:].getT().astype(int)).astype(int)
    print state0


    #print ""
    #rot90 = R90CW.dot(state0m.getT()).getT()
    #print rot90
    #print ""
    #flipH = FH.dot(rot90.getT()).getT()
    #print flipH

    #print ""

    print "\nstate0 rotated clockwise"
    print np.round_(R90CW.dot(state0m.getT())[:2:].getT()).astype(int)


    state1 = np.round_(R90CW_FH.dot(state0m.getT())[:2:].getT()).astype(int)
    print "\nstate1"
    print state1



    #print np.allclose(flipH,state1)



    state2 = np.round_(R90CW_FV.dot(state0m.getT())[:2:].getT()).astype(int)

    print "\nstate2"
    print state2

    print "\nstate3"


    state3 = np.round_(R180CW.dot(state0m.getT())[:2:].getT()).astype(int)

    print state3

    #pad test
    print "\npad test"
    print pad_state2d(state3)

    CURVESTATES = np.array([state0,state1, \
    state2, state3])
    #print "CURVESTATES = "
    #print CURVESTATES


    #base4 for incrementing a capkey
    print np.base_repr(23,4,10)
    #get order of capkey,capkey to base10, increment,
    #result to base4 padded to the order

    #import pdb;pdb.set_trace()

    #make capkey from point and order
    #   face_point = sphere to cube face 2d
    #   face point encoded

    print "CURVESTATES2D = create_curve_states_2D()"
    print CURVESTATES2D.tostring()


    print "\nstate 2 permutated"
    state2_permed =  permutate_curve_state_matrix(CURVESTATES2D[2])
    print state2_permed

    print "\nindex of states in state2_permed"
    for curve_state in state2_permed:
        print index_of_state(curve_state)

    print "\nindex of states in CURVESTATES2D"
    for curve_state in CURVESTATES2D:
        print index_of_state(curve_state)

    print ""

    print "PERMUTATION_KEYS_2D = ",
    print PERMUTATION_KEYS_2D



    ageopoint = GeoPt(36.33459, -94.09610)
    print ageopoint
    capkeyA = get_capkey(ageopoint)
    print capkeyA
    cfcF = CubeFaceAndCoords.from_capkey(capkeyA)
    print "from capkey to position ",
    print cfcF.to_position3()
    print "SPHERE GEOPOINT ",cfcF.to_sphere_geopt()
    cfc = sphere_to_cube_face_2D(ageopoint)
    capkeyB = index_point_at_orders_2D(cfc.position,MAX_ORDER)
    print capkeyB
    print "".join((str(ck) for ck in capkeyB)) == capkeyA[1:]
    #import pdb;pdb.set_trace()
    print cfc.to_position3()
    print "SPHERE GEOPOINT ",cfc.to_sphere_geopt()


    mpoly =  capkey_boundary_multipolygon(ageopoint,order=10)

    print mpoly

    print get_surrounding_capkeys(ageopoint,order=10)


    points = capkey_boundary_multipolygon(ageopoint,10)
    multipolygon = [[[ [point.lon, point.lat] for point in points]]]
    geojsonable ={ "type": "MultiPolygon", "coordinates": multipolygon}
    bbox = BoundingBox_2D.from_points(
        [np.array([point.lat,point.lon]) for point in points])

    center_point = bbox.get_center()
    center_geopt = GeoPt(center_point[0],center_point[1])
    print geojsonable
    print geoPtToGeojsonable(center_geopt)




if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()

