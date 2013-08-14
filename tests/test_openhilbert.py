import unittest
import sys,os

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))
sys.path.append(ROOT_PATH)

__author__ = 'john'

from openhilbert import *

try:
    from google.appengine.api.datastore_types import GeoPt
except ImportError:
    from geopt_alt import GeoPt

state0m = STATE0M
R90CW,FH,FV,R90CW_FH,R90CW_FV,R180CW = get_standard_matrixes()





class TestHilbertTables(unittest.TestCase):
    def setUp(self):
        self.STATE0M = np.matrix([[0, 0, 0, 1],
            [0, 1, 0, 1],
            [1, 1, 0, 1],
            [1, 0, 0, 1]])



    def test_permutations(self):
        self.assertEqual(PERMUTATION_KEYS_2D,create_permutation_keys())


    def test_state_zero(self):
        good_state_zero = np.array([[0,0],
                                    [0,1],
                                    [1,1],
                                    [1,0]])



        state0 = np.round_(self.STATE0M.getT()[:2:].getT().astype(int)).astype(int)
        result = state0 == good_state_zero
        for val in (row for row in result.tolist()):
            self.assertTrue(val)


def clean_bounds(boundsinput):

    bounds = boundsinput
    bounds = bounds.split(',')
    bbox = None
    if len(bounds) != 4:
        raise ValueError("bounds must be 4 values")
    try:

        bbox = BoundingBox_2D(np.array([float(bounds[1]),float(bounds[0])]),\
            np.array([float(bounds[3]),float(bounds[2])]))
    except Exception, e:
        logging.exception("invalid bounds")
        raise
    return bbox

def bounding_box_curvepoint_range_sets(bbox,to_order=None,added_precision=0):
    curve_point_ranges_sets = []
    corner_geopts = bbox.corner_geopts()
    fpSW = sphere_to_cube_face_2D(corner_geopts[0])
    fpNW = sphere_to_cube_face_2D(corner_geopts[1])
    fpNE = sphere_to_cube_face_2D(corner_geopts[2])
    fpSE = sphere_to_cube_face_2D(corner_geopts[3])
    if fpSW.face == fpNW.face and fpNW.face == fpNE.face\
    and fpNE.face == fpSE.face:

        query_bbox = boundingBox_from_points_2D([ fpSW.get_point(),
                                                  fpNW.get_point(),
                                                  fpNE.get_point(),
                                                  fpSE.get_point()])
        curve_point_ranges_sets.append({'face':fpSW.face,\
                                        'ranges':index_ranges_in_bounding_box_2D(query_bbox,\
                                            to_order=to_order,added_precision=added_precision)})

    else:
        pass

    return curve_point_ranges_sets


def capkey_ranges_from_curve_point_range_sets(curve_point_ranges_sets):
    capkey_ranges = []
    for range_set in curve_point_ranges_sets:
        capkey_ranges.extend(capkey_ranges_from_curve_point_ranges(\
            range_set['face'],\
            range_set['ranges']))
    return capkey_ranges


class TestHilbertQuery(unittest.TestCase):
    def setUp(self):
        #http://127.0.0.1:8080/zip_roads_info/bounds?bounds=36.44251723118404,-94.12539839744568,36.44596949306726,-94.10810351371765&z=17&vc=5&lp=
        bbox = clean_bounds("36.44251723118404,-94.12539839744568,36.44596949306726,-94.10810351371765")
        zoom = int("17")
        to_order = None
        if zoom > 3:
            to_order = zoom - 3
        self.bbox = bbox
        self.to_order = to_order
        self.accurate_ranges = [('3022202021000320', '3022202021001000'),
            ('3022202021000002', '3022202021000020'),
            ('3022202012333003', '3022202012333012'),
            ('3022202012333322', '3022202012333331'),
            ('3022202012333000', '3022202012333001')]


    def test_ranges(self):
        logging.debug("start finding curve_point_ranges")
        curve_point_ranges_sets = bounding_box_curvepoint_range_sets(self.bbox,
            to_order=self.to_order)

        logging.debug("start capkey_ranges_from_curve_point_range_sets")
        capkey_ranges = capkey_ranges_from_curve_point_range_sets(\
            curve_point_ranges_sets)


        #import pdb;pdb.set_trace()
        for index,val in enumerate(capkey_ranges):
            self.assertEqual(val,self.accurate_ranges[index])

    def test_boundary_to_geojson(self):
        mpgeojsons = [capkey_boundary_from_encoded_capkey('3022202012333012') \
        for index in range(0,200)]
        self.assertEqual(len(mpgeojsons),200)




if __name__ == '__main__':
    unittest.main()