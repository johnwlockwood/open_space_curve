import unittest

__author__ = 'john'

import sys,os
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))
sys.path.append(ROOT_PATH)


import openhilbert
try:
    from google.appengine.api.datastore_types import GeoPt
except ImportError:
    from geopt_alt import GeoPt



class TestGet_capkey(unittest.TestCase):
    def setUp(self):
        self.latlon_matrix = [
            [36.334590023456,-94.096012345678],
            [51.554569329689, -121.95183684596],
            [62.165573, -145.898569],
            [-35.433468, 149.930205],
            [-35.433724, 149.930272],
            [-35.433667, 149.930039],
            [-15.812874, 144.557934],
            [-8.795315, -39.334067],
            [-52.080891, -60.397449],
            [67.995442, 27.869721],
            [7.711542, -2.285379],
            [70.700128, 137.614534],
            [-77.235310, -34.452937],
            [69.209495, -157.192715],
            [69.2094948967,-157.192712799]

        ]

    def latlon_to_geopt_to_capkey_to_geopt(self,lat,lon,order = 29):
        startgeopt = GeoPt(lat, lon)
        capkeyA = openhilbert.get_capkey(startgeopt,order=order)
        #print startgeopt
        #print capkeyA
        cfcF = openhilbert.CubeFaceAndCoords.from_capkey(capkeyA)
        endgeopoint = cfcF.to_sphere_geopt()

        return startgeopt, endgeopoint

    def convert_geopt(self,lat,lon,order = 29,precision_delta = 0.0000001016):
        startgeopt, endgeopoint = self.latlon_to_geopt_to_capkey_to_geopt(lat,lon,order=order)

#        print endgeopoint.lat - startgeopt.lat
#        print endgeopoint.lon - startgeopt.lon
        #0.0000001 == 1.11 centimeters, times 2.1 = 0.917716535 inches
        #0.00001 == 1.11 meters, times 2.6 = 113.622047 inches

        self.assertAlmostEqual(startgeopt.lat,endgeopoint.lat,delta=precision_delta)
        self.assertAlmostEqual(startgeopt.lon,endgeopoint.lon,delta=precision_delta)



    def test_get_capkey(self):
        for lat,lon in self.latlon_matrix:
            self.convert_geopt(lat,lon,order=29)

    def test_capkey_reconvert_precision(self):
        order = 30
        self.assertEqual(order,30)
        for lat,lon in self.latlon_matrix:
            startgeopt, endgeopoint = self.latlon_to_geopt_to_capkey_to_geopt(lat,lon,
                order = order)
#            print "precision test"
            startgeopt2, endgeopoint2 = self.latlon_to_geopt_to_capkey_to_geopt(endgeopoint.lat,
                endgeopoint.lon,order = order)
#            print startgeopt
#            print endgeopoint
#            print endgeopoint2
#            print "lat diff: ",endgeopoint2.lat - startgeopt2.lat
#            print "lon diff: ",endgeopoint2.lon - startgeopt2.lon
            self.convert_geopt(endgeopoint2.lat,endgeopoint2.lon,order=order)


if __name__ == '__main__':
    unittest.main()