""" IO classes for the DTU WAsP file types

Copyright (C) 2013 DTU Wind Energy

Authors: Pierre-Elouan Rethore
Email: pire@dtu.dk
Last revision: 31/10/2013

License: Apache v2.0, http://www.apache.org/licenses/LICENSE-2.0
"""
from numpy.ma import argmin, array

from we_file_io import WEFileIO, TestWEFileIO
from xml.dom.minidom import parseString
import xml.etree.cElementTree as ET
from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement, dump
import unittest

from matplotlib import pylab as plt

import numpy as np


#### Semi private classes --------------------------------------------------------------------------


def print_elt(e_):
    print e_.tag
    print e_.text
    print [i.tag for i in e_]
    print e_.attrib
    print e_.tail
    
children = lambda e_: [i.tag for i in e_]
find_attr = lambda e_, attrib_, text_: filter(lambda e__: e__.get(attrib_) == text_, e_.getiterator())
find_tag = lambda e_, text_: filter(lambda e__: e__.tag == text_, e_.getiterator())
findall_attr = lambda e_, attrib_: map(lambda e__: e__.get(attrib_),
                               filter(lambda e__: e__.get(attrib_) != None, e_.getiterator()))

sectordata = lambda g: np.array(list(map(lambda d: [float(d.get(i)) for i in ('CentreAngleDegrees', 'SectorFrequency', 'WeibullA', 'WeibullK')], g)))
datapoint = lambda g: np.array(list(map(lambda d: [float(d.get(i)) for i in ('WindSpeed', 'PowerOutput', 'ThrustCoEfficient')], g)))
location = lambda g: np.array([float(g.get(i)) for i in ('x-Location', 'y-Location')])

def get_wt(wtg):
    a = datapoint( find_tag(wtg, 'DataPoint'))
    ## Sorting with respect of wind speed
    turbine = {}
    turbine['name'] = wtg.get('Description')
    wtg2 = wtg.find('MemberData/WindTurbineGenerator')
    turbine['data'] = a[np.argsort(a[:,0]),:]
    turbine['density'] = float(wtg2.find('PerformanceTable').get('AirDensity'))
    turbine['rotor_diameter'] = float(wtg2.get('RotorDiameter'))
    turbine['hub_height'] = wtg2.findall('SuggestedHeights/Height')[0].text
    turbine['manufacturer'] = wtg2.get('ManufacturerName')
    return turbine


def generate_vl(arr_):
    """Iterator for populating a list of VectLines"""
    n1 = 0
    while n1 < len(arr_):
        vl = VectLine(arr_[n1:])
        n1 += vl.n_end
        yield vl

class VectLine(object):
    """Contains a list of points and their respective height"""
    def __init__(self, arr_):
        self.h = arr_[0]
        self.n = int(arr_[1])
        self.n_end = 2+2*self.n
        self.points = arr_[2:self.n_end].reshape([-1,2])
        
    def plot(self, scale=1000.0, colmap=plt.cm.jet, **kwargs):
        """Plot the vectorline"""
        plt.plot(self.points[:,0], self.points[:,1], color=colmap(self.h/scale), **kwargs)
        
    def add_to_wasp(self, wasp_core):
        """Add the vectorline to the wasp core"""
        wasp_core.addorographicline(self.h, self.n, self.points.T)

    def write(self, fid):
        fid.write('%f  %d\n'%(self.h, self.n))
        fid.write(' '.join([str(i) for i in self.points.flatten()]) + '\n')
            

### File I/O Classes ------------------------------------------------------------------------------


class MAP(WEFileIO):
    """WAsP MAP File."""

    ### Private methods to be implemented in the subclasses --------------------
    def _read(self):
        """ Read the map file. Place a list of vector lines in self.data"""
        with open(self.filename, 'r') as f:
            map_str=f.readlines()
        self._finish_reading(map_str)

    def _finish_reading(self, map_str):    
        self.header = map_str[:4]
        arr1 = np.array(''.join(map_str[4:]).split(), dtype='float')
        ### Data contains a list of vector lines
        self.data = list(generate_vl(arr1))
        ### Maximum height of all the vector lines
        self.max_height = np.array([v.h for v in self.data]).max()        

    def _write(self):
        """ Write a file, with the same header"""
        with open(self.filename, 'w') as f:
            f.write(''.join(self.header))
            for v in self.data:
                v.write(f)

    def plot(self, **kwargs):
        """Plot all the vector lines. Scale their color with the height. 
        Returns a list of all the plot handles.
        """
        return [vl.plot(scale=self.max_height, **kwargs) for vl in self.data]
            
    def add_to_wasp(self, wasp_core):
        """Add all the vector lines to the wasp core"""
        for v in self.data:
            v.add_to_wasp(wasp_core)

class ZipMAP(MAP):
    def __init__(self, filename=None, sub_filename=None):
        """ Initialized the classe using the filename

        Parameters:
        ----------
        filename : string (optional)
                   The file name to read and write
        """
        if sub_filename:
            self.sub_filename = sub_filename
            
        super(ZipMAP, self).__init__(filename)

    def _read(self):
        with ZipFile(self.filename, 'r') as f:
            map_str=f.read(self.sub_filename).splitlines()
        self._finish_reading(map_str)


class WTG(WEFileIO):
    """WAsP Turbine File."""

    
    ### Private methods to be implemented in the subclasses --------------------
    def _read(self):
        """ Read the file."""
        xml = ET.parse(self.filename).getroot()
        self.xml = xml

        # There might be several power curves for different air density
        data_sets = []
        densities = []
        for pt in xml.findall('PerformanceTable'):
            densities.append(float(pt.get('AirDensity')))
            data_sets.append(datapoint(pt.findall('DataTable/DataPoint')))

        # Looking for the density closest to normal atmospheric conditions: 1.225
        id = argmin(abs(array(densities) - 1.225))
        self.density = densities[id]
        a = data_sets[id]

        ## Sorting with respect of wind speed
        self.data = a[np.argsort(a[:,0]),:]
        self.rotor_diameter = float(xml.get('RotorDiameter'))
        self.hub_height = float(xml.findall('SuggestedHeights/Height')[0].text)
        self.manufacturer = xml.get('ManufacturerName')

    def _write(self):
        """ Write a file"""
        ### You are going to replace this code when you inherit from this class
        dtable = self.xml.find('PerformanceTable/DataTable')
        #remove the existing datapoints
        dtable.clear()
        for ws, p, ct in self.data:
            SubElement(dtable, 'DataPoint', WindSpeed=str(ws), PowerOutput=str(p), ThrustCoEfficient=str(ct))

        with open(self.filename, 'w' ) as f:
            f.write( '<?xml version="1.0"?>' )
            f.write( ElementTree.tostring( self.xml ) )

class POW(WEFileIO):
    """WAsP POW Turbine File."""

    
    ### Private methods to be implemented in the subclasses --------------------
    def _read(self):
        """ Read the file."""
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        self.name = lines[0].split('\n')[0]
        self.hub_height = float(lines[1].split()[0])
        self.rotor_diameter = float(lines[1].split()[1])
        self.data = np.array([l.split()[:3] for l in lines[3:]], dtype='float32')
        factors = np.array(lines[2].split()[:2], dtype='float32')
        ### Multiplying by the factors to get the right units
        self.data[:,0:2] *= factors
        self.factors = factors

    def _write(self):
        # """ Write a file"""
        # ### You are going to replace this code when you inherit from this class
        with open(self.filename, 'w') as f:
            f.write(self.name + '\n')
            f.write('%f %f\n'%(self.hub_height, self.rotor_diameter))
            f.write('%f %f\n'%(self.factors[0], self.factors[1]))
            for i in range(self.data.shape[0]):
                f.write('%f %f %f\n'%(self.data[i,0] / self.factors[0], 
                                      self.data[i,1] / self.factors[1], 
                                      self.data[i,2]))

class WWF(WEFileIO):
    """WAsP Wind Farm Site."""

    ### Private methods to be implemented in the subclasses --------------------
    def _read(self):
        """ Read the file."""
        xml = ET.parse(self.filename).getroot()
        self.xml = xml
        SiteSummary = lambda d: (d.get('Label'),[float(d.get(i)) for i in ('XLocation', 'YLocation', 'HeightAGL', 'SiteElevation')])
        SectorData = lambda d: [float(d.get(i)) for i in ('CentreAngle', 'Frequency', 'WeibullA', 'Weibullk')]
        SectorWiseData = lambda d: (d.find('SiteSummary').get('Label'), 
                                    np.array(list(map(SectorData, d.findall('PredictedWindClimate/SectorWiseData/SectorData')))))

        self.data = dict(list(map(SiteSummary, xml.findall('TurbineSite/SiteSummary'))))
        self.windroses = dict(list(map(SectorWiseData, xml.findall('TurbineSite'))))
        self.pos = np.array([(d[0], d[1]) for d in self.data.values()])

    def _write(self):
        """ Write a file"""
        ### You are going to replace this code when you inherit from this class
        raise NotImplementedError("This method hasn't been implemented yet")

    def plot(self):
        """ Plot the position of the wind turbines """
        plt.plot(self.pos[:,0], self.pos[:,1] , '.')


from zipfile import ZipFile
class WWH(WEFileIO):
    """WAsP Workspace file .wwh"""

    def _read(self):
        """Unzip and read the file"""
        with ZipFile(self.filename, 'r') as f:
            e = ET.fromstring(f.read('Inventory.xml'))
        self.e_ = e


        wind_turbines = {}
        self.e_positions = {}
        for wt in find_attr(e, 'ClassDescription', 'Turbine site'):
            wt_name = wt.get('Description')
            height = find_tag(wt, 'SiteInformation')[0].get('WorkingHeightAgl')
            position = location(find_tag(wt,'Location')[0]).tolist()
            wind_turbines[wt_name] = {}
            ### Combine the position and the height into one vector
            wind_turbines[wt_name]['position'] = np.array(position + [float(height)])
            wind_climates = find_tag(wt, 'WeibullWind')
            if len(wind_climates) == 0:
                print 'no predicted wind climates'
            else:
                a = sectordata(wind_climates)
                wind_turbines[wt_name]['wind_rose'] = a[np.argsort(a[:,0]),:]
            self.e_positions[wt_name] = find_tag(wt,'Location')


            #wind_turbines[wt_name] = 
            
        ### Look for how many wind turbine generators there is:
        wtgs = find_attr(e, 'ClassDescription', 'Wind turbine generator')
        turbine_descriptions = {}
        if len(wtgs) == 0:
            print 'No WTG present in this file'
        else:
            for wtg in wtgs:
                turbine = get_wt(wtg)
                turbine_descriptions[turbine['name']] = turbine
                print turbine['name']


        for wt_name in wind_turbines.keys():
            wind_turbines[wt_name]['type'] = turbine['name']


        site_groups = find_attr(e, 'ClassDescription', 'Turbine site group')
        wt_groups = {}

        if len(site_groups)>1:
            # Find the main turbine site group
            for ts in site_groups:
                sub_site_groups =  filter(lambda i: i.get('Description') != ts.get('Description'), find_attr(ts, 'ClassDescription', 'Turbine site group'))
                if len(sub_site_groups) > 0:
                    # It's the main site
                    for sts in sub_site_groups:
                        group_name = sts.get('Description')
                        wt_groups[group_name] = {}
                        wt_groups[group_name]['wind_turbines'] = []
                        for wt in find_attr(sts, 'ClassDescription', 'Turbine site'):
                            wt_name = wt.get('Description')
                            wt_groups[group_name]['wind_turbines'].append(wt_name)
                            # Look for existing turbines in the sub_groups
                            wtgs = find_attr(e, 'ClassDescription', 'Wind turbine generator')
                            if len(wtgs) == 0:
                                wind_turbines[wt_name]['type'] = turbine['name']
                            else:
                                wind_turbines[wt_name]['type'] = get_wt(wtgs[0])['name']
                        # Look for Generalised wind climate
                        gwcs = find_attr(sts, 'ClassDescription', 'Generalised wind climate')
                        general_windrose = {}
                        for gwc in gwcs:
                            for h in find_tag(gwc,'WindAtlasWeibullWindRose'):
                                height = h.get('ReferenceHeight')
                                wind_climates = h.findall('WeibullWind')
                                a = sectordata(wind_climates)
                                general_windrose[height] = a[np.argsort(a[:,0]),:]
                        wt_groups[group_name]['general_wind_roses'] = general_windrose
                        
        #map_name = find_tag(e, 'ExternalArchive')[0].get('ArchiveTagID')
        #map_file = ZipMAP(wwh_file, map_name)

        self.wt_groups = wt_groups
        self.wind_turbines = wind_turbines
        self.data = wind_turbines
        self.turbine_descriptions = turbine_descriptions
        self.pos = np.array([d['position'][0:2] for d in wind_turbines.values()])



class old_WWH(WEFileIO):
    """WAsP Workspace file .wwh"""

    def _read(self):
        """Unzip and read the file"""
        with ZipFile(self.filename, 'r') as f:
            e = ET.fromstring(f.read('Inventory.xml'))

        ### Simple function to look for a specific Class Descriptor. Return a list.
        xmlf = lambda exml, keyword, address: filter(lambda x: x.get('ClassDescription') == keyword, exml.findall(address))

        ### Find the turbine sites
        turbine_sites = xmlf(e, 'Turbine site group', 
                             'WaspHierarchyMember/ChildMembers/WaspHierarchyMember/ChildMembers/WaspHierarchyMember')

        self.general_windroses = {}
        self.data = {}
        self.windroses = {}
        self.sites = {}

        for ts in turbine_sites:
            site_name = i.get('Description')
            self.sites[site_name] = {}
            try:
                ### Get the wind turbine description
                turbine_xml = xmlf(ts, 'Wind turbine generator', 'ChildMembers/WaspHierarchyMember')[0].find('MemberData/WindTurbineGenerator')
                a = datapoint(turbine_xml.findall('PerformanceTable/DataTable/DataPoint'))
                ## Sorting with respect of wind speed
                self.turbine = {}
                self.turbine['data'] = a[np.argsort(a[:,0]),:]
                self.turbine['density'] = float(turbine_xml.find('PerformanceTable').get('AirDensity'))
                self.turbine['rotor_diameter'] = float(turbine_xml.get('RotorDiameter'))
                self.turbine['hub_height'] = turbine_xml.findall('SuggestedHeights/Height')[0].text
                self.turbine['manufacturer'] = turbine_xml.get('ManufacturerName')
                self.sites[site_name]['turbine'] = self.turbine
            except:
                print 'No wind turbine generator found'



            ### Data contains the label name of the turbine as keys, and for each one a list of x,y,h
            ts_xml = xmlf(ts, 'Turbine site', 'ChildMembers/WaspHierarchyMember')

            
            
            SectorData = lambda d: [float(d.get(i)) for i in ('CentreAngleDegrees', 
                                                              'SectorFrequency', 
                                                              'WeibullA', 
                                                              'WeibullK')]

            for i in list(ts_xml):
                label = i.get('Description')
                site_info = i.find('MemberData/SiteInformation')
                location = i.find('MemberData/SiteInformation/Location')
                self.data[label] = [float(location.get('x-Location')),
                                    float(location.get('y-Location')),
                                    float(site_info.get('WorkingHeightAgl'))]
                wind_climates = i.findall('MemberData/CalculationResults/PredictedWindClimate/RveaWeibullWindRose/WeibullWind')
                if len(wind_climates) == 0:
                    print 'no predicted wind climates'
                else:
                    a = np.array(list(map(SectorData, wind_climates)))
                    self.windroses[label] = a[np.argsort(a[:,0]),:]


            ts_xml = xmlf(ts, 'Generalised wind climate', 'ChildMembers/WaspHierarchyMember')
            for i in list(ts_xml):
                label = i.get('Description')
                self.general_windroses[label] = {}
                for h in i.findall('MemberData/RveaGeneralisedMeanWindClimate/WindAtlasWeibullWindRose'):
                    height = h.get('ReferenceHeight')
                    wind_climates = h.findall('WeibullWind')
                    a = np.array(list(map(SectorData, wind_climates)))
                    self.general_windroses[label][height] = a[np.argsort(a[:,0]),:]

            ### 2D array of x,y for each turbine
            self.pos = np.array([(d[0], d[1]) for d in self.data.values()])

        
## Do Some testing -------------------------------------------------------
class TestWAsP(TestWEFileIO):
    """ Test class for MyFileType class """

    test_wtg = 'test/wasp/V80-2MW-offshore.wtg'
    test_pow = 'test/wasp/bonus450.pow'
    test_wwf = 'test/wasp/hornsrev1.wwf'
    test_map = 'test/wasp/WaspMap.map'

    def test_WTG_duplication(self):
        self._test_duplication_array(WTG, self.test_wtg)

    def test_WTG_several_power_curves(self):
        wtg = WTG(self.test_wtg)
        # Make sure that we don't get several power curves bundled together (when there are several power curves defined
        # for different air densities)
        self.assertNotAlmostEqual(wtg.data[0,0], wtg.data[1,0])

    def test_POW_duplication(self):
        self._test_duplication_array(POW, self.test_pow)

    def test_MAP_dupplication(self):
        original_file, new_file = self._duplicate(MAP, self.test_map)
        for of, nf in zip(original_file.data, new_file.data):
            self.assertTrue(np.linalg.norm(of.points-nf.points)<1.0E-8)

## Main function ---------------------------------------------------------
if __name__ == '__main__':
    """ This is the main fuction that will run the tests automatically

    $> python my_file_type.py

    ----- SOME_PRINT_STATEMENTS -----
    .
    ----------------------------------------------------------------------
    Ran X test in XXXs

    OK
    """
    unittest.main()