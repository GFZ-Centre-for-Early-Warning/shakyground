import warnings; warnings.filterwarnings("ignore")
import h5py
import numpy as np
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from openquake.hazardlib.geo import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.imt import PGA, PGV, SA
import shakemap_lite as sml
import quakeml
import shakeml
import pandas
import argparse

class rpe:
    '''
    returns epi-central distance based on simple gmpe
    '''
    def __init__(self,c0,c1,c2):
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2

    def get_repi(self,mw,pga):
        return np.exp((np.log(pga)-self.c0-self.c1*mw)/self.c2)

def get_distance(mw,pga):
    '''
    get epicentral distance using
    simple GMPE
    ln(pga) = c0 + c1*magnitude + c2*ln(repi)
    derived from NGA-West2 data
    '''
    c0 = -6.553383603995254
    c1 =  1.6606638130487288
    c2 = -1.7235855714775898
    instance = rpe(c0,c1,c2)
    return instance.get_repi(mw,pga)

def filter_sites(sites,roi):
    '''
    filters sites based on roi
    '''
    #use pandas dataframe
    sites = pandas.DataFrame(sites)
    #filter
    return sites[(sites.lon >= roi[0]) & (sites.lon <= roi[1]) & (sites.lat >= roi[2]) & (sites.lat <= roi[3])].to_dict('list')

def event_shakemap(quakemlfile,imt = "PGA",gmpe = "BindiEtAl2014Rjb", sites=None, roi=None, pgamin=None):
    '''
    Takes an event defined in quakeml: quakemlfile (can be file or string)
    and returns a shakemap calculated for
    an intensity type(available in OpenQuake): imt
    using a GMPE(available in OpenQuake): gmpe
    at sites that can be defined through:
     sites: a shakeml sites-file at least columns longitude, latitude, vs30, (optional: vs30measured, z1pt0, z2pt5, backarc)
     roi  : a region of interest [lonmin,lonmax,latmin,latmax]
            (if no additional sites are provided, will use global vs30 grid from USGS), can be used to filter provided sites
     pgamin: sites can be selected from usgs vs30 grid based on
             the epicentral distance, where pga typically drops below pgamin (using a simple global GMPE from the NGA-West2)
             Note: This is only a rough estimate and might be very different from the actual shakemap which may use a more sophisticated GMPE
     All three ways can be combined also and the intersection of all sets will be used
     if none of the above is provided the latter is used with a default pgamin of 0.1g
    '''

    #override the imt parameter - TEMPORARY MODIFICATION
    imt = ["PGA","SA(0.3)","SA(1.0)"]
    #####################

    #read event
    event = quakeml.quakeml2events(quakemlfile,provider='GFZ')
    #if more than one limit to first
    event = event.iloc[0]
    #create oq event
    eid = "EQ001"
    if any(np.isnan(list(event[['longitude','latitude','depth']]))):
        raise Exception('Need at lest hypocenter of event:\n{}'.format(event))
    if any(np.isnan(list(event[['strike','dip','rake']]))):
        oqevent = sml.Event(eid, event.longitude, event.latitude, hypo_depth=event.depth, mag=event.magnitude)
        rupture = False
    else:
        oqevent = sml.Event(eid, event.longitude, event.latitude, hypo_depth=event.depth, mag=event.magnitude, strike=event.strike, dip=event.dip, rake=event.rake)
        rupture = True
    #generate rupture
    rupt = oqevent.get_rupture()
    
    # Setup the shakemap tool with the desired GMPEs and - only needs to be done once.
    # Note that modification to consider different sets of GMPEs organised by tectonic region is simple
    shakemap = sml.ShakemapLite("dbfile.hdf5", [gmpe], imt)
    # Deletes any existing database of the same name
    shakemap.reset()


    if pgamin!=None:
        val=pgamin
    else:
        val=0.1
    #get some idea about distance to be considered for sites
    #d = 10**np.ceil(np.log10(get_distance(event.magnitude,pgamin)))
    d = get_distance(event.magnitude,val) #km
    dlat = d/111.2 #deg
    dlon = d/111.2*np.cos(event.latitude/180*np.pi) #deg
    #roi rounded to 0.0001 (gmt precision for grdcut)
    #get the bounding box of the rupture
    rupt_bb = rupt.surface.get_bounding_box()
    #extend the bounding box by the computed distance to get the shakemap extent
    droi = [rupt_bb.west-dlon,rupt_bb.east+dlon,rupt_bb.south-dlat,rupt_bb.north+dlat]
    #droi = [event.longitude-dlon,event.longitude+dlon,event.latitude-dlat,event.latitude+dlat]
    #droi = [round(v,3) for v in droi]
    #droi = [np.floor((event.longitude-dlon)*10.)/10.,np.ceil((event.longitude+dlon)*10.)/10.,np.floor((event.latitude-dlat)*10.)/10.,np.ceil((event.latitude+dlat)*10.)/10.]

    #check for sites
    if sites==None:
        regular_grid=True
        #get sites form USGS topo
        if roi!=None:
            sites = sml.get_vs30_sites_from_bbox(roi)
            #maybe in addition pgamin
            if pgamin!=None:
                sites = filter_sites(sites,droi)
        else:
            #no sites, no roi
            #slice sites from vs30
            sites = sml.get_vs30_sites_from_bbox(droi)
    else:
        #sites provided
        #if sites == pandas.core.frame.DataFrame:
        #    #convert to dictionary
        #    sites = sites.to_dict('list')
        #elif type(sites) == dict:
        #    pass
        #else:
        #    raise Exception('Unexpected data type for provided sites: {}'.format(type(sites)))
        #keep copy
        sites_params = sites
        #convert to dictionary (lower case strings)
        sites = sites_params.sites
        sites.columns=[s.lower() for s in list(sites.columns)]
        sites = sites.to_dict('list')
        regular_grid = sites.regular_grid
        #TODO: ENSURE UNITS ARE CORRECT!!

        #filter in case
        if roi!=None:
            sites = filter_sites(sites,roi)
        if pgamin!=None:
            #slice sites from vs30
            sites = sml.get_vs30_sites_from_bbox(droi)

    # Runs the shakemap for a given event and site configuration
    shakemap(oqevent, sites)

    #generate a quakemap
    sm = pandas.DataFrame()
    sm["LON"] = sites["lon"]
    sm["LAT"] = sites["lat"]
    for im in imt:
        sm[im] = shakemap["{}/{}/{}/{}".format(eid,gmpe,im,'median')]
        sm["STD"+im] = shakemap["{}/{}/{}/{}".format(eid,gmpe,im,'sigma')]


    #units
    #FIXME: Check how OQ can return unit for GMPE
    units = pandas.DataFrame()
    units["LON"]=['dd']
    units["LAT"]=['dd']
    units["PGA"]=['g']
    units["STDPGA"]=['g']
    for im in imt:
        if im.startswith("SA("):
            units[im]=['g']
            units["STD"+im]=['g']
    units = units.iloc[0]

    #FIXME: Handle somehow (even necessary)!?
    event_specific_uncertainty=pandas.DataFrame()
    event_specific_uncertainty["name"]=['pga', 'pgv', 'mi', 'psa03', 'psa10', 'psa30']
    event_specific_uncertainty["value"]=[0., 0., 0., 0., 0., 0.]
    event_specific_uncertainty["numsta"]=["","","","","",""]

    quakemap=shakeml.quakemap(event,event_specific_uncertainty,sm,units,regular_grid)

    #convert to shakeml format and return
    return shakeml.quakemap2shakeml(quakemap,provider='GFZ')

def main(quakemlfile):
    print(event_shakemap(quakemlfile))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Return shakemap for quakeml event')
    parser.add_argument('quakemlfile', help='Quakemlfile')
    args = parser.parse_args()

    main(args.quakemlfile)




#Get some sites ... here I'm just slicing from the USGS Vs30 topography data, but can be any site model.
#Should be a dictionary of `{"lon": ..., "lat": ..., "vs30": ...}`

# Get some sites
#sites = sml.get_vs30_sites_from_bbox([5., 7., 47., 49.])
#print(sites)
#plt.figure(figsize=(10,9))
#plt.scatter(sites["lon"], sites["lat"], c=sites["vs30"], edgecolor="None")
#plt.colorbar()
#plt.axis('equal')
#plt.xlim(5., 7.)
#plt.ylim(47., 49)
#plt.title("Vs30 m/s")
#plt.show()

# Build an event - e.g. M 5.5
#eq_lon = 6.0
#eq_lat = 48.0
#eq_depth = 7.0
#mag = 5.5
#FIXME:MICHAEL- add rupture plane --> strike,dip,rake
#if not provided i.e. if NAN use hypocentral distance!
#event = sml.Event("EQ001", eq_lon, eq_lat, eq_depth, mag)
#print(event)

## The rupture generation can take place inside the shakemap - or can be input
#rupt = event.get_rupture()
#rlons, rlats = rupt.surface.get_surface_boundaries()
#plt.figure(figsize=(10,9))
#plt.scatter(sites["lon"], sites["lat"], c=sites["vs30"], edgecolor="None")
#plt.plot(rlons[0], rlats[0], 'k-', lw=2)
#plt.plot(eq_lon, eq_lat, "ws")
#plt.colorbar(label="vs30")
#plt.axis('equal')
#plt.xlim(5., 7.)
#plt.ylim(47., 49)
#plt.title("Vs30 m/s")
#plt.show()

# Run the shakemap - for 3 IMTs and 3 GMPEs
#imts = ["PGA", "SA(0.3)", "SA(1.0)"]
#gmpes = ["BindiEtAl2014Rjb", "CauzziEtAl2014", "AkkarEtAlRjb2014"]
# Setup the shakemap tool with the desired GMPEs and - only needs to be done once.
# Note that modification to consider different sets of GMPEs organised by tectonic region is simple
#shakemap = sml.ShakemapLite("dbfile.hdf5", gmpes, imts)
## Deletes any existing database of the same name
#shakemap.reset()
## Runs the shakemap for a given event and site configuration
#shakemap(event, sites)


# To retrieve a particular GMPE and IMT
#bindi_pga = shakemap["EQ001/BindiEtAl2014Rjb/PGA/median"]
#
#plt.figure(figsize=(10,9))
#plt.scatter(sites["lon"], sites["lat"], c=bindi_pga, edgecolor="None", norm=LogNorm(vmin=0.01, vmax=0.3))
#plt.plot(eq_lon, eq_lat, "ws")
#plt.colorbar(label="PGA (g)")
#plt.axis('equal')
#plt.xlim(5., 7.)
#plt.ylim(47., 49)
#plt.show()

