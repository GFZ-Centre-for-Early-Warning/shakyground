"""
Prototype of a "Shakemap-Lite" system for retreiving expected GMPE
values from an earthquake
"""
import os
import subprocess
import h5py
import numpy as np
import pandas as pd
from openquake.hazardlib.geo import Point, PlanarSurface, NodalPlane, MultiSurface, Mesh
from openquake.hazardlib.gsim import get_available_gsims
from openquake.hazardlib.imt import PGA, PGV, SA, from_string
from openquake.hazardlib.contexts import (ContextMaker, get_distances,
                                          SitesContext, DistancesContext)
from openquake.hazardlib.source.rupture import ParametricProbabilisticRupture
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.scalerel.wc1994 import WC1994
from openquake.hazardlib import const
from openquake.hazardlib.pmf import PMF
import synthetic_rupture_generator as srg


GSIM_SET = get_available_gsims()


def create_planar_surface(top_centroid, strike, dip, area, aspect):
    """
    Given a central location, create a simple planar rupture
    :param top_centroid:
        Centroid of trace of the rupture, as instance of :class:
            openquake.hazardlib.geo.point.Point
    :param float strike:
        Strike of rupture(Degrees)
    :param float dip:
        Dip of rupture (degrees)
    :param float area:
        Area of rupture (km^2)
    :param float aspect:
        Aspect ratio of rupture

    :returns: Rupture as an instance of the :class:
        openquake.hazardlib.geo.surface.planar.PlanarSurface
    """
    rad_dip = dip * pi / 180.
    width = sqrt(area / aspect)
    length = aspect * width
    # Get end points by moving the top_centroid along strike
    top_right = top_centroid.point_at(length / 2., 0., strike)
    top_left = top_centroid.point_at(length / 2.,
                                     0.,
                                     (strike + 180.) % 360.)
    # Along surface width
    surface_width = width * cos(rad_dip)
    vertical_depth = width * sin(rad_dip)
    dip_direction = (strike + 90.) % 360.

    bottom_right = top_right.point_at(surface_width,
                                      vertical_depth,
                                      dip_direction)
    bottom_left = top_left.point_at(surface_width,
                                    vertical_depth,
                                    dip_direction)

    # Create the rupture
    return PlanarSurface(strike, dip, top_left, top_right,
                         bottom_right, bottom_left)


def vs30_to_z1pt0_cy14(vs30, japan=False):
    """
    Returns the estimate depth to the 1.0 km/s velocity layer based on Vs30
    from Chiou & Youngs (2014) California model

    :param numpy.ndarray vs30:
        Input Vs30 values in m/s
    :param bool japan:
        If true returns the Japan model, otherwise the California model
    :returns:
        Z1.0 in m
    """
    if japan:
        c1 = 412. ** 2.
        c2 = 1360.0 ** 2.
        return np.exp((-5.23 / 2.0) *
                      np.log((np.power(vs30,2.) + c1) / (c2 + c1)))
    else:
        c1 = 571 ** 4.
        c2 = 1360.0 ** 4.
        return np.exp((-7.15 / 4.0) * np.log((vs30 ** 4. + c1) / (c2 + c1)))


def vs30_to_z2pt5_cb14(vs30, japan=False):
    """
    Converts vs30 to depth to 2.5 km/s interface using model proposed by
    Campbell & Bozorgnia (2014)

    :param vs30:
        Vs30 values (numpy array or float)

    :param bool japan:
        Use Japan formula (True) or California formula (False)

    :returns:
        Z2.5 in km
    """
    if japan:
        return np.exp(5.359 - 1.102 * np.log(vs30))
    else:
        return np.exp(7.089 - 1.144 * np.log(vs30))


def get_vs30_sites_from_bbox(bbox, isep="\t"):
    """
    Returns a basic site dictionary from a bbox [llon, ulon, llat, ulat]
    """
    filepath=os.path.dirname(__file__)
    site_data_path = os.path.join(filepath, "global_vs30.grd")
    tempfile = "tempfile.grd"
    # Call grdcut
    cutstring = "/".join([str(loc) for loc in bbox])
    subprocess.run(["gmt", "grdcut", site_data_path,
                    "-G{:s}".format(tempfile),
                    "-R{:s}".format(cutstring)])
    # Call grd2xyz
    tempxyz = os.path.join(filepath,"tempfile.xyz")
    subprocess.run(["gmt", "grd2xyz", tempfile, "-sa", ">", tempxyz])
    # Use pandas to read in the xyzdata
    site_data = pd.read_csv(tempxyz, sep=isep)
    sites = {"lon": site_data.iloc[:, 0].values,
             "lat": site_data.iloc[:, 1].values,
             "vs30": site_data.iloc[:, 2].values}
    os.remove(tempfile)
    os.remove(tempxyz)
    return sites


class Event(object):
    """
    Shakemap event object - requires a minimum of an ID, longitude, latitude,
    depth and magnitude. Other attributes con be input to control the rupture
    orientation and mechanism.

    Can input a rupture geometry directly
    """
    def __init__(self, i_d, lon, lat, hypo_depth, mag, strike=0.0, dip=90.,
                 rake=0.0, aspect=1.0, msr=srg.Stafford2014(), usd=0.0,
                 lsd=1000.0, rupture=None):

        self.id = i_d
        self.lon = lon
        self.lat = lat
        self.depth = hypo_depth
        self.mag = mag
        self.strike = strike
        self.dip = dip
        self.rake = rake
        #self.mechanism = []
        self.mechanism = [({"strike": strike, "dip": dip, "rake": rake}, 1.)]
            #({"strike": strike, "dip": dip, "rake": rake}, 1.)]
        self.aspect = aspect
        self.msr = msr
        if rupture:
            assert isinstance(rupture, ParametricProbabilisticRupture)
            self.rupture = rupture
        else:
            self.rupture = None
        self.usd = usd
        self.lsd = lsd
        self.generator = srg.FiniteRuptureSampler(msr, self.usd, self.lsd)


    def __repr__(self):
        return "{:s}|{:.5f}|{:.5f}|{:.5f}|{:.2f}".format(self.id,
                                                         self.lon,
                                                         self.lat,
                                                         self.depth,
                                                         self.mag)
    def get_rupture(self):
        """
        If a rupture is provided then it is returned, otherwise it will
        use a synthetic rupture generator
        """
        if self.rupture:
           return rupture
        else:
            if len(self.mechanism) >= 1:
                # Define a set of rupture mechanisms with weight
                mechanisms = []
                for mech, weight in self.mechanism:
                    npd = NodalPlane(mech["strike"],
                                     mech["dip"],
                                     mech["rake"])
                    mechanisms.append((weight, npd))
                mechanisms = PMF(mechanisms)
            else:
                if self.rake:
                    if self.rake >= 45.0 and self.rake <= 135.0:
                        mechanisms = "R"
                    elif self.rake >= -135.0 and self.rake <= -45.:
                        mechanisms = "N"
                    else:
                        mechanisms = "SS"
                else:
                    mechanisms = "U"
            # Build planar rupture
            planar_surface = self.generator.sample_ruptures(
                Point(self.lon, self.lat, self.depth),
                      self.mag, nsamples=1, mechanisms=mechanisms,
                      dimensions=None)[0]
            return ParametricProbabilisticRupture(
                self.mag, self.rake, None,
                Point(self.lon, self.lat, self.depth),
                planar_surface, 1.0, None)


class ShakemapLite(object):
    """
    Class to implement a lightweight OpenQuake-based shakemap
    """
    def __init__(self, database, gmpes, imts):
        """
        Shakemap results are stored to an hdf5 databse
        :param str database:
            Path to database
        :param list gmpes:
            List of gmpes (as strings)
        :param list imts:
            List of IMTs (as strings)
        """

        self.gmpes = [GSIM_SET[gmpe]() for gmpe in gmpes]
        self.imts = [from_string(imt) for imt in imts]
        self.db_file = database
        self.context = ContextMaker(self.gmpes)
        # Determine the site attributes
        self.site_attribs = []
        for gmpe in self.gmpes:
            for site_attrib in gmpe.REQUIRES_SITES_PARAMETERS:
                if not site_attrib in self.site_attribs:
                    self.site_attribs.append(site_attrib)

    def __call__(self, event, sites, vs30measured=False, backarc=False):
        """
        Execute a shakemap
        """
        # Build the contexts
        sctx, dctx = self.build_sites_distances_contexts(event, sites,
                                                         vs30measured,
                                                         backarc)
        # Connect to database
        fle = h5py.File(self.db_file)
        grp = fle.create_group(event.id)
        # Store distance calculations to database
        dists_grp = grp.create_group("Distances")
        for key, distance in list(dctx.__dict__.items()):
            dset = dists_grp.create_dataset(key, distance.shape, dtype="f")
            dset[:] = distance
        # Get the rupture
        rctx = event.get_rupture()
        # Run the calculations
        for gmpe in self.gmpes:
            gmpe_grp = grp.create_group(gmpe.__class__.__name__)
            for imt in self.imts:
                imt_grp = gmpe_grp.create_group(str(imt))
                # Get ground motion values
                mean, [sigma] = gmpe.get_mean_and_stddevs(
                    sctx, rctx, dctx, imt, [const.StdDev.TOTAL])
                mean_dset = imt_grp.create_dataset("median",
                                                   mean.shape,
                                                   dtype="f")
                mean_dset[:] = np.exp(mean)
                sigma_dset = imt_grp.create_dataset("sigma",
                                                    sigma.shape,
                                                    dtype="f")
                sigma_dset[:] = sigma
        fle.close()

    def reset(self):
        """
        Remove the database file then delete the object itself
        """
        if os.path.exists(self.db_file):
            os.remove(self.db_file)

    def __getitem__(self, i):
        """
        Returns a given data set on a path
        """
        if not os.path.exists(self.db_file):
            raise AttributeError("Database %s does not exist" % self.db_file)
        # Connect to the database
        fle = h5py.File(self.db_file)
        try:
            data = fle[i][:]
            return data
        except:
            # Data not found
            fle.close()
            raise AttributeError("%s not found in database" % i)

    def build_sites_distances_contexts(self, event, sites, vs30measured=False,
                                       backarc=False):
        """
        Builds the contexts from the event and sites
        """
        # Get distances
        mesh = Mesh(sites["lon"], sites["lat"])
        mshape = mesh.lons.shape
        dctx = DistancesContext()
        for param in self.context.REQUIRES_DISTANCES:
            setattr(dctx, param, get_distances(event.get_rupture(),
                                               mesh, param))
        # Get sites context
        sctx = SitesContext()
        for key in self.site_attribs:
            if key.startswith("z1pt0") and not "z1pt0" in sites:
                setattr(sctx, "z1pt0",
                        vs30_to_z1pt0_cy14(sites["vs30"], japan=False))

            elif key.startswith("z2pt5") and not "z5pt5" in sites:
                setattr(sctx, "z2pt5",
                        vs30_to_z2pt5_cb14(sites["vs30"], japan=False))

            elif key.startswith("vs30measured") and not "vs30measured" in sites:
                if vs30measured:
                    setattr(sctx, "vs30measured", np.ones(mshape, dtype=bool))
                else:
                    setattr(sctx, "vs30measured", np.zeros(mshape, dtype=bool))
            elif key.startswith("vs30measured") and not "vs30measured" in sites:
                if vs30measured:
                    setattr(sctx, "vs30measured", np.ones(mshape, dtype=bool))
                else:
                    setattr(sctx, "vs30measured", np.zeros(mshape, dtype=bool))
            elif key.startswith("backarc") and not "backarc" in sites:
                if vs30measured:
                    setattr(sctx, "backarc", np.ones(mshape, dtype=bool))
                else:
                    setattr(sctx, "backarc", np.zeros(mshape, dtype=bool))
            else:
                setattr(sctx, key, sites[key])
        return sctx, dctx

#
#    def build_site_collection(self, sites, vs30measured=False, backarc=False):
#        """
#        Sites should be a dictionary of site parameter
#        {lon, lat, vs30, vs30measured, z1pt0, z2pt5, backarc}
#        """
#        for key in ["lon", "lat", "vs30"]:
#            # check in dict
#            assert key in sites
#        nsites = len(sites["lon"])
#        if not "z1pt0" in sites or sites["z1pt0"] is None:
#            sites["z1pt0"] = vs30_to_z1pt0_cy14(sites["vs30"], japan=False)
#
#        if not "z2pt5" in sites or sites["z2pt5"] is None:
#            sites["z2pt5"] = vs30_to_z2pt5_cb14(sites["vs30"], japan=False)
#
#        if not "vs30measured" in sites or sites["vs30measured"] is None:
#            if vs30measured:
#                sites["vs30measured"] = np.ones(nsites, dtype=bool)
#            else:
#                sites["vs30measured"] = np.zeros(nsites, dtype=bool)
#
#        if not "backarc" in sites or sites["backarc"] is None:
#            sites["backarc"] = backarc
#        site_col = []
#        for i in range(len(sites["lon"])):
#            site = Site(Point(sites["lon"][i], sites["lat"][i], 0.0),
#                        sites["vs30"][i], sites["vs30measured"][i],
#                        sites["z1pt0"][i], sites["z2pt5"][i],
#                        sites["backarc"][i])
#            site_col.append(site)
#        return SiteCollection(site_col)
