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

#Get some sites ... here I'm just slicing from the USGS Vs30 topography data, but can be any site model.
#Should be a dictionary of `{"lon": ..., "lat": ..., "vs30": ...}`

# Get some sites
sites = sml.get_vs30_sites_from_bbox([5., 7., 47., 49.])
print(sites)
plt.figure(figsize=(10,9))
plt.scatter(sites["lon"], sites["lat"], c=sites["vs30"], edgecolor="None")
plt.colorbar()
plt.axis('equal')
plt.xlim(5., 7.)
plt.ylim(47., 49)
plt.title("Vs30 m/s")
plt.show()

# Build an event - e.g. M 5.5
eq_lon = 6.0
eq_lat = 48.0
eq_depth = 7.0
mag = 5.5
event = sml.Event("EQ001", eq_lon, eq_lat, eq_depth, mag)
print(event)

# The rupture generation can take place inside the shakemap - or can be input
rupt = event.get_rupture()
rlons, rlats = rupt.surface.get_surface_boundaries()
plt.figure(figsize=(10,9))
plt.scatter(sites["lon"], sites["lat"], c=sites["vs30"], edgecolor="None")
plt.plot(rlons[0], rlats[0], 'k-', lw=2)
plt.plot(eq_lon, eq_lat, "ws")
plt.colorbar(label="vs30")
plt.axis('equal')
plt.xlim(5., 7.)
plt.ylim(47., 49)
plt.title("Vs30 m/s")
plt.show()

# Run the shakemap - for 3 IMTs and 3 GMPEs
imts = ["PGA", "SA(0.3)", "SA(1.0)"]
gmpes = ["BindiEtAl2014Rjb", "CauzziEtAl2014", "AkkarEtAlRjb2014"]
# Setup the shakemap tool with the desired GMPEs and - only needs to be done once.
# Note that modification to consider different sets of GMPEs organised by tectonic region is simple
shakemap = sml.ShakemapLite("dbfile.hdf5", gmpes, imts)
# Deletes any existing database of the same name
shakemap.reset()
# Runs the shakemap for a given event and site configuration
shakemap(event, sites)


# To retrieve a particular GMPE and IMT
bindi_pga = shakemap["EQ001/BindiEtAl2014Rjb/PGA/median"]

plt.figure(figsize=(10,9))
plt.scatter(sites["lon"], sites["lat"], c=bindi_pga, edgecolor="None", norm=LogNorm(vmin=0.01, vmax=0.3))
plt.plot(eq_lon, eq_lat, "ws")
plt.colorbar(label="PGA (g)")
plt.axis('equal')
plt.xlim(5., 7.)
plt.ylim(47., 49)
plt.show()

