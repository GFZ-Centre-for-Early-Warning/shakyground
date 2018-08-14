"""
Generates a set of synthetic ruptures based on a point source configuration
"""
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
from math import radians, sin, exp, sqrt, pi, log
from scipy.stats import norm, truncnorm, chi2, multivariate_normal
#import emcee
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.geo import (NodalPlane, PlanarSurface,
                                     Point, Mesh, geodetic)


class PlanarSurfaceAlt(PlanarSurface):
    """
    Alternative to the main planar surface object to calculate ry0 without
    generating the mesh
    """
    def get_ry0_distance(self, mesh):
        """
        :param mesh:
            :class:`~openquake.hazardlib.geo.mesh.Mesh` of points to calculate
            Ry0-distance to.
        :returns:
            Numpy array of distances in km.

        See also :meth:`superclass method <.base.BaseSurface.get_ry0_distance>`
        for spec of input and result values.

        This method uses an average strike direction to compute ry0.
        """
        # This computes ry0 by using an average strike direction
        #top_edge = self.get_mesh()[0:1]
        #strike = self.get_strike()

        dst1 = geodetic.distance_to_arc(self.top_left.longitude,
                                        self.top_left.latitude,
                                        (self.strike + 90.) % 360,
                                        mesh.lons, mesh.lats)

        dst2 = geodetic.distance_to_arc(self.top_right.longitude,
                                        self.top_right.latitude,
                                        (self.strike + 90.) % 360,
                                        mesh.lons, mesh.lats)
        # Find the points on the rupture

        # Get the shortest distance from the two lines
        idx = np.sign(dst1) == np.sign(dst2)
        dst = np.zeros_like(dst1)
        dst[idx] = np.fmin(np.abs(dst1[idx]), np.abs(dst2[idx]))
        return dst


class BaseScaleRel(object):
    """
    More generalised form of the scaling relation class
    """   
    
    def get_rupture_dimensions(self, magnitude, dip, rake, lsd):
        pass

    
class Stafford2014(BaseScaleRel):
    """
    Implements Stafford (2014) scaling relation model
    """
    COEFFS = {
        # Coefficients from Table 1
        1: {"U": {"a0": -27.4922, "a1": 4.6656, "a2": -0.2033},
            "SS": {"a0": -30.8395, "a1": 5.4184, "a2": -0.3044},
            "N": {"a0": -36.9770, "a1": 6.3070, "a2": -0.1696},
            "R": {"a0": -35.8239, "a1": 5.0680, "a2": -0.0457}
           },
        # Coefficients from Table 2
        2: {"U": {"b0": -2.300, "b1": 0.7167, "sigma": 0.2337},
            "SS": {"b0": -2.300, "b1": 0.7167, "sigma": 0.2337},
            "N": {"b0": -4.1055, "b1": 1.0370, "sigma": 0.2509},
            "R": {"b0": -3.8300, "b1": 0.9982, "sigma": 0.2285}
        },
        # Coefficients from Table 3
        3: {"U":{"gamma0": -9.3137, "sigma": 0.3138, "rho": 0.3104},
            "SS": {"gamma0": -9.3137, "sigma": 0.3138, "rho": 0.3104},
            "N": {"gamma0": -9.2483, "sigma": 0.3454, "rho": 0.4336},
            "R": {"gamma0": -9.2749, "sigma": 0.2534, "rho": 0.1376},
           },
        # Coefficients from Table 4
        4: {"U": 0.7574, "SS": 0.7574, "N": 0.8490, "R": 0.7496}
    }
    
    def get_rupture_dimensions(self, magnitude, dip, rake, lsd, epsilon=0.0):
        """
        Gets the rupture dimension from for the given magnitude
        """
        sof = self._get_sof(rake)
        # Get mean and standard deviation of rupture width
        mu_rw, sigma_rw, rw_max, p_i = self.get_rupture_width(magnitude,
                                                              dip, sof, lsd)
        # Get mean and standard deviation of rupture area
        mu_ra, sigma_ra = self.get_rupture_area(magnitude, sof,
                                                rw_max, sigma_rw)
        F_rw_max_norm = norm.cdf(log(rw_max), loc=mu_rw, scale=sigma_rw)
        ncdf_epsilon = norm.cdf(epsilon)
        target = ncdf_epsilon / F_rw_max_norm
        if target > 1:
            target = ncdf_epsilon
        epsilon_rw = norm.ppf(target)
        #print F_rw_max_norm, target, epsilon_rw
        #print p_i, F_rw_max_norm, target, epsilon_rw
        r_w = exp(mu_rw + epsilon_rw * sigma_rw)
        if r_w > rw_max:
            r_w = rw_max
        epsilon_ra = self.COEFFS[4][sof] * epsilon
        r_a = exp(mu_ra + epsilon_ra * sigma_ra)
        r_l = r_a / r_w
        aspect = r_l / r_w
        return r_a, r_l, r_w, aspect,\
            multivariate_normal.pdf(epsilon_rw, epsilon_ra,
                                    cov=self.COEFFS[4][sof])
    
    def _get_sof(self, rake):
        """
        Returns the Style of faulting from the Rake
        """
        if rake is None:
            return "U"
        
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike slip
            return "SS"
        elif rake > 0:
            # thrust/reverse
            return "R"
        else:
            # normal
            return "N"
    
    def get_rupture_width(self, magnitude, dip, sof, lsd):
        """
        Returns the rupture width using 
        """
        # Gets the probability of a full width rupture
        rw_max = lsd / sin(radians(dip))
        z_i = self.COEFFS[1][sof]["a0"] + \
            self.COEFFS[1][sof]["a1"] * magnitude +\
            self.COEFFS[1][sof]["a2"] * rw_max
        p_i = 1.0 / (1.0 + exp(-z_i))
        ln_rw = self.COEFFS[2][sof]["b0"] + \
            self.COEFFS[2][sof]["b1"] * magnitude 
        phi_rw = (log(rw_max) - ln_rw) / self.COEFFS[2][sof]["sigma"]
        phi_rw_ncdf =  norm.cdf(phi_rw)
        ln_rw_trunc = ln_rw - self.COEFFS[2][sof]["sigma"] *\
            (norm.pdf(phi_rw) / phi_rw_ncdf)
        mean_rw = p_i * log(rw_max) + (1.0 - p_i) * ln_rw_trunc
        
        stddev_rw = self._get_rupture_width_sigma(self.COEFFS[2][sof]["sigma"],
                                                  phi_rw,
                                                  phi_rw_ncdf,
                                                  p_i)
        return mean_rw, stddev_rw, rw_max, p_i

    def _get_rupture_width_sigma(self, sigma, phi_rw, phi_rw_ncdf, p_i):
        """
        Returns the variabiliy in the rupture width
        """
        denom = sqrt(2.0 * pi) * phi_rw_ncdf
        if phi_rw_ncdf >= 0.0:
            elem1 = sqrt(pi / 2.) * (1.0 + chi2.cdf(phi_rw, 3))
        else:
            elem1 = sqrt(pi / 2.) * (1.0 - chi2.cdf(phi_rw, 3))
        elem2 = exp(-(phi_rw ** 2)) / denom
        sigma_truncated = sqrt(((sigma ** 2.) / denom) * (elem1 - elem2))
        return (1.0 - p_i) * sigma_truncated

    def get_rupture_area(self, magnitude, sof, rw_max, sigma_lnrw):
        """
        Returns the rupture area
        """
        mw_crit = (log(rw_max) - self.COEFFS[2][sof]["b0"]) /\
            self.COEFFS[2][sof]["b1"]
        ln_ra = self.COEFFS[3][sof]["gamma0"] + log(10.) * magnitude
        if magnitude > mw_crit:
            ln_ra -= ((log(10.) / 4.) * (magnitude - mw_crit))
        # Get the sigma log rupture area
        sigma = self.COEFFS[3][sof]["sigma"]
        sigma_ra = sqrt(
            (sigma ** 2.) + (sigma_lnrw ** 2.) +
            (2.0 * self.COEFFS[3][sof]["rho"] * sigma_lnrw * sigma))
        return ln_ra, sigma_ra


NGAWEST_HYPO_DISTRIBUTIONS = {
    "All": {"Mean": np.array([0.46851357, 0.67460261]),
            "Sigma": np.array([0.20691535,0.22040932]),
            "COV": np.array([[0.04293594, 0.00103987],
                             [0.00103987, 0.04871868]]),
            "DipRange": [30.0, 90.0]},
    "SS": {"Mean": np.array([0.48545206, 0.64942746]),
           "Sigma": np.array([ 0.22415657,  0.21677068]),
           "COV": np.array([[0.05064814, -0.00603003],
                            [-0.00603003,  0.04736544]]),
           "DipRange": [80.0, 90.0]},
    "R": {"Mean": np.array([0.4674859, 0.58483914]),
          "Sigma": np.array([ 0.16275562, 0.22017015]),
          "COV": np.array([[0.02673021, 0.0113362],
                           [0.0113362, 0.04891558]]),
          "DipRange": [25.0, 50.0]},
    "N": {"Mean": np.array([0.50887254, 0.82404]),
          "Sigma": np.array([0.22416128, 0.13647917]),
          "COV": np.array([[0.05085368, -0.00332741],
                           [-0.00332741, 0.01885098]]),
          "DipRange": [50.0, 80.0]}
}


NGAWEST_HYPO_DISTRIBUTIONS["TF"] = NGAWEST_HYPO_DISTRIBUTIONS["R"]
NGAWEST_HYPO_DISTRIBUTIONS["U"] = NGAWEST_HYPO_DISTRIBUTIONS["All"]


def lnprob_trunc_norm(x, mean, bounds, C):
    if np.any(x < bounds[:,0]) or np.any(x > bounds[:,1]):
        return -np.inf
    else:
        return -0.5*(x-mean).dot(np.linalg.inv(C)).dot(x-mean)


class FiniteRuptureSampler(object):
    """
    Module to produce a set of finite rupture samples for an event
    Process:
    1. From scaling relation and seismogenic thicknesses sample rupture
       parameters (area, length, width)
    2. Sample hypocentre location from multivariate Gaussian distribution
    3. a. If mechanism is known, use strike and dip
       b. If mechanism is not known sample strike and dip from distribution
    4. For each sample calculate Rjb, Rrup, Rx and Ry0
    5. For each distance calculate the median value
    6. For each rupture define a penalty function 
    
    P = \sum_{i=rjb, rrup, ry0, rx} (Ri - median(Ri)) ** 2.
    Take the rupture with the lowest penalty function!
    """
    def __init__(self, scalerel, usd, lsd):
        self.msr = scalerel
        self.usd = usd
        self.lsd = lsd
        self.thickness = self.lsd - self.usd
    
    def get_finite_rupture_distances(self, nsamples, centroid, magnitude,
                                     site_lons, site_lats, mechanisms,
                                     dimensions=None, rjb=[], rrup=[],
                                     rdim=0.0, weights=[0.25, 0.25, 0.25, 0.25]):
        """
        Gets the finite rupture distances
        """
        # Sample the ruptures
        ruptures = self.sample_ruptures(centroid, magnitude, nsamples,
                                        mechanisms, dimensions)
        # Calculate the distances
        distances, rhypo = self.get_distances(centroid, ruptures,
                                              site_lons, site_lats)
        # Get the ruptures distances and central ruptures
        if len(rjb) or len(rrup):
            # Has Rrup or Rjb pre-calculated, find the best fitting rupture
            central_rupture, central_distances =\
                self.get_best_matching_rupture(rjb, rrup, ruptures, distances,
                                               rhypo, rdim)
            for i, row in enumerate(central_distances):
                if not len(rjb):
                    comp_rjb = ""
                else:
                    comp_rjb = "{:.4f}".format(rjb[i])
                if not len(rrup):
                    comp_rrup = ""
                else:
                    comp_rrup = "{:.4f}".format(rrup[i])
                print("%.5fE, %.5fN: Rjb = %.4f (orig. %s), "
                      "Rrup = %.4f (orig. %s)" % (site_lons[i], site_lats[i],
                                                  row[0],
                                                  comp_rjb,
                                                  row[1],
                                                  comp_rrup))
                    
        else:
            # Get the most central rupture
            central_rupture, central_distances = self.get_central_rupture(
                ruptures, distances, rhypo, rdim, weights
                )
        return central_rupture, central_distances, np.mean(distances, axis=0)
    
    def sample_ruptures(self, centroid, magnitude, nsamples, mechanisms="U",
                        dimensions=None):
        """
        Mechanisms should be a list of Nodal plane objects
        """
        strikes, dips, rakes = self.sample_mechanisms(nsamples, mechanisms)
        # If dimensions
        if dimensions:
            areas = dimensions["area"] * np.ones(nsamples)
            lengths = dimensions["length"] * np.ones(nsamples)
            widths = dimensions["width"] * np.ones(nsamples)
        else:
            # Get the physical rupture dimensions by sampling from the
            # Stafford (2014) model
            areas, lengths, widths = self.sample_rupture_dimensions(
                nsamples, magnitude, strikes, dips, rakes)
        sofs = [self.msr._get_sof(rake) for rake in rakes]
        # Get hypocentres
        hypo_locs = self.sample_hypocentres(nsamples, sofs)
        # Build the ruptures
        return self.build_ruptures(centroid, strikes, dips,
                                   lengths, widths, hypo_locs)
    
    def sample_rupture_dimensions(self, nsamples, magnitude, strikes,
                                  dips, rakes):
        """
        Samples the rupture dimension from the magnitude scaling relation
        """
        # Sample epsilon values
        msr_epsilons = truncnorm.rvs(-3.0, 3.0,
                                     loc=0, scale=1.0,
                                     size=nsamples)
        # Generate rupture dimensions
        lengths = np.zeros(nsamples)
        widths = np.zeros(nsamples)
        areas = np.zeros(nsamples)
        for i, (dip, rake, epsilon) in enumerate(
            zip(dips, rakes, msr_epsilons)):
            areas[i], lengths[i], widths[i], _, _ = \
                self.msr.get_rupture_dimensions(magnitude, dip, rake,
                                                self.thickness, epsilon)
        #print np.column_stack([msr_epsilons, lengths, widths, areas])
        return areas, lengths, widths
    
    def sample_mechanisms(self, nsamples, mechanisms):
        """
        Samples a set of mechanisms from either a PMF or a style of
        faulting string
        """
        # Sample rupture dimensions
        if isinstance(mechanisms, PMF):
            # Sample mechanisms
            mechanism_samples = mechanisms.sample(nsamples)
            strikes, dips, rakes = zip(*[(mech.strike, mech.dip, mech.rake)
                                         for mech in mechanism_samples])
            strikes = np.array(strikes)
            dips = np.array(dips)
            rakes = np.array(rakes)
        else:
            if not mechanisms in NGAWEST_HYPO_DISTRIBUTIONS:
                mechanisms = "All"
            # Sample strike dip and rake freely
            strikes = np.random.uniform(0., 360., nsamples)
            dips = np.random.uniform(
                NGAWEST_HYPO_DISTRIBUTIONS[mechanisms]["DipRange"][0],
                NGAWEST_HYPO_DISTRIBUTIONS[mechanisms]["DipRange"][1],
                nsamples)
            rakes = [None for i in range(nsamples)]
        return strikes, dips, rakes
    
    def sample_hypocentres(self, nsamples, sofs):
        """
        Samples the hypocentral location from 
        """
        # Count the unique styles of faulting
        sof_counter = {}
        for key in NGAWEST_HYPO_DISTRIBUTIONS:
            sof_counter[key] = sofs.count(key)
        hypo_locs = np.zeros([nsamples, 2])
        for key in sof_counter:
            if not sof_counter[key]:
                continue
            mean = NGAWEST_HYPO_DISTRIBUTIONS[key]["Mean"]
            C = NGAWEST_HYPO_DISTRIBUTIONS[key]["COV"]
            pos = np.random.multivariate_normal(mean, C, sof_counter[key])
            all_valid = np.logical_and(
                np.logical_and(pos[:, 0] >= 0., pos[:, 0] <= 1.),
                np.logical_and(pos[:, 1] >= 0., pos[:, 1] <= 1.)
                )
            while np.sum(all_valid) < sof_counter[key]:
                idx = np.logical_not(all_valid)
                tpos = np.random.multivariate_normal(mean, C, sum(idx))
                pos[idx] = tpos
                all_valid = np.logical_and(
                    np.logical_and(pos[:, 0] >= 0., pos[:, 0] <= 1.),
                    np.logical_and(pos[:, 1] >= 0., pos[:, 1] <= 1.)
                )
            cntr = 0
            for i in range(sof_counter[key]):
                if sofs[i] == key:
                    hypo_locs[i, :] = pos[cntr, :]
                    cntr += 1
        return hypo_locs
            
    def build_ruptures(self, centroid, strikes, dips,
                       lengths, widths, hypo_locs):
        """
        Builds the rupture set
        """
        ruptures = []
        rstrikes = strikes * np.pi / 180.
        rdips = dips * np.pi / 180.
        for (strike, rstrike, dip, rdip, length, width, hypo_loc) in zip(
            strikes, rstrikes, dips, rdips, lengths, widths, hypo_locs):
            # Deal with the dip part first
            area = length * width
            updip_width = hypo_loc[1] * width
            downdip_width = (1.0 - hypo_loc[1]) * width
            updip_depth_change = updip_width * np.sin(rdip)
            downdip_depth_change = downdip_width * np.sin(rdip)
            if centroid.depth < updip_depth_change:
                # This would move the rupture above the top surface
                # readjust
                offset = updip_depth_change - centroid.depth
                updip_depth_change = centroid.depth
                downdip_depth_change += offset
            if downdip_depth_change > (self.lsd - centroid.depth):
                if (updip_depth_change + downdip_depth_change) >\
                    self.thickness:
                    # Determine excess width
                    offset = (centroid.depth + downdip_depth_change) - self.lsd
                    offset_area = length * (offset / np.sin(rdip))
                    rw_max = self.thickness / np.sin(rdip)
                    length += (offset_area / rw_max)
                    updip_depth_change = centroid.depth
                    downdip_depth_change = self.lsd - centroid.depth
                else:
                    # This would move the rupture below the lower surface
                    offset = (centroid.depth + downdip_depth_change) - self.lsd
                    downdip_depth_change = self.lsd - centroid.depth
                    updip_depth_change += offset
            updip_surface_length = updip_depth_change / np.tan(rdip)
            downdip_surface_length = downdip_depth_change / np.tan(rdip)
            # Now deal with strike parts
            left_length = hypo_loc[0] * length
            right_length = (1.0 - hypo_loc[0]) * length
            # Build corner points
            downdip_dir = (dip + 90.0) % 360
            updip_dir = (dip - 90.0) % 360
            mid_left = centroid.point_at(left_length,
                                         0.0,
                                         (strike + 180.) % 360.)
            mid_right = centroid.point_at(right_length, 0.0, strike)
            top_left = mid_left.point_at(updip_surface_length,
                                         -updip_depth_change,
                                         updip_dir)
            top_right = mid_right.point_at(updip_surface_length,
                                           -updip_depth_change,
                                           updip_dir)
            bottom_left = mid_left.point_at(downdip_surface_length,
                                            downdip_depth_change,
                                            downdip_dir)
            bottom_right = mid_right.point_at(downdip_surface_length,
                                              downdip_depth_change,
                                              downdip_dir)
            try:
                surface = PlanarSurface.from_corner_points(top_left,
                                                           top_right,
                                                           bottom_right,
                                                           bottom_left)
            except:
                print(strike, dip, length, width, hypo_loc)
                print(top_left, top_right, bottom_right, bottom_left)
                raise ValueError
            ruptures.append(surface)
        return ruptures
    
    def get_distances(self, centroid, ruptures, lons, lats):
        """
        Calculates the source to site distances
        """
        nsites = len(lons)
        mesh = Mesh(lons, lats)
        distances = np.zeros([len(ruptures), nsites, 4])
        for i, rup in enumerate(ruptures):
            distances[i, :, 0] = rup.get_joyner_boore_distance(mesh)
            distances[i, :, 1] = rup.get_min_distance(mesh)
            distances[i, :, 2] = rup.get_rx_distance(mesh)
            distances[i, :, 3] = rup.get_ry0_distance(mesh)
        rhypo = centroid.distance_to_mesh(mesh)
        return distances, rhypo
    
    def get_central_rupture(self, ruptures, distances,
                            hypocentral_distances,
                            rdim=0, weights=[0.25, 0.25, 0.25, 0.25]):
        """
        Returns the rupture closest to the median distances
        """
        nsites = distances.shape[1]
        site_weights = 1./ (hypocentral_distances ** rdim)
        site_weights = site_weights / np.sum(site_weights)

        median_distances = np.zeros([4, nsites])
        for i in range(4):
            median_distances[i, :] = np.percentile(distances[:, :, i],
                                                   50, axis=0)
        penalty_function = np.zeros(distances.shape[0])
        for i in range(distances.shape[2]):
            site_penalty = np.zeros(distances.shape[0])
            for k in range(distances.shape[1]):
                site_penalty += (site_weights[k] * (
                    np.sqrt((distances[:, k, i] - median_distances[i, k]) ** 2.)
                    ))
            penalty_function += (weights[i] * site_penalty)
        min_loc = np.argmin(penalty_function)
        return ruptures[min_loc], distances[min_loc, :, :]

    def get_best_matching_rupture(self, rjb, rrup, ruptures, distances,
                                  hypocentral_distances, rdim=0.0):
        """
        If a Joyner-Boore distance and/or rupture distance is given find the
        rupture that matches these most closely
        """
        nsites = distances.shape[1]
        site_weights = 1./ (hypocentral_distances ** rdim)
        site_weights = site_weights / np.sum(site_weights)

        if len(rjb) and not len(rrup):
            # Match to Rjb only
            assert len(rjb) == distances.shape[1]
            idx = np.array([0])
            weights = np.array([1.0])
            target_distances = np.zeros([1, nsites])
            target_distances[0, :] = rjb
        elif len(rrup) and not len(rjb):
            # Match to Rrup only
            assert len(rrup) == distances.shape[1]
            idx = np.array([1])
            weights = np.array([1.0])
            target_distances = np.zeros([1, nsites])
            target_distances[0, :] = rrup
        else:
            # Match to Rjb and Rrup evenly
            assert len(rrup) == len(rjb)
            assert len(rrup) == distances.shape[1]
            idx = np.array([0, 1])
            weights = np.array([0.5, 0.5])
            target_distances = np.zeros([2, nsites])
            target_distances[0, :] = rjb
            target_distances[1, :] = rrup
        penalty_function = np.zeros(distances.shape[0])
        for i in range(len(idx)):
            site_penalty = np.zeros(distances.shape[0])
            for k in range(distances.shape[1]):
                site_penalty += (site_weights[k] * 
                    np.sqrt((distances[:, k, idx[i]] -
                             target_distances[idx[i], k]) ** 2.0))
            penalty_function += (weights[i] * site_penalty)
        min_loc = np.argmin(penalty_function)
        return ruptures[min_loc], distances[min_loc, :, :]

