#! /usr/bin/env python
# -*- coding: utf-8 -*-

# import os
# import pandas as pd
import sep
import numpy as np
import matplotlib.pyplot as plt

from astropy import units
# from pymage import panstarrs
from matplotlib.patches import Ellipse
from astrobject.baseobject import get_target
from astrobject.utils.tools import flux_to_mag
from astrobject.photometry import get_photopoint
from astrobject.collections.photodiagnostics import get_massestimator

# PanSTARRS zero point according to
# https://coolwiki.ipac.caltech.edu/index.php/
# Central_wavelengths_and_zero_points
ps1_zp = {'g': 0.4810,
          'r': 0.6170,
          'i': 0.7520,
          'z': 0.8660,
          'y': 0.9620}

# From astrobject/instruments/panstarrs.py
PANSTARRS_INFO = {'g': {'lbda': 4866.457871, 'ABmag0': 25.0, 'band': 'ps1.g'},
                  'r': {'lbda': 6214.623038, 'ABmag0': 25.0, 'band': 'ps1.r'},
                  'i': {'lbda': 7544.570357, 'ABmag0': 25.0, 'band': 'ps1.i'},
                  'z': {'lbda': 8679.482571, 'ABmag0': 25.0, 'band': 'ps1.z'},
                  'y': {'lbda': 9633.284241, 'ABmag0': 25.0, 'band': 'ps1.y'}}

# Pymage
try:
    from pymage import panstarrs
    _HAS_PYMAGE = True
except ImportError:
    _HAS_PYMAGE = False


# =========================================================================== #
#                                                                             #
#                                 MASS  CLASS                                 #
#                                                                             #
# =========================================================================== #

class MassMeasure(object):
    '''
    Usage:
    mm = MassMeasure(ra, dec, z, cutout=None)

    mm.get_hostgalaxy_params()
    '''

    # =================================================================== #
    #                               Initial                               #
    # =================================================================== #

    def __init__(self, ra, dec, z, cutout=None, **kwargs):
        '''
        Sets the object's ra, dec coordinates and redshift,
        then instantiates a PS1Target object from `pymage.panstarrs`,
        downloading the needed cutout if not given as an entry parameter.
        '''
        if not _HAS_PYMAGE:
            raise ImportError('This method needs pymage: `pip install pymage`')

        self.target = panstarrs.PS1Target.from_coord(ra, dec)
        self.ra = self.target.coordinate.ra.deg
        self.dec = self.target.coordinate.dec.deg
        self.z = z
        if cutout is None:
            self.target.download_cutout(load_weight=True, **kwargs)
        else:
            self.target._cutout = cutout
        self._cutout = self.target.imgcutout
        self.bands = list(self.cutout.keys())
        sep_obj = self.target.imgcutout['r'].sep_extract(returnobjects=True)
        self.ellipses = sep_obj.get_ellipse_values().T
        self.un_nan(self)

    def update(self, **kwargs):
        '''To redownload a cutout if size was incorrect'''
        self.target.download_cutout(load_weight=True, **kwargs)
        self._cutout = self.target.imgcutout
        sep_obj = self.target.imgcutout['r'].sep_extract(returnobjects=True)
        self.ellipses = sep_obj.get_ellipse_values().T
        self.get_hostgalaxy()
        self.get_hostgalaxy_params()

    def un_nan(self):
        '''Sets the NaNs data to 0 and NaNs var to 1'''
        for band in self.bands:
            mask = self.cutout[band].data == self.cutout[band].data
            self.cutout[band].data[~mask] = 0
            mask = self.cutout[band].weightimage.data\
                == self.cutout[band].weightimage.data
            self.cutout[band].weightimage.data[~mask] = 1

    # =================================================================== #
    #                               Methods                               #
    # =================================================================== #

    # ------------------------------------------------------------------- #
    #                               EXTFUNC                               #
    # ------------------------------------------------------------------- #

    @staticmethod
    def get_DLR(x_sn, y_sn, x_gal, y_gal, a, b, theta):
        '''
        Computes the DLR from the target position for a galaxy defined
        by x, y coordinates, `a` and `b` semimajor and semiminor axes
        '''
        cxx = np.cos(theta)**2/a**2 + np.sin(theta)**2/b**2
        cyy = np.sin(theta)**2/a**2 + np.cos(theta)**2/b**2
        cxy = 2*np.sin(theta)*np.cos(theta)*(1/a**2 + 1/b**2)
        r_gal = np.sqrt(
            cxx*(x_sn - x_gal)**2 +
            cyy*(y_sn - y_gal)**2 +
            cxy*(x_sn - x_gal)*(y_sn - y_gal)
        )
        return r_gal

    def count_to_flux(self, count):
        '''Gives a flux from given count'''
        if not hasattr(self, '_cutout'):
            raise AttributeError(
                "No cutout loaded yet. " +
                "Run self.download() or set self._cutout")
        else:
            return self.target.imgcutout['r'].count_to_flux(count)

    # ------------------------------------------------------------------- #
    #                               SETTER                                #
    # ------------------------------------------------------------------- #

    def set_hostgalaxy(self, ind):
        '''
        If something goes wrong and you want to choose another HG than
        the one determined by the DLR using self.get_hostgalaxy()
        '''
        self.hg_ellipse = self.ellipses[ind]

    def set_septhresh(self, thresh):
        '''
        If something goes wrong with the sep_extract and the threshold needs to
        be changed
        '''
        self.target.sep_extract(thresh=thresh)
        self.target.set_sep()
        self.ellipses = self.target._sep.sepobjects.get_ellipse_values().T

    # ------------------------------------------------------------------- #
    #                               GETTER                                #
    # ------------------------------------------------------------------- #

    def get_hostgalaxy(self):
        '''Automatically computes each ellipse's DLR to the SN location, then
        gives x, y, a, b, theta of host galaxy determined to be the object with
        smallest DLR'''
        self.x_sn, self.y_sn = np.array(np.shape(self.cutout['r'].data))/2
        self.DLR_list = [self.get_DLR(self.x_sn, self.y_sn,
                                      self.ellipses[i][0], self.ellipses[i][1],
                                      self.ellipses[i][2], self.ellipses[i][3],
                                      self.ellipses[i][4])
                         for i in range(len(self.ellipses))]
        self.DLR_list = np.asarray(self.DLR_list)
        self.minDLR = np.min(self.DLR_list[self.ellipses[:, 2] > 3])
        hg_ind = np.where(self.DLR_list == self.minDLR)[0]
        self.hg_ellipse = self.ellipses[hg_ind][0]

    def get_hostgalaxy_params(self, scaleup=2.5):
        '''
        Gives the count, flux, mag (not used) and then mass of the selected
        host.
        '''
        if not hasattr(self, "hg_ellipse"):
            raise AttributeError(
                "No host galaxy detected yet. " +
                "Run self.get_hostgalaxy() or set self.hg_ellipse")
        self.var_count = {band:
                          self.target.imgcutout[band].weightimage.data**(-2)
                          for band in self.bands}

        counts_res = {band: sep.sum_ellipse(self.cutout[band].data,
                                            self.hg_ellipse[0],
                                            self.hg_ellipse[1],
                                            2*self.hg_ellipse[2]*scaleup,
                                            2*self.hg_ellipse[3]*scaleup,
                                            self.hg_ellipse[4],
                                            var=self.var_count[band])
                      for band in self.bands}
        self.hg_count = {band: [counts_res[band][0], counts_res[band][1]]
                         for band in self.bands}

        self.hg_flux = {band: [self.count_to_flux(self.hg_count[band][0]),
                               self.count_to_flux(self.hg_count[band][1])]
                        for band in self.bands}

        self.var_flux = {band: self.hg_flux[band][-1]**2
                         for band in self.bands}

        self.hg_mag = {band: np.array(flux_to_mag(self.hg_flux[band][0],
                                                  self.hg_flux[band][1],
                                                  zp=ps1_zp[band]))
                       for band in self.bands}

        phtpt = {band: get_photopoint(self.hg_flux[band][0],
                                      self.var_flux[band],
                                      lbda=PANSTARRS_INFO[band]['lbda'],
                                      bandname=PANSTARRS_INFO[band]['band'])
                 for band in ['g', 'i']}
        MassEstimate_obj = get_massestimator(photopoints=[phtpt['g'],
                                                          phtpt['i']])
        MassEstimate_obj.set_target(get_target(zcmb=self.z))
        # Gives [mass, mass_err_sown, mass_err_up]
        self.hg_mass = MassEstimate_obj.get_estimate()

        return self.hg_mass

    # ------------------------------------------------------------------- #
    #                               PLOTTER                               #
    # ------------------------------------------------------------------- #

    def show(self, ax=None, band="r",
             ellipse=True, scaleup=2.5,
             ell_color=(1, 0, 0, 0.3)):
        """ """
        fig = plt.figure(figsize=[5, 5])
        if ax is not None:
            ax = ax
        else:
            ax = fig.add_axes([0.12, 0.12, 0.8, 0.8])

        self.cutout[band].show(show_sepobjects=True, logscale=False, ax=ax)

        if ellipse:
            hg_patch = Ellipse([self.hg_ellipse[0], self.hg_ellipse[1]],
                               2*self.hg_ellipse[2]*scaleup,
                               2*self.hg_ellipse[3]*scaleup,
                               self.hg_ellipse[4]*units.rad.to('deg'),
                               color=ell_color)

            ax.add_patch(hg_patch)

        return fig

    # =================================================================== #
    #                              Properties                             #
    # =================================================================== #

    @property
    def coordinates(self):
        """ """
        return [self.ra, self.dec]

    # - Images
    @property
    def cutout(self):
        """ """
        if not hasattr(self, "_cutout"):
            raise AttributeError(
                "No cutout loaded yet. " +
                "Run self.download() or set self._cutout")
        return self._cutout
