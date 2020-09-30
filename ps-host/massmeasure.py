#! /usr/bin/env python
# -*- coding: utf-8 -*-

# import os
# import pandas as pd
import sep
import numpy as np
import matplotlib.pyplot as mpl

from matplotlib.patches import Ellipse
from astrobject.utils.tools import flux_to_mag

# PanSTARRS zero point according to
# https://coolwiki.ipac.caltech.edu/index.php/
# Central_wavelengths_and_zero_points
ps1_zp = {'g': 0.4810,
          'r': 0.6170,
          'i': 0.7520,
          'z': 0.8660,
          'y': 0.9620}

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
    mm = MassMeasure(ra, dec)
    mm.download()  # if necessary; else run mm._cutout = cutout

    # Optional
    mm.show()

    mm.get_ellipses()  # also defines the associated host galaxy
                       # alternatively get HG with mm.get_hostgalaxy()

    # Optional
    mm.show()  # If ellipses have been extracted they are shown

    mm.get_hostgalaxy()
    mm.get_hostmass()

    # Optional:
    mm.test_pixels(n=100)'''

    # =================================================================== #
    #                               Initial                               #
    # =================================================================== #

    def __init__(self, ra, dec, cutout=None):
        '''
        Sets the object's ra, dec coordinates and instantiates a PS1Target
        object from `pymage.panstarrs`, downloading the needed cutout if
        not given as an entry parameter.
        '''
        if not _HAS_PYMAGE:
            raise ImportError('This method needs pymage: `pip install pymage`')

        self.target = panstarrs.PS1Target.from_coord(ra, dec)
        self.ra = self.target.coordinate.ra.deg
        self.dec = self.target.coordinate.dec.deg
        if cutout is None:
            self.target.download_cutout(load_weight=True)
        else:
            self.target._cutout = cutout
        self._cutout = self.target.imgcutout
        self.x, self.y = self.target.imgcutout['r'].coords_to_pixel(
            self.ra, self.dec)
        sep_obj = self.target.sep_extract(returnobjects=True)
        self.ellipses = sep_obj.get_ellipses_values().T

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
        if not hasattr(self, _cutout):
            raise AttributeError(
                "No cutout loaded yet. " +
                "Run self.download() or set self._cutout")
        else:
            return self.target.imgcutout['r'].count_to_flux(count)

    # ------------------------------------------------------------------- #
    #                               GETTER                                #
    # ------------------------------------------------------------------- #

    def get_hostgalaxy_params(self):
        '''
        Gives x, y, a, b, theta of host galaxy
        determined to be the object with smallest DLR from target
        '''
        DLR_list = [self.get_DLR(self.x, self.y,
                                 self.ellipses[i][0], self.ellipses[i][1],
                                 self.ellipses[i][2], self.ellipses[i][3],
                                 self.ellipses[i][4])
                    for i in range(len(self.ellipses))]
        hg_ind = np.argmin(DLR_list)
        self.hg_ellipse = self.ellipses[hg_ind]

        bands = list(self.cutout.keys())
        counts_res = [sep.sum_ellipse(self.cutout[band],
                                      self.hg_ellipse[0], self.hg_ellipse[1],
                                      self.hg_ellipse[2], self.hg_ellipse[3],
                                      self.hg_ellipse[4],
                                      var=self.cutout[band])
                      for band in bands]
        self.hg_count = {band: [counts_res[band][0], counts_res[band][1]]
                         for band in bands}

        self.hg_flux = {band: [self.count_to_flux(self.hg_count[band][0]),
                               self.count_to_flux(self.hg_count[band][1])]
                        for band in bands}

        self.hg_mag = {band: [flux_to_mag(self.hg_flux[band][0],
                                          self.hg_flux[band][1],
                                          zp=ps1_zp[band])[0],
                              flux_to_mag(self.hg_flux[band][0],
                                          self.hg_flux[band][1],
                                          zp=ps1_zp[band])[1]]
                       for band in bands}

    # ------------------------------------------------------------------- #
    #                               PLOTTER                               #
    # ------------------------------------------------------------------- #

    def show(self, ax=None, band="r", show_coord=None,
             source='sep', ellipse=True, ell_color="k",
             coord_color="C1", scaleup=3, **kwargs):
        """ """
        if self.has_cutout():
            img = self.cutout[band]
            ax = img.show(ax=ax, show_sepobjects=False)['ax']
            inpixel = True
            has_img = True
        elif ax is None:
            fig = mpl.figure(figsize=[5, 5])
            ax = fig.add_subplot(111)
            inpixel = False
            has_img = False

        # - RA,Dec or x,y
        if inpixel:
            x, y = img.coords_to_pixel(
                self.coordinate.ra.deg, self.coordinate.dec.deg)
        else:
            x, y = self.coordinate.ra.deg, self.coordinate.dec.deg

        ax.scatter(x, y, **{**dict(marker="+", color='k', s=100), **kwargs})

        if ellipse:
            x, y, a, b, t = self.get_nearest_ellipse(
                source=source, inpixel=inpixel)
            ax.add_patch(Ellipse([x, y],
                                 2*a*scaleup, 2*b*scaleup,
                                 t*units.rad.to("deg"),
                                 facecolor="None", edgecolor=ell_color, lw=2))
            if not has_img:
                ax.scatter(x, y, marker=".", color=ell_color)
                ax.set_xlim(x-3*a, x+3*a)
                ax.set_ylim(y-3*a, y+3*a)

        if show_coord is not None:
            if inpixel:
                x, y = img.coords_to_pixel(*show_coord)
            else:
                x, y, = show_coord

            ax.scatter(x, y, marker="x", color=pt_color, s=100)

        return ax.figure

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

    @property
    def ellipse(self):
        """ """
        if not hasattr(self, "_ellipse") or self._ellipse is None:
            self.derive_ellipse()
        return self._ellipse

    @property
    def catdata(self):
        """ """
        if not hasattr(self, "_catdata"):
            self._catdata = None
        return self._catdata

    @property
    def _extended_cat_set(self):
        """ """
        if not hasattr(self, "_is_extended_cat_set"):
            self._is_extended_cat_set = False
        return self._is_extended_cat_set

    @property
    def sep(self):
        """ """
        if not hasattr(self, "_sep"):
            self._sep = None
        return self._sep

    @property
    def galcat(self):
        """ """
        if not hasattr(self, "_galcat"):
            self._galcat = None
        return self._galcat
