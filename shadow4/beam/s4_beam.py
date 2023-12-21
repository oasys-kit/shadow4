"""
Defines the shadow4 beam and associated methods and utilities.

VERY IMPORTANT: Column 11 (index 10) is wavenumber (cm^-1) as internally in Shadow.

                Photon energy in eV is now column 26 (index 25).

"""
import numpy
import scipy.constants as codata
from numpy.testing import assert_almost_equal

import h5py
import time
import os

from syned.beamline.shape import Rectangle, Ellipse, TwoEllipses, Circle

A2EV = 2.0 * numpy.pi / (codata.h * codata.c / codata.e * 1e2)

class S4Beam(object):
    """
    Implements a beam. Internally it is an array of N rays and 18 columns.

    Parameters
    ----------
    N : int, optional
        The number of rays.

    array numpy array, optional
        The numpy array (N,18) with the data.

    See Also
    --------
    shadow4.S4Beam.column_names :
        columns contents.

    """
    def __init__(self, N=1000, array=None):
        if array is not None:
            N, ncol = array.shape
            if ncol != 18:
                raise Exception ("Bad array shape: must be (npoints,18)")
            self.rays = array.copy()
        else:
            self.rays = numpy.zeros((N,18))

    @classmethod
    def initialize_from_array(cls, array):
        """
        Creates an S4Beam instance from an array.

        Parameters
        ----------
        array : numpy array
            array to initialize the S4beam.

        Returns
        -------
            an instance of S4Beam.

        """
        if array.shape[1] != 18:
            raise Exception("Bad array shape: must be (npoints,18)")
        return S4Beam(array=array)

    @classmethod
    def initialize_as_pencil(cls, N=1000):
        """

        Parameters
        ----------
        N : int, optional
            Number of rays.

        Returns
        -------
        S4Beam instance
            A pencil beam (zero cross section, zero divergence).

        """
        beam = S4Beam(N)
        beam.set_column(5,1.0) # Vy
        beam.set_column(7,1.0) # Es
        beam.set_column(10,1.0) # flag
        beam.set_column(11,2*numpy.pi/1e-8) # wavenumber (1 A)
        beam.set_column(12,numpy.arange(N,dtype=float)) # index

        return beam


    def duplicate(self):
        """
        Duplicates a S4Beam instance.

        Returns
        -------
        S4Beam instance
            A copy of the S4Beam instance.

        """
        return S4Beam.initialize_from_array(self.rays.copy())

    #
    # getters
    #
    def get_rays(self, nolost=0):
        """
        Returns a numpy array (N,18) with the rays.

        Parameters
        ----------
        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        Returns
        -------
        numpy array (npoints,18)
            The rays.

        """

        if nolost == 0:
            return self.rays.copy()
        elif nolost == 1:
            f  = numpy.where(self.rays[:,9] > 0.0)
            if len(f[0])==0:
                print ('S4Beam.get_rays: no GOOD rays, returning empty array')
                return numpy.empty(0)
            else:
                return self.rays[f[0],:].copy()
        elif nolost == 2:
            f  = numpy.where(self.rays[:,9] < 0.0)
            if len(f[0])==0:
                print ('S4Beam.get_rays: no BAD rays, returning empty array')
                return numpy.empty(0)
            else:
                return self.rays[f[0],:].copy()

    @property
    def rays_good(self):
        """
        Returns a numpy array with the good rays.

        Returns
        -------
        numpy array (ngoodrays,18)
            The good rays.

        """
        return self.get_rays(nolost=1)

    @property
    def rays_bad(self):
        """
        Returns a numpy array with the bad rays.

        Returns
        -------
        numpy array (nbadrays,18)
            The bad rays.

        """
        return self.get_rays(nolost=2)

    def get_number_of_rays(self, nolost=0):
        """
        Returns the number of rays.

        Parameters
        ----------
        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        Returns
        -------
        int
            The number of rays.

        """
        try:
            w = self.get_column(10)
        except Exception:
            print("Error: Empty beam...")
            return 0

        if nolost == 0:
            return w.size
        if nolost == 1:
            return numpy.array(numpy.where(w >= 0)).size
        if nolost == 2:
            return numpy.array(numpy.where(w < 0)).size

        return self.rays.shape[0]

    @property
    def N(self):
        """
        Returns the total number of rays.

        Returns
        -------
        int
            The total number of rays.

        """
        return self.get_number_of_rays(nolost=0)

    @property
    def Ngood(self):
        """
        Returns the number of good rays.

        Returns
        -------
        int
            The number of good rays.

        """
        return self.get_number_of_rays(nolost=1)

    @property
    def Nbad(self):
        """
        Returns the number of bad rays.

        Returns
        -------
        int
            The number of bad rays.

        """
        return self.get_number_of_rays(nolost=2)

    def get_photon_energy_eV(self, nolost=0):
        """
        Returns a numpy array with the photon energy of the rays in eV.

        Parameters
        ----------
        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        Returns
        -------
        numpy array (npoints)
            The photon energies (in eV) of the rays.

        """
        return self.get_column(11, nolost=nolost) / A2EV

    def get_photon_wavelength(self, nolost=0):
        """
        Returns a numpy array with the photon energy of the rays in eV.

        Parameters
        ----------
        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        Returns
        -------
        numpy array (npoints)
            The wavelengths (in m) of the rays.

        """
        return 2 * numpy.pi / self.get_column(11, nolost=nolost) * 1e-2

    def get_intensity(self, nolost=0, polarization=0):
        """
        Returns a numpy array with the intensities of the rays.

        Parameters
        ----------
        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        polarization : int, optional
            0=total, 1=sigma, 2=pi.

        Returns
        -------
        numpy array (npoints)
            The intensities of the rays.

        """
        if polarization == 0:
            w = self.get_column(23,nolost=nolost)
        elif polarization == 1:
            w = self.get_column(24, nolost=nolost)
        elif polarization == 2:
            w = self.get_column(25, nolost=nolost)
        return w.sum()

    def get_column(self, column, nolost=0):
        """
        Returns a numpy array with the values of the rays in a given column.

        The column number can be 1:18 (stored data) or > 18 for other parameters
        calculated from the 18 stored ones.

        Parameters
        ----------
        column : int
            Number of column (starting with 1). Possible choice for column are:

             1   X spatial coordinate [user's unit]
             2   Y spatial coordinate [user's unit]
             3   Z spatial coordinate [user's unit]
             4   Xp direction or divergence [rads]
             5   Yp direction or divergence [rads]
             6   Zp direction or divergence [rads]
             7   X component of the electromagnetic vector (s-polariz)
             8   Y component of the electromagnetic vector (s-polariz)
             9   Z component of the electromagnetic vector (s-polariz)
            10   Lost ray flag
            11   wavenumber (2 pi / lambda[cm]) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            12   Ray index
            13   Optical path length
            14   Phase (s-polarization) in rad
            15   Phase (p-polarization) in rad
            16   X component of the electromagnetic vector (p-polariz)
            17   Y component of the electromagnetic vector (p-polariz)
            18   Z component of the electromagnetic vector (p-polariz)

            19   Wavelength [A]
            20   R= SQRT(X^2+Y^2+Z^2)
            21   angle from Y axis
            22   the magnitude of the Electromagnetic vector
            23   |E|^2 (total intensity)
            24   total intensity for s-polarization
            25   total intensity for p-polarization
            26   photon energy in eV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            27   K = 2 pi / lambda * col4 [A^-1]
            28   K = 2 pi / lambda * col5 [A^-1]
            29   K = 2 pi / lambda * col6 [A^-1]
            30   S0-stokes = |Ep|^2 + |Es|^2
            31   S1-stokes = |Ep|^2 - |Es|^2
            32   S2-stokes = 2 |Es| |Ep| cos(phase_s-phase_p)
            33   S3-stokes = 2 |Es| |Ep| sin(phase_s-phase_p)
            34   Power = intensity(col 23) * energy (col 11)
            35   Angle-X with Y: |arcsin(X')|
            36   Angle-Z with Y: |arcsin(Z')|
            37   Angle-X with Y: |arcsin(X') - mean(arcsin(X'))|
            38   Angle-Z with Y: |arcsin(Z') - mean(arcsin(Z'))|
            39   Phase difference in rad: Phase (s-polarization) - Phase (p-polarization)
            40   Complex amplitude of the electric vector (s-polarization)
            41   Complex amplitude of the electric vector (p-polarization)

            -11: column 26

        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        Returns
        -------
        numpy array
            the required array (a copy).

        """
        if column == -11: column = 26

        if column <= 18:
            out = self.rays[:,column-1]
        else:
            column_index = column - 1
            ray = self.rays

             # if colu mn_index==10: out =  ray[:,column_index]/A2EV
            if column == 19: out =  2*numpy.pi*1.0e8/ray[:,10]
            if column == 20: out =  numpy.sqrt(ray[:,0]*ray[:,0]+ray[:,1]*ray[:,1]+ray[:,2]*ray[:,2])
            if column == 21: out =  numpy.arccos(ray[:,4])
            if column == 22: out =  numpy.sqrt(numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8,15,16,17] ]),axis=0))
            if column == 23: out =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8,15,16,17] ]),axis=0)
            if column == 24: out =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8] ]),axis=0)
            if column == 25: out =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [15,16,17] ]),axis=0)
            if column == 26: out =  ray[:,10]/A2EV
            if column == 27: out =  ray[:,3]*ray[:,10]*1.0e8
            if column == 28: out =  ray[:,4]*ray[:,10]*1.0e8
            if column == 29: out =  ray[:,5]*ray[:,10]*1.0e8
            if column == 30:
                E2s = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8] ]),axis=0)
                E2p = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [15,16,17] ]),axis=0)
                out =  E2p+E2s
            if column ==31:
                E2s = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8] ]),axis=0)
                E2p = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [15,16,17] ]),axis=0)
                out =  E2p-E2s
            if column == 32:
                E2s = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8] ]),axis=0)
                E2p = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [15,16,17] ]),axis=0)
                Cos = numpy.cos(ray[:,13]-ray[:,14])
                out =  2*numpy.sqrt(E2s*E2p)*Cos
            if column == 33:
                E2s = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8] ]),axis=0)
                E2p = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [15,16,17] ]),axis=0)
                Sin = numpy.sin(ray[:,13]-ray[:,14])
                out =  2*numpy.sqrt(E2s*E2p)*Sin

            if column == 34:
                out =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8,15,16,17] ]),axis=0) *\
                ray[:,10]/A2EV

            if column == 35:
                out = numpy.abs(numpy.arcsin(ray[:,3]))
            if column == 36:
                out = numpy.abs(numpy.arcsin(ray[:,5]))
            if column == 37:
                f =  ray[:, 10 - 1]
                w =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8,15,16,17] ]),axis=0)
                xp = ray[:, 4 - 1]
                if nolost == 1:
                    findices  = numpy.where(f > 0.0)
                    if len(findices[0])==0:
                        col_mean = numpy.average(xp, weights=w)
                    else:
                        col_mean = numpy.average(xp[findices], weights=w[findices])
                else:
                    col_mean = numpy.average(xp, weights=w)
                out = numpy.abs(numpy.arcsin(xp - col_mean))
            if column == 38:
                f  = ray[:, 10 - 1]
                w  = numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8,15,16,17] ]),axis=0)
                zp = ray[:, 6 - 1]
                if nolost == 1:
                    findices  = numpy.where(f > 0.0)
                    if len(findices[0])==0:
                        col_mean = numpy.average(zp, weights=w)
                    else:
                        col_mean = numpy.average(zp[findices], weights=w[findices])
                else:
                    col_mean = numpy.average(zp, weights=w)

                out = numpy.abs(numpy.arcsin(zp - col_mean))
            if column == 39:
                out = ray[:, 14 - 1] - ray[:, 15 - 1]
            if column == 40:
                out =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [6,7,8] ]),axis=0)
                out *= numpy.exp(1j * ray[:, 14 - 1])
            if column == 41:
                out =  numpy.sum(numpy.array([ ray[:,i]*ray[:,i] for i in [15,16,17] ]),axis=0)
                out *= numpy.exp(1j * ray[:, 15 - 1])

        if nolost == 0:
            return out.copy()

        if nolost == 1:
            f  = numpy.where(self.rays[:,9] > 0.0)
            if len(f[0])==0:
                print ('Beam.get_column: no GOOD rays, returning empty array')
                return numpy.empty(0)
            return out[f].copy()

        if nolost == 2:
            f  = numpy.where(self.rays[:,9] < 0.0)
            if len(f[0])==0:
                print ('Beam.get_column: no BAD rays, returning empty array')
                return numpy.empty(0)
            return out[f].copy()

    def get_columns(self, columns, nolost=0):
        """
        Returns a numpy array with the values of the rays several selected column.

        The column numbers can be 1:18 (stored data) or > 18 for other parameters
        calculated from the 18 stored ones.

        Parameters
        ----------
        columns : list
            The number of the columns (column numbers start from 1).

        nolost : int, optional
            0=return all rays, 1=Return only good rays (non-lost rays), 2=Return only lost rays.

        Returns
        -------
        numpy array
            the required array (N, len(columns)).

        See Also
        --------
        shadow4.S4Beam.get_column

        """
        ret = []
        if isinstance(columns, int): return self.get_column(columns,nolost=nolost)
        for c in columns:
            ret.append(self.get_column(c,nolost=nolost))
        return numpy.array(tuple(ret))

    def get_standard_deviation(self,col, nolost=1, ref=0):
        """
        Returns the standard deviation of one variable in the beam.

        Parameters
        ----------
        col : int
            The column number.


        nolost : int, optional
            0=use all rays, 1=use only good rays (non-lost rays), 2=use only lost rays.

        ref : int, optional
            ref: 0 = no weight, 1=weight with intensity (col23)

        Returns
        -------
        float
            the st dev.

        """
        x = self.get_column(col,nolost=nolost)
        if ref == 0:
            return x.std()
        else:
            w = self.get_column(23,nolost=nolost)
            average = numpy.average(x, weights=w)
            variance = numpy.average( (x-average)**2, weights=w)
            return(numpy.sqrt(variance))

    def intensity(self, nolost=0):
        """
        Returns the intensity of the beam.

        Parameters
        ----------

        nolost : int, optional
            0=use all rays, 1=use only good rays (non-lost rays), 2=use only lost rays.


        Returns
        -------
        float
            total intensity.

        """
        w = self.get_column(23,nolost=nolost)
        return w.sum()

    def get_good_range(self, icol, nolost=0):
        """
        Computes a good limits for a given column. Typically used for plots.

        Parameters
        ----------
        icol: int
            the column number (SHADOW convention, starting from 1)

        nolost : int, optional
            0=use all rays, 1=use only good rays (non-lost rays), 2=use only lost rays.

        Returns
        -------
        list
            [rmin,rmax] the selected range

        """
        col = self.get_column(icol, nolost=nolost)
        if col.size == 0:
            return [-1, 1]
        rmin = min(col)
        rmax = max(col)
        if rmin > 0.0:
            rmin = rmin * 0.95
        else:
            rmin = rmin * 1.05
        if rmax < 0.0:
            rmax = rmax * 0.95
        else:
            rmax = rmax * 1.05
        if rmin == rmax:
            rmin = rmin * 0.95
            rmax = rmax * 1.05
        if rmin == 0.0 and rmax == 0:
            rmin = -1.0
            rmax = 1.0
        return [rmin, rmax]

    def histo1(self, col, xrange=None, nbins=50, nolost=0, ref=0,
               write=None, factor=1.0, calculate_widths=1, calculate_hew=0):
        """
        Calculate the histogram of a column, simply counting the rays, or weighting with another column.

        It returns a dictionary which contains the histogram data.

        Parameters
        ----------
        col : int
            the number of the chosen column.

        xrange : 2 elements tuple or list, optional
            tuple or list describing the interval of interest for x, the data read from the chosen column.
                      (default: None, thus using min and max of the array)

        nbins : int, optional
            number of bins of the histogram.

        nolost : int, optional
            0=use all rays, 1=use only good rays (non-lost rays), 2=use only lost rays.

        ref : int (or str), optional
                 0, None, "no", "NO" or "No":   only count the rays.
                 23, "Yes", "YES" or "yes":     weight with intensity (look at col=23 |E|^2 total intensity).
                 other value: use that column as weight.

        write : str, optional
                file name with the histogram (default=None, do not write any file).

        factor : float, optional
            a scalar factor to multiply the selected column before histogramming
            (e.g., for changing scale from cm to um then factor=1e4).

        calculate_widths : int, optional
            0: do not calculate full-width at half-maximum (FWHM). 1=Do calculate FWHM.

        calculate_hew : int, optional
            0: do not calculate Half-Energy Width (HEW). 1=Do calculate HEW.

        Returns
        -------
        dict
            a python dictionary with the calculated histogram. Some of the keys are:
                 'error', 'col', 'write', 'nolost', 'nbins', 'xrange', 'factor'
                 'histogram', 'bins', 'histogram_sigma', 'bin_center', 'bin_left', 'bin_right',
                 'intensity', 'fwhm', 'nrays', 'good_rays'.

        """
        # initialize return value
        ticket = {'error': 1}

        # coli = col - 1
        if ref == None: ref = 0
        if ref == "No": ref = 0
        if ref == "NO": ref = 0
        if ref == "no": ref = 0

        if ref == "Yes": ref = 23
        if ref == "YES": ref = 23
        if ref == "yes": ref = 23

        if ref == 1:
            print(
                "Shadow.S4Beam.histo1: Warning: weighting with column 1 (X) [not with intensity as may happen in old versions]")

        # copy the inputs
        ticket['col'] = col
        ticket['write'] = write
        ticket['nolost'] = nolost
        ticket['nbins'] = nbins
        ticket['xrange'] = xrange

        ticket['factor'] = factor
        ticket['ref'] = ref

        if ref == 0:
            x = self.get_column(col, nolost=nolost)
            w = numpy.ones(len(x))
        else:
            x, w = self.get_columns((col, ref), nolost=nolost)

        if factor != 1.0: x *= factor

        if len(x) == 0: # no rays
            ticket['error'] = 0
            ticket['histogram'] = ticket['bins'] = ticket['bin_center'] = \
                ticket['histogram_path'] =  ticket['bin_path'] = numpy.empty(0)
            ticket['histogram_sigma'] = 0.0
            ticket['bin_left'] = ticket['bin_right'] = ticket['intensity'] = numpy.nan
            ticket['xrange'] = xrange
            ticket['nrays'] = self.get_number_of_rays(nolost=0)
            ticket['good_rays'] = 0
            ticket['fwhm'] = None

            if calculate_widths > 0:
                ticket['fwhm_subpixel'] = None
                ticket['fwhm_coordinates'] = ticket['fwhm_subpixel_coordinates'] =(numpy.nan, numpy.nan)
            if calculate_widths == 2: ticket["fw25%m"] = ticket["fw75%m"] = None
            if calculate_hew:         ticket["hew"] = numpy.nan
        else:
            if xrange == None: xrange = [x.min(), x.max()]

            h, bins = numpy.histogram(x, bins=nbins, range=xrange, weights=w)
            # evaluate the histogram with squares of the weight for error calculations
            h2, _ = numpy.histogram(x, bins=nbins, range=xrange, weights=(w * w))

            # Evaluation of histogram error.
            # See Pag 17 in Salvat, Fernandez-Varea and Sempau
            # Penelope, A Code System for Monte Carlo Simulation of Electron and Photon Transport, AEN NEA  (2003)
            #
            # See James, Rep. Prog. Phys., Vol 43 (1980) pp 1145-1189 (special attention to pag. 1184)
            h_sigma = numpy.sqrt(h2 - h * h / float(len(w)))

            bin_center = bins[:-1] + (bins[1] - bins[0]) * 0.5

            if write != None and write != "":
                f = open(write, 'w')
                f.write('#F %s \n' % (write))
                f.write('#C This file has been created using Shadow.Beam.histo1() \n')
                f.write('#C COLUMN 1 CORRESPONDS TO ABSCISSAS IN THE CENTER OF EACH BIN\n')
                f.write('#C COLUMN 2 CORRESPONDS TO ABSCISSAS IN THE THE LEFT CORNER OF THE BIN\n')
                f.write('#C COLUMN 3 CORRESPONDS TO INTENSITY\n')
                f.write('#C COLUMN 4 CORRESPONDS TO ERROR: SIGMA_INTENSITY\n')
                f.write('#C col = %d\n' % (col))
                f.write('#C nolost = %d\n' % (nolost))
                f.write('#C nbins = %d\n' % (nbins))
                f.write('#C ref = %d\n' % (ref), )
                f.write(' \n')
                f.write('#S 1 histogram\n')
                f.write('#N 4\n')
                f.write('#L X1  X2  Y  YERR\n')
                for i in range(len(h)):
                    f.write('%f\t%f\t%f\t%f\n' % ((bins[i] + bins[i + 1]) * 0.5, bins[i], h[i], h_sigma[i]))
                f.close()
                print('histo1: file written to disk: %s' % (write))

            # output
            ticket['error'] = 0
            ticket['histogram'] = h
            ticket['bins'] = bins
            ticket['histogram_sigma'] = h_sigma
            ticket['bin_center'] = bin_center
            ticket['bin_left'] = bins[:-1]
            ticket['bin_right'] = bins[:-1] + (bins[1] - bins[0])
            ticket['xrange'] = xrange
            ticket['intensity'] = self.intensity(nolost=nolost)
            ticket['fwhm'] = None
            ticket['nrays'] = self.get_number_of_rays(nolost=0)
            ticket['good_rays'] = self.get_number_of_rays(nolost=1)

            # for practical purposes, writes the points the will define the histogram area
            tmp_b = []
            tmp_h = []
            for s, t, v in zip(ticket["bin_left"], ticket["bin_right"], ticket["histogram"]):
                tmp_b.append(s)
                tmp_h.append(v)
                tmp_b.append(t)
                tmp_h.append(v)
            ticket['histogram_path'] = numpy.array(tmp_h)
            ticket['bin_path'] = numpy.array(tmp_b)

            if calculate_widths > 0:
                # CALCULATE fwhm
                tt = numpy.where(h >= max(h) * 0.5)
                if h[tt].size > 1:
                    binSize = bins[1] - bins[0]
                    ticket['fwhm'] = binSize * (tt[0][-1] - tt[0][0])
                    ticket['fwhm_coordinates'] = (bin_center[tt[0][0]], bin_center[tt[0][-1]])

                # CALCULATE fwhm with subpixel resolution (as suggested by A Wojdyla)
                ixl_e = tt[0][0]
                ixr_e = tt[0][-1]
                try:
                    xl = ixl_e - (h[ixl_e] - max(h) * 0.5) / (h[ixl_e] - h[ixl_e - 1])
                    xr = ixr_e - (h[ixr_e] - max(h) * 0.5) / (h[ixr_e + 1] - h[ixr_e])
                    ticket['fwhm_subpixel'] = binSize * numpy.abs(xr - xl)
                    ticket['fwhm_subpixel_coordinates'] = \
                        (numpy.interp(xl, range(bin_center.size), bin_center),
                         numpy.interp(xr, range(bin_center.size), bin_center))
                except:
                    ticket['fwhm_subpixel'] = None

            if calculate_widths == 2:
                # CALCULATE FW at 25% HEIGHT
                tt = numpy.where(h >= max(h) * 0.25)
                if h[tt].size > 1:
                    binSize = bins[1] - bins[0]
                    ticket['fw25%m'] = binSize * (tt[0][-1] - tt[0][0])
                else:
                    ticket["fw25%m"] = None

                # CALCULATE FW at 75% HEIGHT
                tt = numpy.where(h >= max(h) * 0.75)
                if h[tt].size > 1:
                    binSize = bins[1] - bins[0]
                    ticket['fw75%m'] = binSize * (tt[0][-1] - tt[0][0])
                else:
                    ticket["fw75%m"] = None

            if calculate_hew:
                # CALCULATE HALF-ENERGY-WIDTH
                cdf = numpy.cumsum(ticket["histogram"])
                cdf /= cdf.max()
                # hew is two times  the x value that has cdf=0.5 (Eq. 9 in https://arxiv.org/pdf/1505.07474v2.pdf)
                hew = 2 * float(bin_center[numpy.argwhere(cdf > 0.5)][0])
                ticket["hew"] = hew

        return ticket

    def histo2(self, col_h, col_v, nbins=25, ref=23, nbins_h=None, nbins_v=None, nolost=0,
               xrange=None, yrange=None, calculate_widths=1):
        """
        Performs 2d histogram. Typically used to prepare data for a plotxy plot.
        It uses histogram2d for calculations.

        Note that this shadow4.S4Beam.histo2 was previously called Shadow.Beam.plotxy in shadow3.

        Parameters
        ----------
        col_h: int
            the horizontal column.

        col_v: int
            the vertical column.

        nbins: int
            The number of bins.

        ref : int (or str), optional
                 0, None, "no", "NO" or "No":   only count the rays.
                 23, "Yes", "YES" or "yes":     weight with intensity (look at col=23 |E|^2 total intensity).
                 other value: use that column as weight.

        nbins_h: int
            number of bins in H.

        nbins_v: int
            number of bins in V.

        nolost : int, optional
            0=use all rays, 1=use only good rays (non-lost rays), 2=use only lost rays.

        xrange: tuple or list:
            range for H.

        yrange: tuple or list
            range for V.

        calculate_widths: int
            0=No, 1=calculate FWHM (default), 2=Calculate FWHM and FW at 25% and 75% if Maximum.

        Returns
        -------
        dict
            a dictionary with the histogram and all the data needed (e.g. for a plot).

        """
        ticket = {'error':1}

        if ref == None: ref = 0
        if ref == "No": ref = 0
        if ref == "NO": ref = 0
        if ref == "no": ref = 0

        if ref == "Yes": ref = 23
        if ref == "YES": ref = 23
        if ref == "yes": ref = 23

        if ref == 1:
              print("shadow4.S4Beam.histo2: Warning: weighting with column 1 (X) [not with intensity]")

        if nbins_h == None: nbins_h = nbins
        if nbins_v == None: nbins_v = nbins

        # copy the inputs
        ticket['col_h'] = col_h
        ticket['col_v'] = col_v
        ticket['nolost'] = nolost
        ticket['nbins_h'] = nbins_h
        ticket['nbins_v'] = nbins_v
        ticket['ref'] = ref

        (col1,col2) = self.get_columns((col_h,col_v),nolost=nolost)

        if len(col1) == 0 or len(col2) == 0:# no rays
            ticket['xrange'] = xrange
            ticket['yrange'] = yrange
            ticket['bin_h_edges'] = ticket['bin_v_edges'] = ticket['bin_h_left'] = \
                ticket['bin_v_left'] = ticket['bin_h_right'] = ticket['bin_v_right'] = \
                ticket['histogram'] = ticket['histogram_h'] = ticket['histogram_v'] = numpy.empty(0)
            ticket['bin_h_center'] = ticket['bin_v_center'] = ticket['intensity'] = numpy.nan
            ticket['nrays'] = self.get_number_of_rays(nolost=0)
            ticket['good_rays'] = 0
            ticket['fwhm_h'] = ticket['fwhm_v'] = None
            if calculate_widths > 0:  ticket['fwhm_coordinates_h'] = ticket['fwhm_coordinates_v'] = (numpy.nan, numpy.nan)
            if calculate_widths == 2: ticket["fw25%m_h"] = ticket["fw75%m_h"] = ticket["fw25%m_v"] = ticket["fw75%m_v"] = None
        else:
            if xrange==None: xrange = self.get_good_range(col_h,nolost=nolost)
            if yrange==None: yrange = self.get_good_range(col_v,nolost=nolost)
            if ref == 0:
                weights = col1*0+1
            else:
                weights = self.get_column(ref,nolost=nolost)

            (hh,xx,yy) = numpy.histogram2d(col1, col2, bins=[nbins_h,nbins_v], range=[xrange,yrange], normed=False, weights=weights)

            ticket['xrange'] = xrange
            ticket['yrange'] = yrange
            ticket['bin_h_edges'] = xx
            ticket['bin_v_edges'] = yy
            ticket['bin_h_left'] = numpy.delete(xx,-1)
            ticket['bin_v_left'] = numpy.delete(yy,-1)
            ticket['bin_h_right'] = numpy.delete(xx,0)
            ticket['bin_v_right'] = numpy.delete(yy,0)
            ticket['bin_h_center'] = 0.5*(ticket['bin_h_left']+ticket['bin_h_right'])
            ticket['bin_v_center'] = 0.5*(ticket['bin_v_left']+ticket['bin_v_right'])
            ticket['histogram'] = hh
            ticket['histogram_h'] = hh.sum(axis=1)
            ticket['histogram_v'] = hh.sum(axis=0)
            ticket['intensity'] = self.intensity(nolost=nolost)
            ticket['nrays'] = self.get_number_of_rays(nolost=0)
            ticket['good_rays'] = self.get_number_of_rays(nolost=1)

            # CALCULATE fwhm

            if calculate_widths > 0:
                h = ticket['histogram_h']
                tt = numpy.where(h>=max(h)*0.5)
                if h[tt].size > 1:
                    binSize = ticket['bin_h_center'][1]-ticket['bin_h_center'][0]
                    ticket['fwhm_h'] = binSize*(tt[0][-1]-tt[0][0])
                    ticket['fwhm_coordinates_h'] = (ticket['bin_h_center'][tt[0][0]],ticket['bin_h_center'][tt[0][-1]])
                else:
                    ticket["fwhm_h"] = None

                h = ticket['histogram_v']
                tt = numpy.where(h>=max(h)*0.5)
                if h[tt].size > 1:
                    binSize = ticket['bin_v_center'][1]-ticket['bin_v_center'][0]
                    ticket['fwhm_v'] = binSize*(tt[0][-1]-tt[0][0])
                    ticket['fwhm_coordinates_v'] = (ticket['bin_v_center'][tt[0][0]],ticket['bin_v_center'][tt[0][-1]])
                else:
                    ticket["fwhm_v"] = None

            if calculate_widths == 2:
                # CALCULATE FW at 25% HEIGHT
                h = ticket['histogram_h']
                tt = numpy.where(h>=max(h)*0.25)
                if h[tt].size > 1:
                    binSize = ticket['bin_h_center'][1]-ticket['bin_h_center'][0]
                    ticket['fw25%m_h'] = binSize*(tt[0][-1]-tt[0][0])
                else:
                    ticket["fw25%m_h"] = None

                h = ticket['histogram_v']
                tt = numpy.where(h>=max(h)*0.25)
                if h[tt].size > 1:
                    binSize = ticket['bin_v_center'][1]-ticket['bin_v_center'][0]
                    ticket['fw25%m_v'] = binSize*(tt[0][-1]-tt[0][0])
                else:
                    ticket["fw25%m_v"] = None

                # CALCULATE FW at 75% HEIGHT
                h = ticket['histogram_h']
                tt = numpy.where(h>=max(h)*0.75)
                if h[tt].size > 1:
                    binSize = ticket['bin_h_center'][1]-ticket['bin_h_center'][0]
                    ticket['fw75%m_h'] = binSize*(tt[0][-1]-tt[0][0])
                else:
                    ticket["fw75%m_h"] = None

                h = ticket['histogram_v']
                tt = numpy.where(h>=max(h)*0.75)
                if h[tt].size > 1:
                    binSize = ticket['bin_v_center'][1]-ticket['bin_v_center'][0]
                    ticket['fw75%m_v'] = binSize*(tt[0][-1]-tt[0][0])
                else:
                    ticket["fw75%m_v"] = None

        return ticket

    #
    # setters
    #
    def set_column(self, column, value):
        """
        Sets a given array in a given column number.

        Parameters
        ----------
        column : int
            The column number.

        value : numpy array (or a float or int scalar)
            The values to be set.

        """
        self.rays[:,column-1] = value

    def set_photon_energy_eV(self, energy_eV):
        """
        Sets photon energies from a given array.

        Parameters
        ----------

        value : numpy array (or a float or int scalar)
            The values of the photon energies in eV.

        """
        self.rays[:,10] = energy_eV * A2EV

    def set_photon_wavelength(self,wavelength):
        """
        Sets photon wavelengths from a given array.

        Parameters
        ----------

        value : numpy array (or a float or int scalar)
            The values of the wavelengths in m.

        """
        self.rays[:,10] =  2*numpy.pi/(wavelength * 1e2)

    #
    # info
    #
    def info(self):
        """
        Returns string containing some information of the beam (min, max, center, stDev) for
        several columns (1:6,26). It also gives the mean photon energy and intensity.

        Returns
        -------
        str
            The text.

        """
        try:
            txt = "col   name         min         max      center       stDev\n"
            for col in [1,2,3,4,5,6,26]:
                val = self.get_column(column=col,nolost=True)
                txt += "%3d  %5s  %10g  %10g  %10g  %10g\n" % \
                  (col, self.column_short_names()[col-1], val.min(), val.max(), val.mean(), val.std())
            txt += "\n"
            txt += "Number of rays: %d \n"%(self.get_number_of_rays())
            txt += "Number of good rays: %d \n"%(self.get_number_of_rays(nolost=1))
            txt += "Number of lost rays: %d \n"%(self.get_number_of_rays(nolost=2))
            txt += "Mean energy: %f eV\n"%(self.get_photon_energy_eV().mean() )
            if self.get_photon_energy_eV().mean() != 0.0:
                txt += "Mean wavelength: %f A\n"%(1e10 * self.get_photon_wavelength().mean() )
            txt += "Intensity (total): %f \n"%( self.get_intensity(nolost=1,polarization=0) )
            txt += "Intensity (s-pol): %f \n" % (self.get_intensity(nolost=1,polarization=1))
            txt += "Intensity (p-pol): %f \n" % (self.get_intensity(nolost=1,polarization=2))
        except:
            txt = ""
        return txt

    #
    # propagation / movements
    #
    def retrace(self, dist, resetY=False):
        """
        Propagates a beam at a given distance (in vacuum).

        Parameters
        ----------
        dist : float
            The distance to be re-traced.

        resetY : boolean, optional
            If True, reset all Y values to zero (which physically is like calculating the
            beam intersections with a plane perpendicular to the optical axis)

        """
        a0 = self.rays
        try:
            tof = (-a0[:,1] + dist)/a0[:,4]
            self.rays[:,0] += tof * self.rays[:,3]
            self.rays[:,1] += tof * self.rays[:,4]
            self.rays[:,2] += tof * self.rays[:,5]

            if resetY:
                self.rays[:,1] = 0.0
            #
            # TODO: modify optical path
            #
            self.rays[:,12] += tof
        except AttributeError:
            print ('shadow4.S4Beam.retrace: No rays')

    def translation(self, qdist1):
        """
        Translates spatially a beam by a given vector.

        Parameters
        ----------
        qdist1 : 3 elements list or tuple
            The distances to translate the X,Y and Z components.

        """
        if numpy.array(qdist1).size != 3:
            raise Exception("Input must be a vector [x,y,z]")

        self.rays[:,0] += qdist1[0]
        self.rays[:,1] += qdist1[1]
        self.rays[:,2] += qdist1[2]
        #
        # TODO: update optical path and may be phases of electric vectors
        #


    def rotate(self, theta, axis=1, rad=True):
        """
        Rotates a beam by a given angle along a given axis.

        Parameters
        ----------
        theta: float
            the rotation angle radians (of degress if rad=False).

        axis: int
            The axis number (Shadow's column) for the rotation (i.e, 1:x (default), 2:y, 3:z)

        rad: boolean, optional
            set False if theta is given in degrees.

        """
        if rad:
            theta1 = theta
        else:
            theta1 = theta * numpy.pi / 180

        a1 = self.rays.copy()

        if axis == 1:
            torot = [2,3]
        elif axis == 2:
            torot = [1,3]
        elif axis == 3:
            torot = [1,2]


        costh = numpy.cos(theta1)
        sinth = numpy.sin(theta1)

        tstart = numpy.array([1,4,7,16])

        for i in range(len(tstart)):

            newaxis = axis + tstart[i] - 1
            newaxisi = newaxis - 1
            newtorot = torot + tstart[i] - 1
            newtoroti = newtorot - 1

            self.rays[:,newtoroti[0]] =  a1[:,newtoroti[0]] * costh + a1[:,newtoroti[1]] * sinth
            self.rays[:,newtoroti[1]] = -a1[:,newtoroti[0]] * sinth + a1[:,newtoroti[1]] * costh
            self.rays[:,newaxisi]     =  a1[:,newaxisi]

    def change_to_image_reference_system(self, theta, T_IMAGE, rad=True,
                                         refraction_index=1.0,
                                         apply_attenuation=0,
                                         linear_attenuation_coefficient=0.0, # in SI, i.e. m^-1
                                         ):
        """
        Implements the propagation from the mirror reference frame to the screen (image) reference.
        Mimics IMREF and IMAGE1 subrutines in shadow3

        Parameters
        ----------
        theta: float
            the grazing angle in rad (default) or deg (if deg=False).

        T_IMAGE: float
            the distance o.e. to image in m.

        rad: boolean, optional
            set False if theta is given in degrees.

        """
        if rad:
            theta1 = theta
        else:
            theta1 = theta * numpy.pi / 180

        T_REFLECTION = numpy.pi / 2 - theta1

        a1 = self.rays.copy()

        nrays = self.rays.shape[0]
        UXIM_x = numpy.ones(nrays)
        UXIM_y = numpy.zeros(nrays)
        UXIM_z = numpy.zeros(nrays)

        VZIM_x = numpy.zeros(nrays)
        VZIM_y = numpy.zeros(nrays) - numpy.cos(T_REFLECTION)
        VZIM_z = numpy.zeros(nrays) + numpy.sin(T_REFLECTION)

        VNIMAG_x = numpy.zeros(nrays)
        VNIMAG_y = numpy.zeros(nrays) + numpy.sin(T_REFLECTION)
        VNIMAG_z = numpy.zeros(nrays) + numpy.cos(T_REFLECTION)

        ABOVE = T_IMAGE - a1[:,0] * VNIMAG_x - a1[:,1] * VNIMAG_y - a1[:,2] * VNIMAG_z
        BELOW = VNIMAG_x * a1[:,3] + VNIMAG_y * a1[:,4] + VNIMAG_z * a1[:,5]

        DIST = ABOVE / BELOW

        failure = numpy.argwhere(BELOW == 0)
        if len(failure) > 0:
            a1[failure, 9] = -3.0e-6

        # ! ** Computes now the intersections onto TRUE image plane.
        a1[:, 0]  +=   DIST * a1[:, 3]
        a1[:, 1]  +=   DIST * a1[:, 4]
        a1[:, 2]  +=   DIST * a1[:, 5]

        #!  ** Rotate now the results in the STAR (or TRUE image) reference plane.
        #!  ** Computes the projection of P_MIR onto the image plane versors.
        RIMCEN_x = VNIMAG_x * T_IMAGE
        RIMCEN_y = VNIMAG_y * T_IMAGE
        RIMCEN_z = VNIMAG_z * T_IMAGE

        a1[:, 0]  -=   RIMCEN_x
        a1[:, 1]  -=   RIMCEN_y
        a1[:, 2]  -=   RIMCEN_z

        #! ** Computes now the new vectors for the beam in the U,V,N ref.
        for i in [1,4,7,16]: # position, direction, Es, Ep
            # dot product
            self.rays[:, i - 1 + 0] = a1[:, i - 1 + 0] * UXIM_x   + a1[:, i - 1 + 1] * UXIM_y   + a1[:, i - 1 + 2] * UXIM_z
            self.rays[:, i - 1 + 1] = a1[:, i - 1 + 0] * VNIMAG_x + a1[:, i - 1 + 1] * VNIMAG_y + a1[:, i - 1 + 2] * VNIMAG_z
            self.rays[:, i - 1 + 2] = a1[:, i - 1 + 0] * VZIM_x   + a1[:, i - 1 + 1] * VZIM_y   + a1[:, i - 1 + 2] * VZIM_z

        # optical path col 13
        self.rays[:, 12] += numpy.abs(DIST) * refraction_index

        if apply_attenuation:
            att1 = numpy.sqrt(numpy.exp(-numpy.abs(DIST) * linear_attenuation_coefficient))
            self.rays[:, 7 - 1 ] *= att1
            self.rays[:, 8 - 1 ] *= att1
            self.rays[:, 9 - 1 ] *= att1
            self.rays[:, 16 - 1] *= att1
            self.rays[:, 17 - 1] *= att1
            self.rays[:, 18 - 1] *= att1

    @classmethod
    def get_UVW(self, X_ROT=0, Y_ROT=0, Z_ROT=0): # in radians!!
        """
        returns the rotation matrix given resulting from sequential rotations around the 3 axes X,Y,Z [todo: verify the order]

        Parameters
        ----------
        X_ROT : float
            rotation angle in rad around the X axis.

        Y_ROT : float
            rotation angle in rad around the Y axis.

        Z_ROT : float
            rotation angle in rad around the Z axis.

        Returns
        -------
        tuple
            the 9 matrix elements.

        """
        COSX =  numpy.cos(X_ROT)
        SINX = -numpy.sin(X_ROT)
        COSY =  numpy.cos(Y_ROT)
        SINY = -numpy.sin(Y_ROT)
        COSZ =  numpy.cos(Z_ROT)
        SINZ = -numpy.sin(Z_ROT)
        # ! C
        # ! C Computes the rotation matrix coefficients
        # ! C
        U_MIR_1 =  COSZ * COSY
        V_MIR_1 =  COSZ * SINX * SINY - SINZ * COSX
        W_MIR_1 =  COSZ * SINY * COSX + SINZ * SINX
        U_MIR_2 =  SINZ * COSY
        V_MIR_2 =  SINZ * SINX * SINY + COSZ * COSX
        W_MIR_2 =  SINZ * SINY * COSX - SINX * COSZ
        U_MIR_3 = -SINY
        V_MIR_3 =  COSY * SINX
        W_MIR_3 =  COSY * COSX

        return U_MIR_1, U_MIR_2, U_MIR_3,\
               V_MIR_1, V_MIR_2, V_MIR_3,\
               W_MIR_1, W_MIR_2, W_MIR_3,

    def rot_for(self, OFFX=0, OFFY=0, OFFZ=0, X_ROT=0, Y_ROT=0, Z_ROT=0):
        """
        Applies the roto-translation of the optical movement movements to the beam.

        Parameters
        ----------
        OFFX : float
            translation distance in m along the X axis.

        OFFX : float
            translation distance in m along the Y axis.

        OFFX : float
            translation distance in m along the Z axis.

        X_ROT : float
            rotation angle in rad around the X axis.

        Y_ROT : float
            rotation angle in rad around the Y axis.

        Z_ROT : float
            rotation angle in rad around the Z axis.

        """
        # ! C+++
        # ! C	SUBROUTINE	ROT_FOR
        # ! C
        # ! C	PURPOSE		Applies the roto-translation of the mirror movements
        # ! C			to the beam. This allows a complete decoupling of the system.
        # ! C
        # ! C	ARGUMENT	[ I ]	RAY	: the beam, as computed by RESTART.
        # ! C			    [ O ] 	RAY	: the beam, as seen by a MOVED mirror.
        # ! C
        # ! C---
        P_IN_1 = self.rays[:, 1-1].copy()
        P_IN_2 = self.rays[:, 2-1].copy()
        P_IN_3 = self.rays[:, 3-1].copy()
        V_IN_1 = self.rays[:, 4-1].copy()
        V_IN_2 = self.rays[:, 5-1].copy()
        V_IN_3 = self.rays[:, 6-1].copy()
        A_IN_1 = self.rays[:, 7-1].copy()
        A_IN_2 = self.rays[:, 8-1].copy()
        A_IN_3 = self.rays[:, 9-1].copy()

        # the P component is not implemented in shadow3. Buggy there?
        AP_IN_1 = self.rays[:, 16-1].copy()
        AP_IN_2 = self.rays[:, 17-1].copy()
        AP_IN_3 = self.rays[:, 18-1].copy()

        U_MIR_1, U_MIR_2, U_MIR_3, V_MIR_1, V_MIR_2, V_MIR_3, W_MIR_1, W_MIR_2, W_MIR_3 = \
            self.get_UVW(X_ROT=X_ROT, Y_ROT=Y_ROT, Z_ROT=Z_ROT)

        P_OUT_1 = (P_IN_1 - OFFX) * U_MIR_1 + (P_IN_2 - OFFY) * U_MIR_2 + (P_IN_3 - OFFZ) * U_MIR_3

        P_OUT_2 = (P_IN_1 - OFFX) * V_MIR_1 + (P_IN_2 - OFFY) * V_MIR_2 + (P_IN_3 - OFFZ) * V_MIR_3

        P_OUT_3 = (P_IN_1 - OFFX) * W_MIR_1 + (P_IN_2 - OFFY) * W_MIR_2 + (P_IN_3 - OFFZ) * W_MIR_3

        V_OUT_1 = V_IN_1 * U_MIR_1 + V_IN_2 * U_MIR_2 + V_IN_3 * U_MIR_3

        V_OUT_2 = V_IN_1 * V_MIR_1 + V_IN_2 * V_MIR_2 + V_IN_3 * V_MIR_3

        V_OUT_3 = V_IN_1 * W_MIR_1 + V_IN_2 * W_MIR_2 + V_IN_3 * W_MIR_3

        A_OUT_1 = A_IN_1 * U_MIR_1 + A_IN_2 * U_MIR_2 + A_IN_3 * U_MIR_3
        AP_OUT_1 = AP_IN_1 * U_MIR_1 + AP_IN_2 * U_MIR_2 + AP_IN_3 * U_MIR_3

        A_OUT_2 = A_IN_1 * V_MIR_1 + A_IN_2 * V_MIR_2 + A_IN_3 * V_MIR_3
        AP_OUT_2 = AP_IN_1 * V_MIR_1 + AP_IN_2 * V_MIR_2 + AP_IN_3 * V_MIR_3

        A_OUT_3 = A_IN_1 * W_MIR_1 + A_IN_2 * W_MIR_2 + A_IN_3 * W_MIR_3
        AP_OUT_3 = AP_IN_1 * W_MIR_1 + AP_IN_2 * W_MIR_2 + AP_IN_3 * W_MIR_3

        self.rays[:, 1-1 ] = P_OUT_1
        self.rays[:, 2-1 ] = P_OUT_2
        self.rays[:, 3-1 ] = P_OUT_3
        self.rays[:, 4-1 ] = V_OUT_1
        self.rays[:, 5-1 ] = V_OUT_2
        self.rays[:, 6-1 ] = V_OUT_3
        self.rays[:, 7-1 ] = A_OUT_1
        self.rays[:, 8-1 ] = A_OUT_2
        self.rays[:, 9-1 ] = A_OUT_3
        self.rays[:, 16-1] = AP_OUT_1
        self.rays[:, 17-1] = AP_OUT_2
        self.rays[:, 18-1] = AP_OUT_3

    def rot_back(self, OFFX=0, OFFY=0, OFFZ=0, X_ROT=0, Y_ROT=0, Z_ROT=0):
        """
        Applies the roto-translation of the optical movement movements to the beam.
        This will bring back the beam in the normal optical-element frame.


        Parameters
        ----------
        OFFX : float
            translation distance in m along the X axis.

        OFFX : float
            translation distance in m along the Y axis.

        OFFX : float
            translation distance in m along the Z axis.

        X_ROT : float
            rotation angle in rad around the X axis.

        Y_ROT : float
            rotation angle in rad around the Y axis.

        Z_ROT : float
            rotation angle in rad around the Z axis.


        """
        # ! C+++
        # ! C	SUBROUTINE	ROT_BACK
        # ! C
        # ! C	PURPOSE		Applies the roto-translation of the mirror movements
        # ! C			to the beam. This will bring bak the beam in the normal MIRROR frame.
        # ! C
        # ! C	ARGUMENT	[ I ]	RAY	: the beam, as computed by MIRROR.
        # ! C               [ O ] 	RAY	: the beam, as seen back in the mirror refernece frame.
        # ! C
        # ! C---
        P_IN_1  = self.rays[:, 1-1 ].copy()
        P_IN_2  = self.rays[:, 2-1 ].copy()
        P_IN_3  = self.rays[:, 3-1 ].copy()
        V_IN_1  = self.rays[:, 4-1 ].copy()
        V_IN_2  = self.rays[:, 5-1 ].copy()
        V_IN_3  = self.rays[:, 6-1 ].copy()
        A_IN_1  = self.rays[:, 7-1 ].copy()
        A_IN_2  = self.rays[:, 8-1 ].copy()
        A_IN_3  = self.rays[:, 9-1 ].copy()
        AP_IN_1 = self.rays[:, 16-1].copy()
        AP_IN_2 = self.rays[:, 17-1].copy()
        AP_IN_3 = self.rays[:, 18-1].copy()


        U_MIR_1, U_MIR_2, U_MIR_3, V_MIR_1, V_MIR_2, V_MIR_3, W_MIR_1, W_MIR_2, W_MIR_3 = \
            self.get_UVW(X_ROT=X_ROT, Y_ROT=Y_ROT, Z_ROT=Z_ROT)

        P_OUT_1 = P_IN_1 * U_MIR_1 + P_IN_2 * V_MIR_1 + P_IN_3 * W_MIR_1 + OFFX

        P_OUT_2 = P_IN_1 * U_MIR_2 + P_IN_2 * V_MIR_2 + P_IN_3 * W_MIR_2 + OFFY

        P_OUT_3 = P_IN_1 * U_MIR_3 + P_IN_2 * V_MIR_3 + P_IN_3 * W_MIR_3 + OFFZ

        V_OUT_1 = V_IN_1 * U_MIR_1 + V_IN_2 * V_MIR_1 + V_IN_3 * W_MIR_1

        V_OUT_2 = V_IN_1 * U_MIR_2 + V_IN_2 * V_MIR_2 + V_IN_3 * W_MIR_2

        V_OUT_3 = V_IN_1 * U_MIR_3 + V_IN_2 * V_MIR_3 + V_IN_3 * W_MIR_3

        A_OUT_1 = A_IN_1 * U_MIR_1 + A_IN_2 * V_MIR_1 + A_IN_3 * W_MIR_1
        AP_OUT_1 = AP_IN_1 * U_MIR_1 + AP_IN_2 * V_MIR_1 + AP_IN_3 * W_MIR_1

        A_OUT_2 = A_IN_1 * U_MIR_2 + A_IN_2 * V_MIR_2 + A_IN_3 * W_MIR_2
        AP_OUT_2 = AP_IN_1 * U_MIR_2 + AP_IN_2 * V_MIR_2 + AP_IN_3 * W_MIR_2

        A_OUT_3 = A_IN_1 * U_MIR_3 + A_IN_2 * V_MIR_3 + A_IN_3 * W_MIR_3
        AP_OUT_3 = AP_IN_1 * U_MIR_3 + AP_IN_2 * V_MIR_3 + AP_IN_3 * W_MIR_3

        self.rays[:, 1-1 ] = P_OUT_1
        self.rays[:, 2-1 ] = P_OUT_2
        self.rays[:, 3-1 ] = P_OUT_3
        self.rays[:, 4-1 ] = V_OUT_1
        self.rays[:, 5-1 ] = V_OUT_2
        self.rays[:, 6-1 ] = V_OUT_3
        self.rays[:, 7-1 ] = A_OUT_1
        self.rays[:, 8-1 ] = A_OUT_2
        self.rays[:, 9-1 ] = A_OUT_3
        self.rays[:, 16-1] = AP_OUT_1
        self.rays[:, 17-1] = AP_OUT_2
        self.rays[:, 18-1] = AP_OUT_3


    #
    # crop & boundaries
    #
    def crop_rectangle(self, x_col, x_min, x_max, y_col, y_min, y_max, negative=False, flag_lost_value=-1):
        """
        Crops the beam to a rectangle along the axes x_col and y_col.

        Parameters
        ----------
        x_col : int
            column 1 of the cropping surface.

        x_min : float
            minimum for x_col.

        x_max : float
            maximum for x_col.

        y_col : int
            column 2 of the cropping surface.

        y_min : float
            minimum for y_col.

        y_max : float
            maximum for y_col.

        negative : boolean, optional
            makes the negative (i.e. obstruction instead of aperture).

        flag_lost_value : float, optional
            the value to be assigned to the flagged bad rays.

        """
        x = self.get_column(x_col)
        y = self.get_column(y_col)
        flag = self.get_column(10)        # numpy.array(a3.getshonecol(10))

        TESTX = (x - x_min) * (x_max - x)
        TESTY = (y - y_min) * (y_max - y)

        indices_out = numpy.where(numpy.logical_or(TESTX < 0.0, TESTY < 0.0))

        if negative:
            window = numpy.zeros_like(flag)
            if len(indices_out) > 0: window[indices_out] = 1
        else:
            window = numpy.ones_like(flag)
            if len(indices_out) > 0: window[indices_out] = 0

        flag[window < 1] = flag_lost_value
        self.rays[:, 9] = flag

        return window

    def crop_ellipse(self, x_col, a1, a2, y_col, b1, b2, negative=False, flag_lost_value=-1):
        """
        Crops the beam to an ellipse along the axes x_col and y_col.

        Parameters
        ----------
        x_col : int
            column 1 of the cropping surface.

        a1 : float
            minimum for x_col.

        a2 : float
            maximum for x_col.

        y_col : int
            column 2 of the cropping surface.

        b1 : float
            minimum for y_col.

        b2 : float
            maximum for y_col.

        negative : boolean, optional
            makes the negative (i.e. obstruction instead of aperture).

        flag_lost_value : float, optional
            the value to be assigned to the flagged bad rays.

        """
        x =   self.get_column(x_col)
        y = self.get_column(y_col)
        flag = self.get_column(10)

        a0 = 0.5 * (a1 + a2)
        b0 = 0.5 * (b1 + b2)
        a11 = a2 - a1
        b11 = b2 - b1

        TESTX = (x - a0) ** 2 / (a11 / 2) ** 2 + (y - b0) ** 2 / (b11 / 2) ** 2 - 1.0
        TESTX = - TESTX
        TESTY = TESTX
        indices_out = numpy.where(numpy.logical_or(TESTX < 0.0, TESTY < 0.0))

        if negative:
            window = numpy.zeros_like(flag)
            if len(indices_out) > 0: window[indices_out] = 1
        else:
            window = numpy.ones_like(flag)
            if len(indices_out) > 0: window[indices_out] = 0

        flag[window < 1] = flag_lost_value
        self.rays[:, 9] = flag

        return window

    def crop_ellipse_with_hole(self, x_col, a1, a2, a3, a4,
                                    y_col, b1, b2, b3, b4, negative=False, flag_lost_value=-1):
        """
        Crops the beam to an elliptical annulus along the axes x_col and y_col.

        Parameters
        ----------
        x_col : int
            column 1 of the cropping surface.

        a1 : float
            inner ellipse minimum for x_col.

        a2 : float
            inner ellipse maximum for x_col.

        a3 : float
            outer ellipse minimum for x_col.

        a4 : float
            outer ellipse maximum for x_col.

        y_col : int
            column 2 of the cropping surface.

        b1 : float
            inner ellipse minimum for y_col.

        b2 : float
            inner ellipse maximum for y_col.

        b3 : float
            outer ellipse minimum for y_col.

        b4 : float
            outer ellipse maximum for y_col.

        negative : boolean, optional
            makes the negative (i.e. obstruction instead of aperture).

        flag_lost_value : float, optional
            the value to be asigned to the flagged bad rays.

        """
        x =   self.get_column(x_col)
        y = self.get_column(y_col)
        flag = self.get_column(10)

        a0 = 0.5 * (a1 + a2)
        b0 = 0.5 * (b1 + b2)
        a11 = a2 - a1
        b11 = b2 - b1

        A0 = 0.5 * (a3 + a4)
        B0 = 0.5 * (b3 + b4)
        A11 = a4 - a3
        B11 = b4 - b3

        TESTX = (x - A0) ** 2 / (A11 / 2) ** 2 + (y - B0) ** 2 / (B11 / 2) ** 2 - 1.0
        TESTX = - TESTX
        TESTY = (x - a0) ** 2 / (a11 / 2) ** 2 + (y - b0) ** 2 / (b11 / 2) ** 2 - 1.0
        indices_out = numpy.where(numpy.logical_or(TESTX < 0.0, TESTY < 0.0))

        if negative:
            window = numpy.zeros_like(flag)
            if len(indices_out) > 0: window[indices_out] = 1
        else:
            window = numpy.ones_like(flag)
            if len(indices_out) > 0: window[indices_out] = 0

        flag[window < 1] = flag_lost_value
        self.rays[:, 9] = flag

        return window


    def apply_boundaries_syned(self, syned_boundary_object, flag_lost_value=-1):
        """
        Crops the beam to a shape defined in a syned object.

        Parameters
        ----------
        syned_boundary_object : instance of syned.beamline.shape.Shape
            The cropping shape.


        flag_lost_value : float, optional
            the value to be assigned to the flagged bad rays.

        See Also
        --------
        syned.beamline.shape.Shape

        --
        """
        if isinstance(syned_boundary_object, type(None)):
            return
        elif isinstance(syned_boundary_object, Rectangle):
            x_left, x_right, y_bottom, y_top = syned_boundary_object.get_boundaries()
            self.crop_rectangle(1, x_left, x_right, 2, y_bottom, y_top, flag_lost_value=flag_lost_value)
        elif isinstance(syned_boundary_object, Ellipse):
            a_axis_min, a_axis_max, b_axis_min, b_axis_max = syned_boundary_object.get_boundaries()
            self.crop_ellipse(1, a_axis_min, a_axis_max, 2, b_axis_min, b_axis_max, flag_lost_value=flag_lost_value)
        elif isinstance(syned_boundary_object, Circle):
            radius, x_center, y_center = syned_boundary_object.get_boundaries()

            a_axis_min = x_center - radius
            a_axis_max = x_center + radius
            b_axis_min = y_center - radius
            b_axis_max = y_center + radius

            self.crop_ellipse(1, a_axis_min, a_axis_max, 2, b_axis_min, b_axis_max, flag_lost_value=flag_lost_value)
        elif isinstance(syned_boundary_object, TwoEllipses):
            a1_axis_min, a1_axis_max, b1_axis_min, b1_axis_max, \
                a2_axis_min, a2_axis_max, b2_axis_min, b2_axis_max = syned_boundary_object.get_boundaries()
            self.crop_ellipse_with_hole(1, a1_axis_min, a1_axis_max, a2_axis_min, a2_axis_max,
                                        2, b1_axis_min, b1_axis_max, b2_axis_min, b2_axis_max,
                                        flag_lost_value=flag_lost_value)
        else:
            raise Exception("Not good mirror boundary")

    def apply_boundaries_shadow(self, fhit_c=0, fshape=1, rlen1=0.0, rlen2=0.0, rwidx1=0.0, rwidx2=0.0,
                                flag_lost_value=-1):
        """
        Apply boundaries using the shadow3 flags and variables.

        Parameters
        ----------
        fhit_c: int, optional
            flag: mirror dimensions finite: yes (1), no(0).

        fshape: int, optional
            for fhit_c=1: mirror shape rectangular (1)
                full ellipse (2)
                ellipse with hole (3).

        rlen1: float, optional
            fshape=1: mirror half length +Y.
            fshape=3: internal minor axis (Y).

        rlen2: float, optional
            fshape=1: mirror half length -Y.
            fshape=2,3: external outline minor

        rwidx1: float, optional
            fshape=1: mirror half width +X.
            fshape=3: internal major axis (X).

        rwidx2: float, optional
            fshape=1: mirror half width -X.
            fshape=2,3: external outline major axis (X).

        flag_lost_value : float, optional
            the value to be assigned to the flagged bad rays.

        """
        if fhit_c == 0:
            return

        x =   self.get_column(1)
        y = self.get_column(2)
        flag = self.get_column(10)        # numpy.array(a3.getshonecol(10))

        if fshape == 1:  # rectangle
            x_min = -rwidx2
            x_max =  rwidx1
            y_min = -rlen2
            y_max =  rlen1
            self.crop_rectangle(1, x_min, x_max, 2, y_min, y_max, negative=False, flag_lost_value=flag_lost_value)
        elif fshape == 2: # ellipse
            self.crop_ellipse(1, -rwidx2/2, rwidx2/2, 2, -rlen2/2, rlen2/2, negative=False, flag_lost_value=flag_lost_value)
        elif fshape == 3:  # hole in ellipse
            raise Exception("Not yet implemented")


    def apply_reflectivity_s(self, Rs):
        """
        Apply sigma-reflectivity to the beam.

        Parameters
        ----------
        Rs : float
            The reflectivity value (real number, if complex, use apply_complex_reflectivity_s()).

        """
        if numpy.iscomplexobj(Rs):
            print(">>> Warning: using complex reflectivities. Use apply_complex_reflectivity_s() instead")
        self.rays[:, 6] *= Rs
        self.rays[:, 7] *= Rs
        self.rays[:, 8] *= Rs

    def apply_reflectivity_p(self, Rp):
        """
        Apply pi-reflectivity to the beam.

        Parameters
        ----------
        Rp : float
            The reflectivity value (real number, if complex, use apply_complex_reflectivity_p()).

        """
        if numpy.iscomplexobj(Rp):
            print(">>> Warning: using complex reflectivities. Use apply_complex_reflectivity_p() instead")

        self.rays[:, 15] *= Rp
        self.rays[:, 16] *= Rp
        self.rays[:, 17] *= Rp

    def apply_reflectivities(self, Rs, Rp):
        """
        Apply sigma- and pi- reflectivities to the beam.

        Parameters
        ----------
        Rs : float
            The reflectivity value (real number, if complex, use apply_complex_reflectivities()).

        Rp : float
            The reflectivity value (real number, if complex, use apply_complex_reflectivities()).

        """
        self.apply_reflectivity_s(Rs)
        self.apply_reflectivity_p(Rp)

    # complex reflectivities...
    def apply_complex_reflectivity_s(self, Rs):
        """
        Apply sigma-reflectivity to the beam.

        Parameters
        ----------
        Rs : complex
            The reflectivity value.

        """
        self.rays[:, 6] *= numpy.abs(Rs)
        self.rays[:, 7] *= numpy.abs(Rs)
        self.rays[:, 8] *= numpy.abs(Rs)
        self.rays[:, 13] += numpy.angle(Rs)

    def apply_complex_reflectivity_p(self, Rp):
        """
        Apply pi-reflectivity to the beam.

        Parameters
        ----------
        Rp : complex
            The reflectivity value.

        """
        self.rays[:, 15] *= numpy.abs(Rp)
        self.rays[:, 16] *= numpy.abs(Rp)
        self.rays[:, 17] *= numpy.abs(Rp)
        self.rays[:, 14] += numpy.angle(Rp)

    def apply_complex_reflectivities(self, Rs, Rp):
        """
        Apply sigma- and pi- reflectivities to the beam.

        Parameters
        ----------
        Rs : complex
            The reflectivity value.

        Rp : complex
            The reflectivity value.

        """
        self.apply_complex_reflectivity_s(Rs)
        self.apply_complex_reflectivity_p(Rp)

    # phases
    def add_phase_s(self, phase):
        """
        Add a sigma-phase to the beam.

        Parameters
        ----------
        phase : float
            phase angle in rad.

        """
        self.rays[:, 13] += phase

    def add_phase_p(self, phase):
        """
        Add a pi-phase to the beam.

        Parameters
        ----------
        phase : float
            phase angle in rad.

        """
        self.rays[:, 14] += phase

    def add_phases(self, phase_s, phase_p):
        """
        Add a sigma and pi phases to the beam.

        Parameters
        ----------
        phase_s : float
            The sigma phase in rad.

        phase_p : float
            the pi phase in rad.

        """
        self.add_phase_s(phase_s)
        self.add_phase_p(phase_p)

    #
    #  interfaces like in shadow3
    #
    def generate_source(self, source_object): #todo: remove?
        """
        Put rays that sample a given source_object. This option mimics the obsolete "source" in shadow3.
        It should not be used, only for compatibility purposes.

        Parameters
        ----------
        source_object : instance of a shadow4.source
            Instance of a source class that implements the get_rays() method.

        Returns
        -------
        S4Beam instance
            the new source. Note that the calling beam is also modified.

        """
        try:
            tmp = source_object.get_beam()
            self.rays = tmp.get_rays()
            return tmp
        except:
            raise Exception("shadow4 source class must implement get_rays method")

    def trace_oe(self, oe_object, n, overwrite=True): #todo: remove?
        """
        Modify rays to trace an optical element. This option mimics the obsolete "trace" in shadow3.
        It should not be used, only for compatibility purposes.

        Parameters
        ----------
        oe_object : instance of a shadow4 optical element.
            Instance of a source class that implements the trace_beam() method.

        n : int
            The element number (not used here)

        """
        try:
            beam_result = oe_object.trace_beam(self)
            if overwrite:
                self.rays = beam_result.rays
            return beam_result
        except:
            raise Exception("shadow4 class used a optical element must implement the trace_beam method")
    #
    # useful tools (labels)
    #

    @classmethod
    def column_names(cls):
        """
        returns a list with the names of shadow4 beam columns (the 18 columns in the beam plus the extended columns).

        Returns
        -------
        list
            The column names.

        """
        return [
            "X spatial coordinate",
            "Y spatial coordinate",
            "Z spatial coordinate",
            "X' direction or divergence",
            "Y' direction or divergence",
            "Z' direction or divergence",
            "X component of the electromagnetic vector (s-polariz)",
            "Y component of the electromagnetic vector (s-polariz)",
            "Z component of the electromagnetic vector (s-polariz)",
            "Lost ray flag",
            "Wavenumber",
            "Ray Index",
            "Optical path length",
            "\u03C6\u209B Phase (s-polarization)",
            "\u03C6\u209A Phase (p-polarization)",
            "X component of the electromagnetic vector (p-polariz)",
            "Y component of the electromagnetic vector (p-polariz)",
            "Z component of the electromagnetic vector (p-polariz)",
            "Wavelength",
            "R = \u221A(X\u00B2+Y\u00B2+Z\u00B2)",
            "Angle from Y axis",
            "Magnitude of the Electromagnetic vector",
            "|E\u209B|\u00B2 + |E\u209A|\u00B2 (total intensity)",
            "|E\u209B|\u00B2 (total intensity for s-polarization)",
            "|E\u209A|\u00B2 (total intensity for p-polarization)",
            "Photon Energy",
            "Kx = 2\u03C0 / \u03BB * X'",
            "Ky = 2\u03C0 / \u03BB * Y'",
            "Kz = 2\u03C0 / \u03BB * Z'",
            "S0-stokes = |E\u209B|\u00B2 + |E\u209A|\u00B2",
            "S1-stokes = |E\u209B|\u00B2 - |E\u209A|\u00B2",
            "S2-stokes = 2 |E\u209B| |E\u209A| cos(\u03C6\u209B-\u03C6\u209A)",
            "S3-stokes = 2 |E\u209B| |E\u209A| sin(\u03C6\u209B-\u03C6\u209A)",
            "Power = Intensity(col 23) * Energy (col 26)",
            "Angle-X with Y = |arcsin(X')|",
            "Angle-Z with Y = |arcsin(Z')|",
            "Angle-X with Y = |arcsin(X') - mean(arcsin(X'))|",
            "Angle-Z with Y = |arcsin(Z') - mean(arcsin(Z'))|",
            "Phase difference \u03C6\u209B-\u03C6\u209A",
            "Complex electric field (s-polariz)",
            "Complex electric field (p-polariz)",
        ]

    @classmethod
    def column_units(cls):
        """
        returns a list with the column units of shadow4 beam columns (the 18 columns in the beam plus the extended columns).

        Returns
        -------
        list
            The column units.

        """
        return [
            "[m]",
            "[m]",
            "[m]",
            "[rad]",
            "[rad]",
            "[rad]",
            "",
            "",
            "",
            "",
            "[cm\u207B\u00B9]",
            "",
            "[m]",
            "[rad]",
            "[rad]",
            "",
            "",
            "",
            "[\u00C5]",
            "[m]",
            "[rad]",
            "",
            "",
            "",
            "",
            "[eV]",
            "[\u00C5\u207B\u00B9]",
            "[\u00C5\u207B\u00B9]",
            "[\u00C5\u207B\u00B9]",
            "",
            "",
            "",
            "",
            "",
            "[rad]",
            "[rad]",
            "[rad]",
            "[rad]",
            "[rad]",
            "",
            "",
        ]


    @classmethod
    def column_short_names(cls):
        """
        returns a list with the short-names of shadow4 beam columns (the 18 columns in the beam plus the extended columns).
        Useful for labeling plots.

        Returns
        -------
        list
            The column short-names.

        """
        return [
                "X","Y","Z",
                "X'", "Z'", "Z'",
                "Ex\u209B", "Ey\u209B", "Ez\u209B",
                "Flag","K","Idx","Opt. Path",
                "Phase\u209B","Phase\u209A",
                "Ex\u209A", "Ey\u209A", "Ez\u209A",
                "\u03BB",
                "R",
                "Angle from Y",
                "|E|",
                "I tot",
                "I s-pol",
                "I p-pol",
                "Energy",
                "Kx",
                "Ky",
                "Kz",
                "S0",
                "S1",
                "S2",
                "S3",
                "Power",
                "Angle X^Y",
                "Angle Z^Y",
                "Angle X^Y",
                "Angle Z^Y",
                "\u03C6\u209B-\u03C6\u209A", # phase_s - phase_p
                "complex E\u209B",
                "complex E\u209A",
        ]

    @classmethod
    def column_names_with_column_number(cls):
        """
        returns a list with the names of shadow4 beam columns including the column number (the 18 columns in the beam plus the extended columns).

        Returns
        -------
        list
            The column names.

        """
        names = cls.column_names()
        for i in range(len(names)): names[i] = "%2d: %s" % (i+1, names[i])
        return names

    @classmethod
    def column_short_names_with_column_number(cls):
        """
        returns a list with the short-names of shadow4 beam columns including the column number (the 18 columns in the beam plus the extended columns).

        Returns
        -------
        list
            The column names.

        """
        names = cls.column_short_names()
        for i in range(len(names)): names[i] = "%2d: %s" % (i+1, names[i])
        return names

    @classmethod
    def column_names_formatted(cls):
        return [
            r"$x$ [user unit]", #"x [user unit]"],
            r"$y$ [user unit]", #"y [user unit]"],
            r"$z$ [user unit]", #"z [user unit]"],
            r"$\dot{x}$ [rads]", #"x' [rads]"],
            r"$\dot{y}$ [rads]", #"y' [rads]"],
            r"$\dot{z}$ [rads]", #"z' [rads]"],
            r"$\mathbf{E}_{\sigma x}$", #"Es_x"],
            r"$\mathbf{E}_{\sigma y}$", #"Es_y"],
            r"$\mathbf{E}_{\sigma z}$", #"Es_z"],
            r"ray flag", #"Ray flag"],
            r"$K = \frac{2 \pi}{\lambda} [A^{-1}]$", #"K magnitude"], ## [r"E [eV]", "Energy"],
            r"Ray index", #"Ray index"],
            r"s", #"Opt. Path"],
            r"$\phi_{\sigma}$", #"phase_s"],
            r"$\phi_{\pi}$", #"phase_p"],
            r"$\mathbf{E}_{\pi x}$", #"Ep_x"],
            r"$\mathbf{E}_{\pi y}$", #"Ep_y"],
            r"$\mathbf{E}_{\pi z}$", #"Ep_z"],
            r"$\lambda$ [$\AA$]", #"wavelength"],
            r"$R= \sqrt{x^2+y^2+z^2}$", #"R [user unit]"],
            r"$\theta$", #"theta"],
            r"$\Vert\mathbf{E_{\sigma}}+\mathbf{E_{\pi}}\Vert$", #"Electromagnetic vector magnitude"],
            r"$\Vert\mathbf{E_{\sigma}}\Vert^2+\Vert\mathbf{E_{\pi}}\Vert^2$", #"intensity (weight column = 23: |E|^2 (total intensity))"],
            r"$\Vert\mathbf{E_{\sigma}}\Vert^2$", #"intensity (weight column = 24: |E_s|^2 (sigma intensity))"],
            r"$\Vert\mathbf{E_{\pi}}\Vert^2$", #"intensity (weight column = 25: |E_p|^2 (pi intensity))"],
            r"E [eV]", #"Energy"], #       [r"$K = \frac{2 \pi}{\lambda} [A^{-1}]$", "K magnitude"]
            r"$K_x = \frac{2 \pi}{\lambda} \dot{x}$ [$\AA^{-1}$]", #"K_x"],
            r"$K_y = \frac{2 \pi}{\lambda} \dot{y}$ [$\AA^{-1}$]", #"K_y"],
            r"$K_z = \frac{2 \pi}{\lambda} \dot{z}$ [$\AA^{-1}$]", #"K_z"],
            r"$S_0 = \Vert\mathbf{E}_{\sigma}\Vert^2 + \Vert\mathbf{E}_{\pi}\Vert^2 $", #"S0"],
            r"$S_1 = \Vert\mathbf{E}_{\sigma}\Vert^2 - \Vert\mathbf{E}_{\pi}\Vert^2 $", #"S1"],
            r"$S_2 = 2 \Vert\mathbf{E}_{\sigma}\Vert \cdot \Vert\mathbf{E}_{\pi}\Vert \cos{(\phi_{\sigma}-\phi_{\pi})}$", #"S2"],
            r"$S_3 = 2 \Vert\mathbf{E}_{\sigma}\Vert \cdot \Vert\mathbf{E}_{\pi}\Vert \sin{(\phi_{\sigma}-\phi_{\pi})}$", #"S3"],
            r"Power [eV/s]", #"Power", "Power",
            r"Angle X^Y [rad]", # "Angle X^Y [rad]"],
            r"Angle Z^Y [rad]", # "Angle Z^Y [rad]"],
            r"Angle X^Y [rad]", # "Angle X^Y [rad]"],
            r"Angle Z^Y [rad]", # "Angle Z^Y [rad]"],
            r"Phase difference \phi_{s}-\phi_{p}",
            r"complex $\mathbf{E}_{s}$",  # "Es_complex"],
            r"complex $\mathbf{E}_{p}$",  # "Es_complex"],
        ]

    @classmethod
    def column_names_formatted_with_column_number(cls):
        """
        returns a list with the formatted-names of shadow4 beam columns including the column number (the 18 columns in the beam plus the extended columns).

        Returns
        -------
        list
            The column names.

        """
        names = cls.column_names_formatted()
        for i in range(len(names)): names[i] = "%2d: %s" % (i+1, names[i])
        return names

    #
    # useful tools (h5 files)
    #
    def write_h5(self, filename, overwrite=True, simulation_name="run001", beam_name="begin"):
        """
        writes a beam in an h5 file.

        Parameters
        ----------
        filename : str
            file name.

        overwrite : boolean, optional
            if True, overwrite existing file (if False, the beam is appended in the existing file).

        simulation_name : str, optional
            a simulation name,

        beam_name : str, optional
            a beam name.

        """
        if overwrite:
            try:
                os.remove(filename)
            except:
                pass
            f = h5py.File(filename, 'w')
        else:
            f = h5py.File(filename, 'a')

        # header
        # point to the default data to be plotted
        f.attrs['default'] = 'entry'
        # give the HDF5 root some more attributes
        f.attrs['file_name'] = filename
        f.attrs['file_time'] = time.time()
        f.attrs['creator'] = "shadow4"
        f.attrs['HDF5_Version'] = h5py.version.hdf5_version
        f.attrs['h5py_version'] = h5py.version.version


        try:
            f1 = f[simulation_name]
        except:
            f1 = f.create_group(simulation_name)

        f1.attrs['NX_class'] = 'NXentry'
        f1.attrs['default'] = "begin"

        rays = self.get_rays()

        f2 = f1.create_group(beam_name)
        f2.attrs['NX_class'] = 'NXdata'
        f2.attrs['signal'] =  b'col03 z'
        f2.attrs['axes'] = b'col01 x'

        column_names = self.column_short_names_with_column_number()

        for i in range(18):
            column_name = column_names[i]
            # Y data
            ds = f2.create_dataset(column_name, data=rays[:,i].copy())
            ds.attrs['long_name'] = "column %s"%(i+1)  # suggested X axis plot label
        f.close()
        print("File written/updated: %s"%filename)

    @classmethod
    def load_h5(cls, filename, simulation_name="run001", beam_name="begin"):
        """
        loads a beam from an h5 file.

        Parameters
        ----------
        filename : str
            file name.

        simulation_name : str, optional
            a simulation name,

        beam_name : str, optional
            a beam name.

        Returns
        -------
        S4beam instance
            The beam

        """
        f = h5py.File(filename, 'r')
        column_names = cls.column_short_names_with_column_number()

        try:
            x = (f["%s/%s/col01 x"%(simulation_name,beam_name)])[:]

            beam = S4Beam(N=x.size)
            rays = numpy.zeros( (x.size,18))
            for i in range(18):
                column_name = column_names[i]
                rays[:,i] = (f["%s/%s/%s"%(simulation_name,beam_name,column_name)])[:]
        except:
            f.close()
            raise Exception("Cannot find data in %s:/%s/%s" % (filename, simulation_name, beam_name))

        f.close()
        return S4Beam.initialize_from_array(rays)


    #
    # useful tools (compare beams)
    #
    def identical(self, beam2):
        """
        Compares two beams

        Parameters
        ----------
        beam2 : S4 instance
            the beam to compare with.

        Returns
        -------
        boolean
            True if the two beams are almost equal.

        """
        try:
            assert_almost_equal(self.rays,beam2.rays)
            return True
        except:
            return False

    def difference(self, beam2):
        """
        Compares two beams and prints the diferences in different columns.

        Parameters
        ----------
        beam2 : S4 instance
            the beam to compare with.


        """
        raysnew = beam2.rays
        fact  = 1.0
        for i in range(18):
            m0 = (raysnew[:, i] * fact).mean()
            m1 = self.rays[:, i].mean()
            if numpy.abs(m1) > 1e-10:
                print("\ncol %d, mean: beam_tocheck %g, beam %g , diff/beam %g: " % (i + 1, m0, m1, numpy.abs(m0 - m1) / numpy.abs(m1)))
            else:
                print("\ncol %d, mean: beam_tocheck %g, beam %g " % (i + 1, m0, m1))

            std0 = (raysnew[:, i] * fact).std()
            std1 = self.rays[:, i].std()
            if numpy.abs(std1) > 1e-10:
                print("col %d, std: beam_tocheck %g, beam %g , diff/beam %g" % (i + 1, std0, std1, numpy.abs(std0 - std1) / numpy.abs(std1)))
            else:
                print("col %d, std: beam_tocheck %g, beam %g " % (i + 1, std0, std1))



if __name__ == "__main__":
    print(S4Beam.column_names_with_column_number())
    print(S4Beam.column_short_names_with_column_number())
    print(S4Beam.column_names_formatted())
    print(S4Beam.column_names_formatted_with_column_number())
    print(S4Beam.get_UVW())

    B = S4Beam.initialize_as_pencil(N=1)
    print(B.info())

    print("number of rays : ", B.N, B.get_number_of_rays())
    print("number of good rays : ", B.Ngood, B.get_number_of_rays(nolost=1))
    print("number of bad rays : ", B.Nbad, B.get_number_of_rays(nolost=2))
    B = S4Beam.initialize_as_pencil(N=100)
    print("all : ", (B.get_rays(nolost=0)).shape, B.rays.shape)
    print("good: ", (B.get_rays(nolost=1)).shape, B.rays_good.shape)
    print("bad : ", (B.get_rays(nolost=2)).shape, B.rays_bad.shape)
