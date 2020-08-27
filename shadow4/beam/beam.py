import numpy
import scipy.constants as codata
from numpy.testing import assert_almost_equal

import h5py
import time
import os

from shadow4.syned.shape import Rectangle, Ellipse, TwoEllipses
# TODO from syned.beamline.shape import Rectangle, Ellipse, TwoEllipses


# IMPORTANT: Column 11 (index 10) is wavenumber (cm^-1) as internally in Shadow

class Beam(object):

    def __init__(self, N=1000, array=None):
        """

        :param N:
        :param array:
        :return:
        """
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

        :param array:
        :return:
        """
        if array.shape[1] != 18:
            raise Exception("Bad array shape: must be (npoints,18)")
        return Beam(array=array)

    @classmethod
    def initialize_as_pencil(cls, N=1000):
        """

        :param array:
        :return:
        """
        beam = Beam(N)
        beam.set_column(5,1.0) # Vy
        beam.set_column(7,1.0) # Es
        beam.set_column(10,1.0) # flag
        beam.set_column(11,2*numpy.pi/1e-8) # wavenumber (1 A)
        beam.set_column(12,numpy.arange(N,dtype=float)) # index

        return beam


    def duplicate(self):
        return Beam.initialize_from_array(self.rays.copy())

    #
    # getters
    #

    def get_rays(self):
        """

        :return: numpy array (npoints,18)
        """
        return self.rays.copy()

    def get_number_of_rays(self,nolost=0):
        """

        :param nolost: flag (default=0)
        :return: number of rays
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

    def get_photon_energy_eV(self,nolost=0):
        """
        returns the array with photon energy

        :param nolost: 0: all rays  1: good rays, 2: bad rays
        :return: array
        """
        A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        return self.get_column(11,nolost=nolost) / A2EV

    def get_photon_wavelength(self):
        return 2*numpy.pi/self.get_column(11) * 1e-2

    def get_intensity(self,nolost=0,polarization=0):
        if polarization == 0:
            w = self.get_column(23,nolost=nolost)
        elif polarization == 1:
            w = self.get_column(24, nolost=nolost)
        elif polarization == 2:
            w = self.get_column(25, nolost=nolost)
        return w.sum()

    def get_column(self,column,nolost=0):
        """
            Possible choice for column are:
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

            -11: column 26
        :param column:
        :return: a copy of the arrays
        """

        if column == -11: column = 26

        if column <= 18:
            out = self.rays[:,column-1]
        else:
            A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
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
                f = self.getshonecol(10)
                w = self.getshonecol(23)
                xp = self.getshonecol(4)
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
                f = self.getshonecol(10)
                w = self.getshonecol(23)
                zp = self.getshonecol(6)
                if nolost == 1:
                    findices  = numpy.where(f > 0.0)
                    if len(findices[0])==0:
                        col_mean = numpy.average(zp, weights=w)
                    else:
                        col_mean = numpy.average(zp[findices], weights=w[findices])
                else:
                    col_mean = numpy.average(zp, weights=w)


                out = numpy.abs(numpy.arcsin(zp - col_mean))




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

    def get_columns(self,columns,nolost=0):
        ret = []
        if isinstance(columns, int): return self.get_column(columns,nolost=nolost)
        for c in columns:
            ret.append(self.get_column(c,nolost=nolost))
        return numpy.array(tuple(ret))



    def get_standard_deviation(self,col, nolost=1, ref=0):
        '''
        returns the standard deviation of one viariable in the beam
        :param col: variable (shadow column number)
        :param nolost: 0 = use all rays, 1=good only, 2= lost only
        :param ref: 0 = no weight, 1=weight with intensity (col23)
        :return:
        '''
        x = self.get_column(col,nolost=nolost)
        if ref == 0:
            return x.std()
        else:
            w = self.get_column(23,nolost=nolost)
            average = numpy.average(x, weights=w)
            variance = numpy.average( (x-average)**2, weights=w)
            return(numpy.sqrt(variance))

    def intensity(self,nolost=0):
        w = self.get_column(23,nolost=nolost)
        return w.sum()

    def get_good_range(self, icol, nolost=0):
        """

        :param icol: the column number (SHADOW convention, starting from 1)
        :param nolost: lost rays flag (0=all, 1=good, 2=losses)
        :return: [rmin,rmax] the selected range
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
        if rmin == 0.0:
            rmin = -1.0
            rmax = 1.0
        return [rmin, rmax]

    def histo2(self,col_h,col_v,nbins=25,ref=23, nbins_h=None, nbins_v=None, nolost=0,xrange=None,yrange=None,
             calculate_widths=1):
        """
        performs 2d histogram to prepare data for a plotxy plot
        It uses histogram2d for calculations
        Note that this Shadow.Beam.histo2 was previously called Shadow.Beam.plotxy
        :param col_h: the horizontal column
        :param col_v: the vertical column
        :param nbins: number of bins
        :param ref      :
                   0, None, "no", "NO" or "No":   only count the rays
                   23, "Yes", "YES" or "yes":     weight with intensity (look at col=23 |E|^2 total intensity)
                   other value: use that column as weight
        :param nbins_h: number of bins in H
        :param nbins_v: number of bins in V
        :param nolost: 0 or None: all rays, 1=good rays, 2=only losses
        :param xrange: range for H
        :param yrange: range for V
        :param calculate_widths: 0=No, 1=calculate FWHM (default), 2=Calculate FWHM and FW at 25% and 75% if Maximum
        :return: a dictionary with all data needed for plot
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
              print("Shadow.Beam.histo2: Warning: weighting with column 1 (X) [not with intensity as may happen in old versions]")

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

    def set_column(self,column,value):
        """
        :param column:
        :param value:
        :return:
        """
        self.rays[:,column-1] = value

    def set_photon_energy_eV(self,energy_eV):
        A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        self.rays[:,10] = energy_eV * A2EV

    def set_photon_wavelength(self,wavelength):
        self.rays[:,10] =  2*numpy.pi/(wavelength * 1e2)


    #
    # info
    #
    def info(self):
        """

        :return:
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

    def retrace(self,dist,resetY=False):
        """

        :param dist:
        :param resetY:
        :return:
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
            print ('Beam.retrace: No rays')



    def translation(self,qdist1):
        """

        :param qdist1: translation vector
        :return:

        """

        if numpy.array(qdist1).size != 3:
            raise Exception("Input must be a vector [x,y,z]")

        self.rays[:,0] += qdist1[0]
        self.rays[:,1] += qdist1[1]
        self.rays[:,2] += qdist1[2]

        #
        # TODO: update optical path and may be phases of electric vectors
        #


    def rotate(self,theta,axis=1,rad=1):
        """

        :param theta: the rotation angle in degrees (default=0)
        :param axis: The axis number (Shadow's column) for the rotation
                    (i.e, 1:x (default), 2:y, 3:z)
        :param file:
        :param rad: set this flag when theta1 is in radiants
        :return:
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



    #
    # crop
    #
    def crop_rectangle(self, x_col, x_min, x_max, y_col, y_min, y_max, negative=False, flag_lost_value=-1):

        x = self.get_column(x_col)
        y = self.get_column(y_col)
        flag = self.get_column(10)        # numpy.array(a3.getshonecol(10))

        TESTX = (x - x_min) * (x_max - x)
        TESTY = (y - y_min) * (y_max - y)

        indices_out = numpy.where(numpy.logical_or(TESTX < 0.0, TESTY < 0.0))

        # if not negative:
        #     window = numpy.ones_like(flag)
        #     lower_window_x = numpy.where(x < x_min)
        #     upper_window_x = numpy.where(x > x_max)
        #     lower_window_y = numpy.where(y < y_min)
        #     upper_window_y = numpy.where(y > y_max)
        #
        #     #
        #     if len(lower_window_x) > 0: window[lower_window_x] = 0
        #     if len(upper_window_x) > 0: window[upper_window_x] = 0
        #     if len(lower_window_y) > 0: window[lower_window_y] = 0
        #     if len(upper_window_y) > 0: window[upper_window_y] = 0
        #
        # else:
        #     window = numpy.ones_like(flag)
        #     window2 = numpy.ones_like(window)
        #     window_x = numpy.where((x_min <= x) & (x <= x_max))
        #     window_y = numpy.where((y_min <= y) & (y <= y_max))
        #
        #     if len(window_x) > 0: window[window_x] = 0.0
        #     if len(window_y) > 0: window2[window_y] = 0.0
        #
        #     window += window2
        #     window_good = numpy.where(window > 0)
        #     if len(window_good) > 0: window[window_good] = 1.0

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
        print(">>>>> apply_boundaries_syned: ", syned_boundary_object, Rectangle)
        if isinstance(syned_boundary_object, type(None)):
            return
        elif isinstance(syned_boundary_object, Rectangle):
            x_left, x_right, y_bottom, y_top = syned_boundary_object.get_boundaries()
            print(">>>>>>>>>>>>>>>>>>>>dimensions ",x_left, x_right, y_bottom, y_top)
            self.crop_rectangle(1, x_left, x_right, 2, y_bottom, y_top)
        elif isinstance(syned_boundary_object, Ellipse):
            a_axis_min, a_axis_max, b_axis_min, b_axis_max = syned_boundary_object.get_boundaries()
            self.crop_ellipse(1, a_axis_min, a_axis_max, 2, b_axis_min, b_axis_max)
        elif isinstance(syned_boundary_object, TwoEllipses):
            a1_axis_min, a1_axis_max, b1_axis_min, b1_axis_max, \
                a2_axis_min, a2_axis_max, b2_axis_min, b2_axis_max = syned_boundary_object.get_boundaries()
            self.crop_ellipse_with_hole(1, a1_axis_min, a1_axis_max, a2_axis_min, a2_axis_max,
                                        2, b1_axis_min, b1_axis_max, b2_axis_min, b2_axis_max,)
        else:
            raise Exception("Not good mirror boundary")

    def apply_boundaries_shadow(self, fhit_c=0, fshape=1, rlen1=0.0, rlen2=0.0, rwidx1=0.0, rwidx2=0.0, flag_lost_value=-1):

        """
        fhit_c	= 0 - flag: mirror dimensions finite: yes (1), no(0).

        fshape	= 0 - for fhit_c=1:
                mirror shape rectangular (1)
                full ellipse (2)
                ellipse with hole (3).

        rlen1	=  0.0
                fshape=1: mirror half length +Y.
                 fshape=3: internal minor axis (Y).
        rlen2	=  0.0
                fshape=1: mirror half length -Y.
                fshape=2,3: external outline minor

        rwidx1	=  0.0
                fshape=1: mirror half width +X.
                fshape=3: internal major axis (X).
        rwidx2	=  0.0
                fshape=1: mirror half width -X.
                fshape=2,3: external outline major axis (X).
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
        self.rays[:, 6] *= Rs
        self.rays[:, 7] *= Rs
        self.rays[:, 8] *= Rs

    def apply_reflectivity_p(self, Rp):
        self.rays[:, 15] *= Rp
        self.rays[:, 16] *= Rp
        self.rays[:, 17] *= Rp

    def apply_reflectivities(self, Rs, Rp):
        self.apply_reflectivity_s(Rs)
        self.apply_reflectivity_p(Rp)

    #
    # file i/o
    #

    # def get_shadow3_beam(self):
    #     #TODO this dump uses now shadow3. To be removed after checking or write using fully python
    #     import Shadow
    #     beam_shadow3 = Shadow.Beam(N=self.get_number_of_rays())
    #     beam_shadow3.rays = self.get_rays().copy()
    #     return beam_shadow3
    #
    #     # beam_shadow3.write(file)
    #     # print("File %s written to disk. "%file)
    #
    # def dump_shadow3_file(self,file):
    #     #TODO this dump uses now shadow3. To be removed after checking or write using fully python
    #     beam3 = self.get_shadow3_beam()
    #     beam3.write(file)
    #     print("File %s written to disk. "%file)


    #
    #  interfaces like in shadow3
    #



    #
    # file i/o
    #

    # def get_shadow3_beam(self):
    #     #TODO this dump uses now shadow3. To be removed after checking or write using fully python
    #     import Shadow
    #     beam_shadow3 = Shadow.Beam(N=self.get_number_of_rays())
    #     beam_shadow3.rays = self.get_rays().copy()
    #     return beam_shadow3
    #
    #     # beam_shadow3.write(file)
    #     # print("File %s written to disk. "%file)
    #
    # def dump_shadow3_file(self,file):
    #     #TODO this dump uses now shadow3. To be removed after checking or write using fully python
    #     beam3 = self.get_shadow3_beam()
    #     beam3.write(file)
    #     print("File %s written to disk. "%file)


    #
    #  interfaces like in shadow3
    #

    def genSource(self,source_object):
        tmp = source_object.get_beam()
        self.rays = tmp.get_rays()
        return tmp

    def traceOE(self,oe_object,n,overwrite=True):
        beam_result = oe_object.trace_beam(self)
        if overwrite:
            self.rays = beam_result.rays
        return beam_result

    #
    # histograms
    #
    # TODO: histo1, histo2

    @classmethod
    def column_names(cls):
        return [
                        "X spatial coordinate [user's unit]",
                        "Y spatial coordinate [user's unit]",
                        "Z spatial coordinate [user's unit]",
                        "Xp direction or divergence [rads]",
                        "Yp direction or divergence [rads]",
                        "Zp direction or divergence [rads]",
                        "X component of the electromagnetic vector (s-polariz)",
                        "Y component of the electromagnetic vector (s-polariz)",
                        "Z component of the electromagnetic vector (s-polariz)",
                        "Lost ray flag",
                        "Energy [eV]",
                        "Wavenumber [cm-1]",
                        "Optical path length",
                        "Phase (s-polarization) in rad",
                        "Phase (p-polarization) in rad",
                        "X component of the electromagnetic vector (p-polariz)",
                        "Y component of the electromagnetic vector (p-polariz)",
                        "Z component of the electromagnetic vector (p-polariz)",
                        "Wavelength [A]",
                        "R= SQRT(X^2+Y^2+Z^2)",
                        "angle from Y axis",
                        "the magnitude of the Electromagnetic vector",
                        "|E|^2 (total intensity)",
                        "total intensity for s-polarization",
                        "total intensity for p-polarization",
                        "photon energy [eV]",
                        "K = 2 pi / lambda * col4 [A^-1]",
                        "K = 2 pi / lambda * col5 [A^-1]",
                        "K = 2 pi / lambda * col6 [A^-1]",
                        "S0-stokes = |Es|^2 + |Ep|^2",
                        "S1-stokes = |Es|^2 - |Ep|^2",
                        "S2-stokes = 2 |Es| |Ep| cos(phase_s-phase_p)",
                        "S3-stokes = 2 |Es| |Ep| sin(phase_s-phase_p)",
                        "Power = intensity(col 23) * energy (col 11)",
                        "Angle-X with Y: |arcsin(X')|",
                        "Angle-Z with Y: |arcsin(Z')|",
                        "Angle-X with Y: |arcsin(X') - mean(arcsin(X'))|",
                        "Angle-Z with Y: |arcsin(Z') - mean(arcsin(Z'))|",
                ]

    @classmethod
    def column_short_names(cls):
        return [
                        "x","y","z",
                        "v_x", "v_y", "v_z",
                        "Es_x", "Es_y", "Es_z",
                        "flag","K [cm-1]","idx","OpPath",
                        "PhaseS","PhaseP",
                        "Ep_x", "Ep_y", "Ep_z",
                        "lambda [A]",
                        "R",
                        "angle from Y",
                        "|E|",
                        "I tot",
                        "I s-pol",
                        "I p-pol",
                        "E[eV]",
                        "Kx",
                        "Ky",
                        "Kz",
                        "S0",
                        "S1",
                        "S2",
                        "S3",
                        "PWR",
                        "Angle X^Y",
                        "Angle Z^Y",
                        "Angle X^Y",
                        "Angle Z^Y",
                ]

    @classmethod
    def column_short_names_with_column_number(cls):
        names = cls.column_short_names()
        for i in range(len(names)):
            names[i] = "col%02d %s" % (i+1, names[i])
        return names

    def write(self,filename,overwrite=True,simulation_name="run001",beam_name="begin"):

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
    def load(cls,filename,simulation_name="run001",beam_name="begin"):

        import h5py

        f = h5py.File(filename, 'r')

        column_names = cls.column_short_names_with_column_number()

        try:
            x = (f["%s/%s/col01 x"%(simulation_name,beam_name)])[:]

            beam = Beam(N=x.size)
            rays = numpy.zeros( (x.size,18))
            for i in range(18):
                column_name = column_names[i]
                print(">>>>>","%s/%s/%s"%(simulation_name,beam_name,column_name))
                rays[:,i] = (f["%s/%s/%s"%(simulation_name,beam_name,column_name)])[:]
        except:
            f.close()
            raise Exception("Cannot find data in %s:/%s/%s" % (filename, simulation_name, beam_name))

        f.close()

        return Beam.initialize_from_array(rays)

    def identical(self,beam2):
        try:
            assert_almost_equal(self.rays,beam2.rays)
            return True
        except:
            return False

    def difference(self,beam2):
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
    # pass

    print(Beam().column_short_names_with_column_number())