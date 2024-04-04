import socket
import getpass
import datetime

from shadow4.beam.s4_beam import S4Beam

try:
    import matplotlib.pyplot as plt
    import matplotlib
except:
    print("Please install matplotlib to allow graphics")

import numpy
import sys
import os

import matplotlib.pylab as plt
from matplotlib import collections

codata_h = numpy.array(6.62606957e-34)
codata_ec = numpy.array(1.602176565e-19)
codata_c = numpy.array(299792458.0)
A2EV = 2.0*numpy.pi/(codata_h*codata_c/codata_ec*1e2)


class ArgsError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NoValueSelectedError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Histo1_Ticket:
    def __init__(self):
        self.histogram = None
        self.bin_center = None
        self.bin_left = None
        self.figure = None
        self.xrange = None
        self.yrange = None
        self.xtitle = None
        self.ytitle = None
        self.title = None
        self.fwhm = None


class plotxy_Ticket:
    def __init__(self):
        self.figure = None
        self.xrange = None
        self.yrange = None
        self.xtitle = None
        self.ytitle = None
        self.title = None
        self.fwhmx = None
        self.fwhmy = None


def ErrorMsg(fromFunc, value):
    print(fromFunc + " called with an error in the arguments.\nCheck help function")
    return ArgsError(value)


def getshonecol_CheckArg(beam, col):
    if not isinstance(beam, S4Beam): raise ErrorMsg('getshonecol', 'beam')
    if not isinstance(col, int): raise ErrorMsg('getshonecol', 'col')
    if col < 1 or col > 33: raise ErrorMsg('getshonecol', 'col')


def getshcol_CheckArg(beam, col):  # the rest of checks are included in the function getshonecol_CheckArg
    if not isinstance(beam, S4Beam): raise ErrorMsg('getshcol', 'beam')
    if not isinstance(col, (int, tuple, list)): raise ErrorMsg('getshcol', 'col')
    if isinstance(col, int):
        if col < 1 or col > 33: raise ErrorMsg('getshcol', 'col')
    else:
        for c in col:
            if not isinstance(c, int): raise ErrorMsg('getshcol', 'col')
            if c < 1 or c > 33: raise ErrorMsg('getshcol', 'col')


def Histo1_CheckArg(beam, col, xrange, yrange, nbins, nolost, ref, write, title, xtitle, ytitle, calfwhm, noplot):
    if not isinstance(beam, S4Beam): raise ErrorMsg('Histo1', 'beam')
    if not isinstance(col, int): raise ErrorMsg('Histo1', 'col')
    if col < 1 or col > 33: raise ErrorMsg('Histo1', 'col')
    # the next 3 lines don't matter, it is a trick to pass the test when None
    if xrange == None: xrange = (1.0, 2.0)
    if yrange == None: yrange = (1.0, 2.0)
    if xtitle == None: xtitle = 'pippo'
    if ytitle == None: ytitle = 'pippo'
    if not isinstance(xrange, (tuple, list)): raise ErrorMsg('Histo1', 'xrange')
    if len(xrange) != 2: raise ErrorMsg('Histo1', 'xrange')
    if not isinstance(xrange[0], (int, float)) or not isinstance(xrange[1], (int, float)): raise ErrorMsg('Histo1',
                                                                                                          'xrange')
    if not isinstance(yrange, (tuple, list)): raise ErrorMsg('Histo1', 'yrange')
    if len(yrange) != 2: raise ErrorMsg('Histo1', 'yrange')
    if not isinstance(yrange[0], (int, float)) or not isinstance(yrange[1], (int, float)): raise ErrorMsg('Histo1',
                                                                                                          'yrange')
    if not isinstance(nbins, int): raise ErrorMsg('Histo1', 'nbins')
    if nbins <= 0: raise ErrorMsg('Histo1', 'nbins')
    if nolost != 0 and nolost != 1 and nolost != 2: raise ErrorMsg('Histo1', 'nolost')
    if ref >= 22 and ref <= 33: ref = 1
    if ref != 0 and ref != 1: raise ErrorMsg('Histo1', 'ref')
    if write != 0 and write != 1: raise ErrorMsg('Histo1', 'write')
    if not isinstance(title, str): raise ErrorMsg('Histo1', 'title')
    if not isinstance(xtitle, str): raise ErrorMsg('Histo1', 'xtitle')
    if not isinstance(ytitle, str): raise ErrorMsg('Histo1', 'ytitle')
    if calfwhm != 0 and calfwhm != 1: raise ErrorMsg('Histo1', 'calfwhm')
    if noplot != 0 and noplot != 1: raise ErrorMsg('Histo1', 'noplot')


def plotxy_CheckArg(beam, cols1, cols2, nbins, nbins_h, level, xrange, yrange, nolost, title, xtitle, ytitle, noplot,
                    calfwhm, contour):
    if not isinstance(beam, S4Beam): raise ErrorMsg('plotxy', 'beam')
    if cols1 < 1 or cols1 > 33: raise ErrorMsg('plotxy', 'cols1')
    if cols2 < 1 or cols2 > 33: raise ErrorMsg('plotxy', 'cols2')
    if not isinstance(nbins, int): raise ErrorMsg('plotxy', 'nbins')
    if nbins <= 0: raise ErrorMsg('plotxy', 'nbins')
    if not isinstance(nbins_h, int): raise ErrorMsg('plotxy', 'nbins_h')
    if nbins_h <= 0: raise ErrorMsg('plotxy', 'nbins_h')
    if not isinstance(level, int): raise ErrorMsg('plotxy', 'level')
    if level <= 0: raise ErrorMsg('plotxy', 'level')
    # the next 4 lines don't matter, it is a trick to pass the test when None
    if xrange == None: xrange = (1.0, 2.0)
    if yrange == None: yrange = (1.0, 2.0)
    if xtitle == None: xtitle = 'pippo'
    if ytitle == None: ytitle = 'pippo'
    if not isinstance(xrange, (tuple, list)): raise ErrorMsg('plotxy', 'xrange')
    if len(xrange) != 2: raise ErrorMsg('plotxy', 'xrange')
    if not isinstance(xrange[0], (int, float)) or not isinstance(xrange[1], (int, float)): raise ErrorMsg('plotxy',
                                                                                                          'xrange')
    if not isinstance(yrange, (tuple, list)): raise ErrorMsg('plotxy', 'yrange')
    if len(yrange) != 2: raise ErrorMsg('plotxy', 'yrange')
    if not isinstance(yrange[0], (int, float)) or not isinstance(yrange[1], (int, float)): raise ErrorMsg('plotxy',
                                                                                                          'yrange')
    if nolost != 0 and nolost != 1 and nolost != 2: raise ErrorMsg('plotxy', 'nolost')
    if not isinstance(title, str): raise ErrorMsg('plotxy', 'title')
    if not isinstance(xtitle, str): raise ErrorMsg('plotxy', 'xtitle')
    if not isinstance(ytitle, str): raise ErrorMsg('plotxy', 'ytitle')
    if noplot != 0 and noplot != 1: raise ErrorMsg('plotxy', 'noplot')
    # if ref!=0 and ref!=1: raise ErrorMsg('plotxy','ref')
    if calfwhm != 0 and calfwhm != 1 and calfwhm != 2: raise ErrorMsg('plotxy', 'calfwhm')
    if not isinstance(contour, int): raise ErrorMsg('plotxy', 'contour')
    if contour < 0 or contour > 6: raise ErrorMsg('plotxy', 'contour')


def setGoodRange(col):
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


def findIndex(xx, n, la, lb):
    return int(numpy.floor((xx - (lb - la) * 0.5 / n - la) * n / (lb - la)))


def calcFWHM(h, binSize):
    t = numpy.where(h > max(h) * 0.5)
    return binSize * (t[0][-1] - t[0][0] + 1), t[0][-1], t[0][0]


def Histo1_write(title, bins, h, w, col, beam, ref):
    if isinstance(beam, S4Beam): usubtitle = "Shadow running in dir " + os.getcwd()
    if isinstance(beam, str): usubtitle = os.getcwd() + beam
    now = str(datetime.datetime.now())
    usubtitle += " " + now + " " + getpass.getuser() + "@" + socket.gethostname()
    file = open(title, 'w')
    print("#F HISTO1", file=file)
    print("#C This file has been created using histo1 (python ShadowTools)", file=file)
    print("#D " + now, file=file)
    print("#UTITLE", file=file)
    print("#USUBTITLE " + usubtitle, file=file)
    print("#UTTEXT", file=file)
    print("#C COLUMN 1 CORRESPONDS TO ABSCISSAS IN THE CENTER OF EACH BIN", file=file)
    print("#C COLUMN 2 CORRESPONDS TO ABSCISSAS IN THE THE LEFT CORNER OF THE BIN", file=file)
    print("#C COLUMN 3 CORRESPONDS TO INTENSITY (COUNTING RAYS)", file=file)
    print("#C COLUMN 4 CORRESPONDS TO INTENSITY (WEIGHTED IF SELECTED)", file=file)
    print(" ", file=file)
    print("#S 1 histogram", file=file)
    print("#N 4", file=file)
    print("#L " + getLabel(col)[1] + "  " + (getLabel(col))[1] + "  " + "intensity (rays)" + "  " + (getLabel(ref))[1],
          file=file)
    for i in range(len(h)):
        print("%f\t%f\t%f\t%f" % ((bins[i] + bins[i + 1]) * 0.5, bins[i], h[i], w[i]), file=file)
    file.close()


def getLabel(col):
    Label = [
        [r"$x$ [user unit]", "x [user unit]"],
        [r"$y$ [user unit]", "y [user unit]"],
        [r"$z$ [user unit]", "z [user unit]"],
        [r"$\dot{x}$ [rads]", "x' [rads]"],
        [r"$\dot{y}$ [rads]", "y' [rads]"],
        [r"$\dot{z}$ [rads]", "z' [rads]"],
        [r"$\mathbf{E}_{\sigma x}$", "Es_x"],
        [r"$\mathbf{E}_{\sigma y}$", "Es_y"],
        [r"$\mathbf{E}_{\sigma z}$", "Es_z"],
        [r"ray flag", "Ray flag"],
        [r"E [eV]", "Energy"],
        [r"Ray index", "Ray index"],
        [r"s", "Opt. Path"],
        [r"$\phi_{\sigma}$", "phase_s"],
        [r"$\phi_{\pi}$", "phase_p"],
        [r"$\mathbf{E}_{\pi x}$", "Ep_x"],
        [r"$\mathbf{E}_{\pi y}$", "Ep_y"],
        [r"$\mathbf{E}_{\pi z}$", "Ep_z"],
        [r"$\lambda$ [$\AA$]", "wavelength"],
        [r"$R= \sqrt{x^2+y^2+z^2}$", "R [user unit]"],
        [r"$\theta$", "theta"],
        [r"$\Vert\mathbf{E_{\sigma}}+\mathbf{E_{\pi}}\Vert$", "Electromagnetic vector magnitude"],
        [r"$\Vert\mathbf{E_{\sigma}}+\mathbf{E_{\pi}}\Vert^2$",
         "intensity (weight column = 23: |E|^2 (total intensity))"],
        [r"$\Vert\mathbf{E_{\sigma}}\Vert^2$", "intensity (weight column = 24: |E_s|^2 (sigma intensity))"],
        [r"$\Vert\mathbf{E_{\pi}}\Vert^2$", "intensity (weight column = 25: |E_p|^2 (pi intensity))"],
        [r"$K = \frac{2 \pi}{\lambda} [A^{-1}]$", "K magnitude"],
        [r"$K_x = \frac{2 \pi}{\lambda} \dot{x}$ [$\AA^{-1}$]", "K_x"],
        [r"$K_y = \frac{2 \pi}{\lambda} \dot{y}$ [$\AA^{-1}$]", "K_y"],
        [r"$K_z = \frac{2 \pi}{\lambda} \dot{z}$ [$\AA^{-1}$]", "K_z"],
        [r"$S_0 = \Vert\mathbf{E}_{\sigma}\Vert^2 + \Vert\mathbf{E}_{\pi}\Vert^2 $", "S0"],
        [r"$S_1 = \Vert\mathbf{E}_{\sigma}\Vert^2 - \Vert\mathbf{E}_{\pi}\Vert^2 $", "S1"],
        [r"$S_2 = 2 \Vert\mathbf{E}_{\sigma}\Vert \cdot \Vert\mathbf{E}_{\pi}\Vert \cos{(\phi_{\sigma}-\phi_{\pi})}$",
         "S2"],
        [r"$S_3 = 2 \Vert\mathbf{E}_{\sigma}\Vert \cdot \Vert\mathbf{E}_{\pi}\Vert \sin{(\phi_{\sigma}-\phi_{\pi})}$",
         "S3"],
        [r"Power [eV/s]", "Power"],
    ]
    return Label[col]

def histo1(beam, col, notitle=0, nofwhm=0,  bar=0,  **kwargs):
    """
    Plot the histogram of a column, as calculated by Shadow.Beam.histo1 using matplotlib

    NOTE: This will replaces the old histo1 still available as histo1_old

    :param beam: a Shadow.Beam() instance, or a file name with Shadow binary file
    :param col: the Shadow column number (start from 1)
    :param notitle: set to 1 to avoid displaying title
    :param nofwhm: set to 1 to avoid labeling FWHM value
    :param bar: 1=bar plot, 0=line plot
    :param kwargs: keywords accepted by Shadow.Beam.histo1()
    :return: the dictionary returned by Shadow.beam.histo1() with some keys added.
    """

    title = "histo1"

    if isinstance(beam,str):
        raise NotImplementedError()

    tk2 = beam.histo1(col, **kwargs)



    h = tk2["histogram"]
    bins = tk2["bin_left"]
    xrange = tk2["xrange"]
    yrange = [0,1.1*numpy.max(h)]
    fwhm = tk2["fwhm"]

    xtitle = "column %d"%tk2["col"]
    ytitle = "counts ("

    if tk2["nolost"] == 0:
        ytitle += " all rays"
    if tk2["nolost"] == 1:
        ytitle += " good rays"
    if tk2["nolost"] == 2:
        ytitle += " lost rays"

    if tk2["ref"] == 0:
        ytitle += " = weight: number of rays"
    else:
        if tk2["ref"] == 23:
            ytitle += " - weight: intensity"
        else:
            ytitle += " - weight column: %d"%(tk2["ref"])

    ytitle += ")"


    if fwhm != None: print ("fwhm = %g" % fwhm)

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)

    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    if notitle != 1: ax.set_title(title)
    ax.set_xlim(xrange[0],xrange[1])
    ax.set_ylim(yrange[0],yrange[1])
    ax.grid(True)

    if bar:
        l = ax.bar(bins, h, 1.0*(bins[1]-bins[0]),color='blue') #,error_kw=dict(elinewidth=2,ecolor='red'))
    else:
        l = plt.plot(tk2["bin_path"], tk2["histogram_path"], color='blue') #,error_kw=dict(elinewidth=2,ecolor='red'))

    if tk2["fwhm"] != None:
        hh = 0.5*numpy.max(tk2["histogram"])
        lines = [ [ (tk2["fwhm_coordinates"][0],hh), \
                    (tk2["fwhm_coordinates"][1],hh) ]]
        lc = collections.LineCollection(lines,color='red',linewidths=2)
        ax.add_collection(lc)
        if nofwhm != 1:
            if tk2["fwhm_coordinates"][0] < 0:
                shift1 = 0.9
            else:
                shift1 = 1.0
            ax.annotate('FWHM=%f'%tk2["fwhm"], xy=(shift1*tk2["fwhm_coordinates"][0],1.01*tk2["fwhm_coordinates"][0]))

    plt.show()
    return tk2


def plotxy(beam,col_h,col_v, nofwhm=1, title="", **kwargs):
  """

  plotxy implementation using matplotlib.
  Calculations are done using Shadow.beam.histo2()

  :param beam: it can be a SHADOW binary file, an instance of Shadow.Beam() or a dictionary from Shadow.Beam.histo2
  :param col_h: The column for the H coordinate in the plot (irrelevant of beam is a dictionary)
  :param col_v: The column for the H coordinate in the plot (irrelevant of beam is a dictionary)
  :param nofwhm: set to 0 to label the FWHM value in the plot (default do not label)
  :param kwargs: keywrods passed to Shadow.Beam.histo2
  :return: the dictionary returned by Shadow.beam.histo2() with some added keys.
  """
  if title == "":
      title = "plotxy"

  if isinstance(beam,dict):
    tkt = beam
    col_h = tkt["col_h"]
    col_v = tkt["col_v"]
  else:
    if isinstance(beam,str):
      raise NotImplementedError()
    tkt = beam.histo2(col_h,col_v,**kwargs)


  xtitle = "Column %d"%tkt["col_h"]
  ytitle = "Column %d"%tkt["col_v"]

  figure = plt.figure(figsize=(12,8),dpi=96)

  ratio = 8.0/12.0

  rect_scatter = [0.10*ratio, 0.10, 0.65*ratio, 0.65]
  rect_histx =   [0.10*ratio, 0.77, 0.65*ratio, 0.20]
  rect_histy =   [0.77*ratio, 0.10, 0.20*ratio, 0.65]
  rect_text =    [1.00*ratio, 0.10, 1.20*ratio, 0.65]

  #
  #main plot
  #
  axScatter = figure.add_axes(rect_scatter)
  axScatter.set_xlabel(xtitle)
  axScatter.set_ylabel(ytitle)

  # axScatter.set_xlim(tkt["xrange"])
  # axScatter.set_ylim(tkt["yrange"])

  axScatter.axis(xmin=tkt["xrange"][0],xmax=tkt["xrange"][1])
  axScatter.axis(ymin=tkt["yrange"][0],ymax=tkt["yrange"][1])
  #axScatter.pcolor(tkt["bin_h_edges"], tkt["bin_v_edges"], tkt["histogram"].T)
  axScatter.pcolormesh(tkt["bin_h_edges"], tkt["bin_v_edges"], tkt["histogram"].T)

  for tt in axScatter.get_xticklabels():
    tt.set_size('x-small')
  for tt in axScatter.get_yticklabels():
    tt.set_size('x-small')

  #
  #histograms
  #
  axHistx = figure.add_axes(rect_histx, sharex=axScatter)
  axHisty = figure.add_axes(rect_histy, sharey=axScatter)

  #for practical purposes, writes the full histogram path
  tmp_h_b = []
  tmp_h_h = []
  for s,t,v in zip(tkt["bin_h_left"],tkt["bin_h_right"],tkt["histogram_h"]):
    tmp_h_b.append(s)
    tmp_h_h.append(v)
    tmp_h_b.append(t)
    tmp_h_h.append(v)
  tmp_v_b = []
  tmp_v_h = []
  for s,t,v in zip(tkt["bin_v_left"],tkt["bin_v_right"],tkt["histogram_v"]):
    tmp_v_b.append(s)
    tmp_v_h.append(v)
    tmp_v_b.append(t)
    tmp_v_h.append(v)

  axHistx.plot(tmp_h_b,tmp_h_h)
  axHisty.plot(tmp_v_h,tmp_v_b)

  for tl in axHistx.get_xticklabels(): tl.set_visible(False)
  for tl in axHisty.get_yticklabels(): tl.set_visible(False)
  for tt in axHisty.get_xticklabels():
    tt.set_rotation(270)
    tt.set_size('x-small')
  for tt in axHistx.get_yticklabels():
    tt.set_size('x-small')

  if tkt["fwhm_h"] != None:
      hh = 0.5*numpy.max(tkt["histogram_h"])
      lines = [ [ (tkt["fwhm_coordinates_h"][0],hh), \
                  (tkt["fwhm_coordinates_h"][1],hh) ]]
      lc = collections.LineCollection(lines,color='red',linewidths=2)
      axHistx.add_collection(lc)
      if nofwhm != 1:
          if tkt["fwhm_coordinates_h"][0] < 0:
              shift1 = 0.9
          else:
              shift1 = 1.0
          axHistx.annotate('FWHM=%f'%tkt["fwhm_h"], xy=(shift1*tkt["fwhm_coordinates_h"][0],1.01*hh))

  if tkt["fwhm_v"] != None:
      hh = 0.5*numpy.max(tkt["histogram_v"])
      lines = [ [ (hh,tkt["fwhm_coordinates_v"][0]), \
                  (hh,tkt["fwhm_coordinates_v"][1]) ]]
      lc = collections.LineCollection(lines,color='green',linewidths=2)
      axHisty.add_collection(lc)
      if nofwhm != 1:
          if tkt["fwhm_coordinates_v"][0] < 0:
              shift1 = 0.9
          else:
              shift1 = 1.0
          axHisty.annotate('FWHM=%f'%tkt["fwhm_v"], xy=(shift1*tkt["fwhm_coordinates_v"][0],1.01*hh))


  if title!=None:
    axHistx.set_title(title)
  axText = figure.add_axes(rect_text)
  if tkt["nolost"] == 0: axText.text(0.0,0.8,"ALL RAYS")
  if tkt["nolost"] == 1: axText.text(0.0,0.8,"GOOD RAYS")
  if tkt["nolost"] == 2: axText.text(0.0,0.8,"LOST RAYS")


  axText.text(0.0,0.7,"intensity: %8.2f"%(tkt["intensity"]))
  axText.text(0.0,0.6,"total number of rays: "+str(tkt["nrays"]))
  axText.text(0.0,0.5,"total good rays: "+str(tkt["good_rays"]))
  axText.text(0.0,0.4,"total lost rays: "+str(tkt["nrays"]-tkt["good_rays"]))
  calfwhm = 1
  if tkt["fwhm_h"] != None:
    axText.text(0.0,0.3,"fwhm H: "+str(tkt["fwhm_h"]))
  if tkt["fwhm_v"] != None:
    axText.text(0.0,0.2,"fwhm V: "+str(tkt["fwhm_v"]))

  axText.text(0.0,0.1,"from Shadow.Beam instance")
  if tkt["ref"] == 0:
    axText.text(0.0,0.0,"WEIGHT: RAYS")
  else:
      axText.text(0.0,0.0,"WEIGHT: INTENSITY")
  axText.set_axis_off()
  plt.show()
  return tkt
