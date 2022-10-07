
import numpy
import sys
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QWidget, QLabel, QSizePolicy
from oasys.widgets import gui as oasysgui
from orangewidget import gui, widget

from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from mpl_toolkits.mplot3d import Axes3D



class ShowSurfaceShapeDialog(QWidget):


    def __init__(self, parent=None, ccc=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                 x_min=-10, x_max=10, y_min=-10, y_max=10, branch=0,
                 title='O.E. Surface Shape'):

        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        self.branch = branch


        QDialog.__init__(self, parent)
        self.setWindowTitle(title)

        self.setFixedWidth(1000)

        layout = QtWidgets.QGridLayout(self)
        self.layout = layout

        figure = Figure(figsize=(100, 100))
        figure.patch.set_facecolor('white')

        axis = figure.add_subplot(111, projection='3d')

        axis.set_xlabel("X [m]")
        axis.set_ylabel("Y [m]")
        axis.set_zlabel("Z [m]")

        figure_canvas = FigureCanvasQTAgg(figure)
        figure_canvas.setFixedWidth(700)
        figure_canvas.setFixedHeight(650)

        if parent is None:

            pass

        else:
            self.x_min = -parent.dim_x_minus
            self.x_max = parent.dim_x_plus
            self.y_min = -parent.dim_y_minus
            self.y_max = parent.dim_y_plus

        x = numpy.linspace(self.x_min, self.x_max, 100)
        y = numpy.linspace(self.y_min, self.y_max, 100)

        X, Y = numpy.meshgrid(x, y)

        if parent is None:
            self.c1 = ccc[0]
            self.c2 = ccc[1]
            self.c3 = ccc[2]
            self.c4 = ccc[3]
            self.c5 = ccc[4]
            self.c6 = ccc[5]
            self.c7 = ccc[6]
            self.c8 = ccc[7]
            self.c9 = ccc[8]
            self.c10= ccc[9]


        else:
            self.c1 = round(parent.conic_coefficient_0, 10)
            self.c2 = round(parent.conic_coefficient_1, 10)
            self.c3 = round(parent.conic_coefficient_2, 10)
            self.c4 = round(parent.conic_coefficient_3, 10)
            self.c5 = round(parent.conic_coefficient_4, 10)
            self.c6 = round(parent.conic_coefficient_5, 10)
            self.c7 = round(parent.conic_coefficient_6, 10)
            self.c8 = round(parent.conic_coefficient_7, 10)
            self.c9 = round(parent.conic_coefficient_8, 10)
            self.c10= round(parent.conic_coefficient_9, 10)

        c = self.c1*(X**2) + self.c2*(Y**2) + self.c4*X*Y + self.c7*X + self.c8*Y + self.c10
        b = self.c5*Y + self.c6*X + self.c9
        a = self.c3

        if self.branch == 0:
            mysign = -1.0
        else:
            mysign= 1.0
        z_values = (-b + mysign * numpy.sqrt(b**2 - 4*a*c))/(2*a)
        z_values[b**2 - 4*a*c < 0] = numpy.nan

        axis.plot_surface(X, Y, z_values,
                          rstride=1, cstride=1, cmap=cm.autumn, linewidth=0.5, antialiased=True)

        title_head = "Surface from generated conic coefficients:\n"
        title = ""
        max_dim = 40

        if self.c1 != 0: title +=       str(self.c1) + u"\u00B7" + "X" + u"\u00B2"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c2 < 0 or (self.c2 > 0 and title == ""): title +=       str(self.c2) + u"\u00B7" + "Y" + u"\u00B2"
        elif self.c2 > 0                                 : title += "+" + str(self.c2) + u"\u00B7" + "Y" + u"\u00B2"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c3 < 0 or (self.c3 > 0 and title == ""): title +=       str(self.c3) + u"\u00B7" + "Z" + u"\u00B2"
        elif self.c3 > 0                                 : title += "+" + str(self.c3) + u"\u00B7" + "Z" + u"\u00B2"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c4 < 0 or (self.c4 > 0 and title == ""): title +=       str(self.c4) + u"\u00B7" + "XY"
        elif self.c4 > 0                                 : title += "+" + str(self.c4) + u"\u00B7" + "XY"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c5 < 0 or (self.c5 > 0 and title == ""): title +=       str(self.c5) + u"\u00B7" + "YZ"
        elif self.c5 > 0                                 : title += "+" + str(self.c5) + u"\u00B7" + "YZ"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c6 < 0 or (self.c6 > 0 and title == ""): title +=       str(self.c6) + u"\u00B7" + "XZ"
        elif self.c6 > 0                                 : title += "+" + str(self.c6) + u"\u00B7" + "XZ"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c7 < 0 or (self.c7 > 0 and title == ""): title +=       str(self.c7) + u"\u00B7" + "X"
        elif self.c7 > 0                                 : title += "+" + str(self.c7) + u"\u00B7" + "X"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c8 < 0 or (self.c8 > 0 and title == ""): title +=       str(self.c8) + u"\u00B7" + "Y"
        elif self.c8 > 0                                 : title += "+" + str(self.c8) + u"\u00B7" + "Y"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c9 < 0 or (self.c9 > 0 and title == ""): title +=       str(self.c9) + u"\u00B7" + "Z"
        elif self.c9 > 0                                 : title += "+" + str(self.c9) + u"\u00B7" + "Z"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c10< 0 or (self.c10> 0 and title == ""): title +=       str(self.c10)
        elif self.c10> 0                                 : title += "+" + str(self.c10)

        axis.set_title(title_head + title + " = 0")

        print(title_head + title + " = 0")

        figure_canvas.draw()

        axis.mouse_init()

        widget = QWidget(parent=self)

        surface_box = oasysgui.widgetBox(widget, "Conic Coefficients", addSpace=False, orientation="vertical", width=240, height=750)

        label  = "c[1]" + u"\u00B7" + "X" + u"\u00B2" + " + c[2]" + u"\u00B7" + "Y" + u"\u00B2" + " + c[3]" + u"\u00B7" + "Z" + u"\u00B2" + " +\n"
        label += "c[4]" + u"\u00B7" + "XY" + " + c[5]" + u"\u00B7" + "YZ" + " + c[6]" + u"\u00B7" + "XZ" + " +\n"
        label += "c[7]" + u"\u00B7" + "X" + " + c[8]" + u"\u00B7" + "Y" + " + c[9]" + u"\u00B7" + "Z" + " + c[10] = 0"

        gui.label(surface_box, self, label)

        gui.separator(surface_box, 10)

        oasysgui.lineEdit(surface_box, self, "c1" , "c[1]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c2" , "c[2]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c3" , "c[3]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c4" , "c[4]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c5" , "c[5]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c6" , "c[6]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c7" , "c[7]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c8" , "c[8]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c9" , "c[9]" , labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "c10", "c[10]", labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)


        oasysgui.lineEdit(surface_box, self, "x_min", "xmin", labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "x_max", "xmax", labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "y_min", "ymin", labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)
        oasysgui.lineEdit(surface_box, self, "y_max", "ymax", labelWidth=60, valueType=float, orientation="horizontal", callback=self.refresh)

        gui.comboBox(surface_box, self, "branch",
                    label="branch", addSpace=False,
                    items=['-1','+1'],
                    valueType=int, orientation="horizontal", callback=self.refresh)


        layout.addWidget(figure_canvas, 0, 0)
        layout.addWidget(widget, 0, 1)


        self.setLayout(layout)

        self.axis = axis
        self.figure_canvas = figure_canvas

    def refresh(self):
        print("refresh")


        x = numpy.linspace(self.x_min, self.x_max, 100)
        y = numpy.linspace(self.y_min, self.y_max, 100)

        X, Y = numpy.meshgrid(x, y)


        c = self.c1*(X**2) + self.c2*(Y**2) + self.c4*X*Y + self.c7*X + self.c8*Y + self.c10
        b = self.c5*Y + self.c6*X + self.c9
        a = self.c3

        if self.branch == 0:
            mysign = -1.0
        else:
            mysign= 1.0
        z_values = (-b + mysign * numpy.sqrt(b**2 - 4*a*c))/(2*a)
        z_values[b**2 - 4*a*c < 0] = numpy.nan

        self.axis.clear()
        self.axis.set_xlabel("X [m]")
        self.axis.set_ylabel("Y [m]")
        self.axis.set_zlabel("Z [m]")
        self.axis.plot_surface(X, Y, z_values,
                          rstride=1, cstride=1, cmap=cm.autumn, linewidth=0.5, antialiased=True)

        title_head = "Surface from generated conic coefficients:\n"
        title = ""
        max_dim = 40

        if self.c1 != 0: title +=       str(self.c1) + u"\u00B7" + "X" + u"\u00B2"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c2 < 0 or (self.c2 > 0 and title == ""): title +=       str(self.c2) + u"\u00B7" + "Y" + u"\u00B2"
        elif self.c2 > 0                                 : title += "+" + str(self.c2) + u"\u00B7" + "Y" + u"\u00B2"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c3 < 0 or (self.c3 > 0 and title == ""): title +=       str(self.c3) + u"\u00B7" + "Z" + u"\u00B2"
        elif self.c3 > 0                                 : title += "+" + str(self.c3) + u"\u00B7" + "Z" + u"\u00B2"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c4 < 0 or (self.c4 > 0 and title == ""): title +=       str(self.c4) + u"\u00B7" + "XY"
        elif self.c4 > 0                                 : title += "+" + str(self.c4) + u"\u00B7" + "XY"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c5 < 0 or (self.c5 > 0 and title == ""): title +=       str(self.c5) + u"\u00B7" + "YZ"
        elif self.c5 > 0                                 : title += "+" + str(self.c5) + u"\u00B7" + "YZ"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c6 < 0 or (self.c6 > 0 and title == ""): title +=       str(self.c6) + u"\u00B7" + "XZ"
        elif self.c6 > 0                                 : title += "+" + str(self.c6) + u"\u00B7" + "XZ"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c7 < 0 or (self.c7 > 0 and title == ""): title +=       str(self.c7) + u"\u00B7" + "X"
        elif self.c7 > 0                                 : title += "+" + str(self.c7) + u"\u00B7" + "X"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c8 < 0 or (self.c8 > 0 and title == ""): title +=       str(self.c8) + u"\u00B7" + "Y"
        elif self.c8 > 0                                 : title += "+" + str(self.c8) + u"\u00B7" + "Y"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c9 < 0 or (self.c9 > 0 and title == ""): title +=       str(self.c9) + u"\u00B7" + "Z"
        elif self.c9 > 0                                 : title += "+" + str(self.c9) + u"\u00B7" + "Z"
        if len(title) >=  max_dim:
            title_head += title + "\n"
            title = ""
        if self.c10< 0 or (self.c10> 0 and title == ""): title +=       str(self.c10)
        elif self.c10> 0                                 : title += "+" + str(self.c10)

        self.axis.set_title(title_head + title + " = 0")

        print(title_head + title + " = 0")

        self.figure_canvas.draw()

        self.axis.mouse_init()

def view_conic(ccc=[1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,-764,0.0], x_min=-10, x_max=10, y_min=-10, y_max=10, branch=0):
    app = QApplication(sys.argv)
    dialog = ShowSurfaceShapeDialog(ccc=ccc, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, branch=0)
    dialog.show()
    app.exec()

def compare_conics(ccc1=[1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,-764,0.0],
                   ccc2=[1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,764,0.0],
                   x_min=-10, x_max=10, y_min=-10, y_max=10,
                   titles=['',''], branch=0):
    app = QApplication(sys.argv)
    dialog = ShowSurfaceShapeDialog(ccc=ccc1, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, title=titles[0], branch=branch)
    dialog.show()
    dialog2 = ShowSurfaceShapeDialog(ccc=ccc2, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, title=titles[1], branch=branch)
    dialog2.show()


    app.exec()

if __name__ == "__main__":
    # view_conic(ccc=[1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,-764,0.0], x_min=-10, x_max=10, y_min=-10, y_max=10)
    compare_conics(titles=["A","B"])
