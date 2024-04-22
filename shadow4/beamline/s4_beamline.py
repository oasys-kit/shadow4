"""
Defines the S4 beamline.

As a reminder, and following the SYNED philosophy, the S4 beamline is a container of the S4 light source and a
list of S4 beamline elements.

"""
from syned.beamline.beamline import Beamline
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
import numpy

from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4Mirror
from shadow4.beamline.optical_elements.crystals.s4_crystal import S4Crystal
from shadow4.beamline.optical_elements.gratings.s4_grating import S4Grating
from shadow4.beamline.optical_elements.refractors.s4_lens import S4Lens
from shadow4.beamline.optical_elements.refractors.s4_transfocator import S4Transfocator
from shadow4.beamline.optical_elements.refractors.s4_crl import S4CRL
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen




class S4Beamline(Beamline):
    """
    Constructor.

    Parameters
    ----------
    light_source : instance of LightSource
        The light source
    beamline_elements_list : list
        The beamline elements (each one an instance of S4BeamlineElement).
    """
    def __init__(self,
                 light_source=None,
                 beamline_elements_list=None):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)

    def duplicate(self):
        """
        Returns a copy of the S4 beamline instance.

        Returns
        -------
        S4BeamlineElement  instance
            A copy of the object instance.

        """
        if self.get_beamline_elements_number() == 0:
            beamline_elements_list = None
        else:
            beamline_elements_list = []
            for beamline_element in self._beamline_elements_list:
                beamline_elements_list.append(beamline_element.duplicate())

        return S4Beamline(light_source=self._light_source,
                          beamline_elements_list = beamline_elements_list)

    def append_beamline_element(self, beamline_element: S4BeamlineElement):
        """
        Appends a S4 beamline element.

        Parameters
        ----------
        beamline_element : instance of S4BeamlineElement.
            The beamline element to append.
        """
        self._beamline_elements_list.append(beamline_element)

    def to_python_code(self, **kwargs):
        """
        Returns the python code to create the beamline.

        Parameters
        ----------
        **kwargs
            Passed arguments

        Returns
        -------
        str
            The python code.
        """
        script = "from shadow4.beamline.s4_beamline import S4Beamline"
        script += "\n\nbeamline = S4Beamline()\n"
        try:
            script += self.get_light_source().to_python_code()
            script += "\n\nbeamline.set_light_source(light_source)"
        except:
            script +=  "\n\n\n# Error getting python code for S4Beamline S4LightSource "

        for i,element in enumerate(self.get_beamline_elements()):
            try:
                script += element.to_python_code()
                script += "\n\nbeamline.append_beamline_element(beamline_element)"
            except:
                script += "\n\n\n# Error getting python code for S4Beamline S4BeamlineElement # %d  :" % (i+1)
                script += "\n#       %s " % (str(element))

        return script

    def run_beamline(self, **params):
        """
        Runs (performs the ray tracing) of the full beamline.

        Parameters
        ----------
        **params
            Passed params.

        Returns
        -------
        tuple
            (output_beam, output_mirr) the traced beam and footprint (after the last beamline element).
        """
        try:
            output_beam = self.get_light_source().get_beam(**params)
            output_mirr = None
        except:
            raise Exception("Error running beamline light source")

        for i, element in enumerate(self.get_beamline_elements()):
            try:
                element.set_input_beam(output_beam)
                output_beam, output_mirr = element.trace_beam(**params)
            except:
                raise Exception("Error running beamline element # %d" % (i+1) )

        return output_beam, output_mirr

    def _get_info_coordinates(self, oe_index):
        coordinates = self.get_beamline_element_at(oe_index).get_coordinates()
        T_SOURCE, T_IMAGE, T_INCIDENCE, T_REFLECTION, ALPHA = coordinates.get_positions()
        print(coordinates)
        txt_coordinates = ""

        txt_coordinates += "Central Axis parameters :                          \n"
        txt_coordinates += "Source Plane Distance                    %f m\n" % T_SOURCE
        txt_coordinates += "Image  Plane                             %f m\n" % T_IMAGE
        txt_coordinates += "Incidence Angle (to normal)              %f deg\n" % numpy.degrees(T_INCIDENCE)
        txt_coordinates += "Reflection/Diffraction Angle (to normal) %f deg\n" % numpy.degrees(T_REFLECTION)
        txt_coordinates += "Grazing Incidence Angle                  %f mrad\n" % (1e3 * (numpy.pi / 2 - T_INCIDENCE))
        txt_coordinates += "Grazing Reflection/Diffraction Angle     %f mrad\n" % (1e3 * (numpy.pi / 2 - T_REFLECTION))

        return txt_coordinates


    def sourcinfo(self):
        """
        Returns the source information (sourcinfo in shadow3).

        Returns
        -------
        str
        """
        from shadow4.sources.source_geometrical.source_grid_cartesian import SourceGridCartesian
        from shadow4.sources.source_geometrical.source_grid_polar import SourceGridPolar

        light_source = self.get_light_source()



        txt_source = \
"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
**************  S O U R C E       D E S C R I P T I O N  **************
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""

        if isinstance(light_source, SourceGridCartesian):
            txt_source += "Grid source (cartesian)\n"
        elif isinstance(light_source, SourceGridPolar):
            txt_source += "Grid source (polar)\n"
        else:
            txt_source += "Random source\n"

        txt_source += "Generated total %d rays.\n" % light_source.get_nrays()

        txt_end = \
"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
***************                 E N D                  ***************
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
        return txt_source + "\n" + self.get_light_source().get_info() + txt_end

    def oeinfo(self, oe_index=None):
        """
        Returns the optical element(s) information (oeinfo or mirinfo in shadow3).

        Returns
        -------
        str
        """
        top_txt1 = \
"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
********************   OPTICAL ELEMENT  DESCRIPTION   ********************
"""
        top_txt2 = \
"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
        bottom_txt = \
"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
***************                 E N D                  ***************
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
        if oe_index is None:
            txt = ""
            for i, element in enumerate(self.get_beamline_elements()):
                txt += top_txt1
                txt += "O.E. %d" % (i+1)
                txt += top_txt2
                txt += self._get_info_coordinates(i)
                txt += element.get_optical_element().get_info()
                txt += bottom_txt
        else:
            txt = ""
            txt += top_txt1
            txt += "O.E. %d" % (oe_index + 1)
            txt += top_txt2
            txt += self._get_info_coordinates(i)
            txt += self.get_beamline_element_at(oe_index).get_optical_element().get_info()
            txt += bottom_txt

        return txt


    def syspositions(self):
        """
        Returns a dictionary with the positions of the source, o.e. and images.

        Returns
        -------
        dict
            dict["source"]  numpy  array (icoor) with the three coordinates of source
            dict["source"]  numpy  array (n_oe,icoor) with the three coordinates of optical elements for all elements
            dict["source"]  numpy  array (n_oe,icoor) with the three coordinates of image positions for all elements
            dict["optical_axis_x"] numpy array with the x coordinates of the optical axis
            dict["optical_axis_y"] numpy array with the y coordinates of the optical axis
            dict["optical_axis_z"] numpy array with the z coordinates of the optical axis
        """
        #
        # This is a translation of the OPTAXIS routine in shadow3
        #
        CENTRAL = numpy.zeros( (25) )
        U_VEC =   numpy.zeros( (4) )
        V_VEC =   numpy.zeros( (4) )
        W_VEC =   numpy.zeros( (4) )
        V_REF =   numpy.zeros( (4) )
        V_PERP =  numpy.zeros( (4) )

        n_elements = self.get_beamline_elements_number()

        mirr = numpy.zeros((3, n_elements))
        star = numpy.zeros_like(mirr)


        for i in range(n_elements):

            oe = self.get_beamline_element_at(i)
            coor = oe.get_coordinates()

            T_SOURCE, T_IMAGE, T_INCIDENCE, T_REFLECTION, ALPHA = coor.get_positions()
            T_IMAGE += oe.get_optical_element().interthickness()

            COSAL = numpy.cos(ALPHA)
            SINAL = numpy.sin(ALPHA)
            PIHALF = numpy.pi/2
            PI = numpy.pi
            DEFLECTION	=   T_INCIDENCE + T_REFLECTION
            if i == 0:
                U_VEC [1]	=   COSAL
                U_VEC [2]	=   0.0
                U_VEC [3]   =   SINAL
                V_VEC [1]   = - numpy.sin(PIHALF - T_INCIDENCE)*SINAL
                V_VEC [2]   =   numpy.cos(PIHALF - T_INCIDENCE)
                V_VEC [3]   =   numpy.sin(PIHALF - T_INCIDENCE)*COSAL
                W_VEC [1]   = - numpy.sin(PI - T_INCIDENCE)*SINAL
                W_VEC [2]   =   numpy.cos(PI - T_INCIDENCE)
                W_VEC [3]   =   numpy.sin(PI - T_INCIDENCE)*COSAL
                V_REF [1]   = - numpy.sin(PI - DEFLECTION)*SINAL
                V_REF [2]   =   numpy.cos(PI - DEFLECTION)
                V_REF [3]   =   numpy.sin(PI - DEFLECTION)*COSAL
                V_PERP [1]  = - numpy.sin(3*PIHALF - DEFLECTION)*SINAL
                V_PERP [2]  =   numpy.cos(3*PIHALF - DEFLECTION)
                V_PERP [3]  =   numpy.sin(3*PIHALF - DEFLECTION)*COSAL

                CENTRAL[1]  =   .0
                CENTRAL[2]  =   .0
                CENTRAL[3]  =   .0
                CENTRAL[4]  =   .0
                CENTRAL[5]  =   T_SOURCE
                CENTRAL[6]  =   .0
                CENTRAL[7] 	=   T_IMAGE*V_REF[1]
                CENTRAL[8] 	=   T_IMAGE*V_REF[2] + T_SOURCE
                CENTRAL[9] 	=   T_IMAGE*V_REF[3]
                CENTRAL[10]	=   U_VEC[1]
                CENTRAL[11]	=   U_VEC[2]
                CENTRAL[12]	=   U_VEC[3]
                CENTRAL[13]	=   V_VEC[1]
                CENTRAL[14]	=   V_VEC[2]
                CENTRAL[15]	=   V_VEC[3]
                CENTRAL[16]	=   W_VEC[1]
                CENTRAL[17]	=   W_VEC[2]
                CENTRAL[18]	=   W_VEC[3]
                CENTRAL[19]	=   V_REF[1]
                CENTRAL[20]	=   V_REF[2]
                CENTRAL[21]	=   V_REF[3]
                CENTRAL[22] =   V_PERP[1]
                CENTRAL[23] =   V_PERP[2]
                CENTRAL[24] =   V_PERP[3]

                source = numpy.zeros(3)
                source [0] = CENTRAL[1]
                source [1] = CENTRAL[2]
                source [2] = CENTRAL[3]

            else:
                #    ! ** Computes now the OLD mirror reference frame in the lab. coordinates
                #    ! ** system. The rotation angle ALPHA of the current mirror is defined in
                #    ! ** this reference frame, as ALPHA measure the angle between the two
                #    ! ** incidence planes (not necessarily the same).
                CENTRAL_OLD = CENTRAL.copy()
                U_OLD = numpy.array(  [0.0,CENTRAL[10],CENTRAL[11],CENTRAL[12]])
                R_OLD = numpy.array(  [0.0,CENTRAL[19],CENTRAL[20],CENTRAL[21]])
                RP_OLD = numpy.array( [0.0,CENTRAL[22],CENTRAL[23],CENTRAL[24]])

                #    ! ** This vector is the NORMAL of the new mirror in the OMRF (U,R_OLD,RP_OLD) **
                V_TEMP = numpy.zeros(4)
                V_TEMP [1]	= - numpy.sin(PI - T_INCIDENCE)*SINAL
                V_TEMP [2]	=   numpy.cos(PI - T_INCIDENCE)
                V_TEMP [3]	=   numpy.sin(PI - T_INCIDENCE)*COSAL

                #    ! ** Rotate it finally to (x,y,z) SRF **
                W_VEC [1]	=    V_TEMP[1]*U_OLD[1] + V_TEMP[2]*R_OLD[1] + V_TEMP[3]*RP_OLD[1]
                W_VEC [2]	=    V_TEMP[1]*U_OLD[2] + V_TEMP[2]*R_OLD[2] + V_TEMP[3]*RP_OLD[2]
                W_VEC [3]	=    V_TEMP[1]*U_OLD[3] + V_TEMP[2]*R_OLD[3] + V_TEMP[3]*RP_OLD[3]

                #    ! ** This vector is the reflected beam from the new mirror in the OMRF **
                V_TEMP[1] = -  numpy.sin(PI - DEFLECTION)*SINAL
                V_TEMP[2] =    numpy.cos(PI - DEFLECTION)
                V_TEMP[3] =    numpy.sin(PI - DEFLECTION)*COSAL

                #    ! ** Express it now in the (x,y,z) SRF
                V_REF[1] = V_TEMP[1] * U_OLD[1] + V_TEMP[2] * R_OLD[1] + V_TEMP[3]*RP_OLD[1]
                V_REF[2] = V_TEMP[1] * U_OLD[2] + V_TEMP[2] * R_OLD[2] + V_TEMP[3]*RP_OLD[2]
                V_REF[3] = V_TEMP[1] * U_OLD[3] + V_TEMP[2] * R_OLD[3] + V_TEMP[3]*RP_OLD[3]

                #    ! ** This is now the perp. vector in the OMRF **
                V_TEMP[1] = - numpy.sin(3*PIHALF - DEFLECTION)*SINAL
                V_TEMP[2] =   numpy.cos(3*PIHALF - DEFLECTION)
                V_TEMP[3] =   numpy.sin(3*PIHALF - DEFLECTION)*COSAL

                #    ! ** Rotate it to the SRF\
                V_PERP[1] = V_TEMP[1]*U_OLD[1] + V_TEMP[2]*R_OLD[1] + V_TEMP[3]*RP_OLD[1]
                V_PERP[2] = V_TEMP[1]*U_OLD[2] + V_TEMP[2]*R_OLD[2] + V_TEMP[3]*RP_OLD[2]
                V_PERP[3] = V_TEMP[1]*U_OLD[3] + V_TEMP[2]*R_OLD[3] + V_TEMP[3]*RP_OLD[3]

                #    ! ** This is the tangent vector in the OMRF **
                V_TEMP[1] = - numpy.sin(PIHALF - T_INCIDENCE)*SINAL
                V_TEMP[2] =   numpy.cos(PIHALF - T_INCIDENCE)
                V_TEMP[3] =   numpy.sin(PIHALF - T_INCIDENCE)*COSAL

                #    ! ** Rotate it to the SRF.
                V_VEC[1] = V_TEMP[1] * U_OLD[1] + V_TEMP[2] * R_OLD[1] + V_TEMP[3] * RP_OLD[1]
                V_VEC[2] = V_TEMP[1] * U_OLD[2] + V_TEMP[2] * R_OLD[2] + V_TEMP[3] * RP_OLD[2]
                V_VEC[3] = V_TEMP[1] * U_OLD[3] + V_TEMP[2] * R_OLD[3] + V_TEMP[3] * RP_OLD[3]

                #    ! ** Last, we generate U_VEC in the OMRF **
                V_TEMP[1] =   COSAL
                V_TEMP[2] =   .0
                V_TEMP[3] =   SINAL

                #    ! ** rotate to SRF
                U_VEC[1] = V_TEMP[1] * U_OLD[1] + V_TEMP[2] * R_OLD[1] + V_TEMP[3] * RP_OLD[1]
                U_VEC[2] = V_TEMP[1] * U_OLD[2] + V_TEMP[2] * R_OLD[2] + V_TEMP[3] * RP_OLD[2]
                U_VEC[3] = V_TEMP[1] * U_OLD[3] + V_TEMP[2] * R_OLD[3] + V_TEMP[3] * RP_OLD[3]

                #    ! ** All done. Write to the array and leave.
                CENTRAL[1]  =   CENTRAL_OLD[7]
                CENTRAL[2]  =   CENTRAL_OLD[8]
                CENTRAL[3]  =   CENTRAL_OLD[9]
                CENTRAL[4]  =   T_SOURCE * R_OLD[1] + CENTRAL[1]
                CENTRAL[5]  =   T_SOURCE * R_OLD[2] + CENTRAL[2]
                CENTRAL[6]  =   T_SOURCE * R_OLD[3] + CENTRAL[3]
                CENTRAL[7]  =   T_IMAGE *  V_REF[1] + CENTRAL[4]
                CENTRAL[8]  =   T_IMAGE *  V_REF[2] + CENTRAL[5]
                CENTRAL[9]  =   T_IMAGE *  V_REF[3] + CENTRAL[6]
                CENTRAL[10] =   U_VEC[1]
                CENTRAL[11] =   U_VEC[2]
                CENTRAL[12] =   U_VEC[3]
                CENTRAL[13] =   V_VEC[1]
                CENTRAL[14] =   V_VEC[2]
                CENTRAL[15] =   V_VEC[3]
                CENTRAL[16] =   W_VEC[1]
                CENTRAL[17] =   W_VEC[2]
                CENTRAL[18] =   W_VEC[3]
                CENTRAL[19] =   V_REF[1]
                CENTRAL[20] =   V_REF[2]
                CENTRAL[21] =   V_REF[3]
                CENTRAL[22] =   V_PERP[1]
                CENTRAL[23] =   V_PERP[2]
                CENTRAL[24] =   V_PERP[3]


            mirr[0,i] = CENTRAL[4]
            mirr[1,i] = CENTRAL[5]
            mirr[2,i] = CENTRAL[6]
            star[0,i] = CENTRAL[7]
            star[1,i] = CENTRAL[8]
            star[2,i] = CENTRAL[9]

        if n_elements > 0:
            x = numpy.append(source[0],mirr[0,:])
            x = numpy.append(x,star[0,-1])
            x = numpy.round(x,10)

            y = numpy.append(source[1],mirr[1,:])
            y = numpy.append(y,star[1,-1])
            y = numpy.round(y,10)

            z = numpy.append(source[2],mirr[2,:])
            z = numpy.append(z,star[2,-1])
            z = numpy.round(z,10)

        # print(">>>>>> source: ", source)
        # print(">>>>>> X", x)
            return {"source":source, "mirr":mirr, "star":star, "optical_axis_x":x, "optical_axis_y":y, "optical_axis_z":z}
        else:
            return {"source": numpy.zeros(3)}

    def sysinfo(self, title="", comment=""):
        """
        Returns the system information (sysinfo in shadow3).

        Returns
        -------
        str
        """
        txt = "\n"

        TOPLIN = '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
        TSCR = {}
        TYPE1 = {}

        TSCR["1"] = 'AFTER Mirror'
        TSCR["2"] = 'BEFORE  Mirror'

        TYPE1["1"]  = 'SPHERICAL   '
        TYPE1["2"]  = 'ELLIPTICAL  '
        TYPE1["3"]  = 'TOROIDAL    '
        TYPE1["4"]  = 'PARABOLICAL '
        TYPE1["5"]  = 'PLANE       '
        TYPE1["6"]  = 'CODLING SLIT'
        TYPE1["7"]  = 'HYPERBOLICAL'
        TYPE1["8"]  = 'CONICAL     '
        TYPE1["9"]  = '            '
        TYPE1["10"] = '            '
        TYPE1["11"] = '            '
        TYPE1["12"] = '            '
        BREAK	='    ----------------\n'

        txt += TOPLIN
        txt += '**************  S Y S T E M      D E S C R I P T I O N  **************\n'
        txt += TOPLIN
        txt += title + "\n"
        txt += comment + "\n"
        txt += TOPLIN

        n_elements = self.get_beamline_elements_number()



        for i in range(n_elements):
            ble = self.get_beamline_element_at(i)
            coor = ble.get_coordinates()
            oe = ble.get_optical_element()

            txt += ' \n'
            txt += 'Optical Element # %d\n'%(i + 1)
            if isinstance(oe, S4Empty):
                TEXT = "EMPTY ELEMENT"
            elif isinstance(oe, S4Mirror):
                TEXT = "MIRROR"
            elif isinstance(oe, S4Crystal):
                TEXT = "CRYSTAL"
            elif isinstance(oe, S4Grating):
                TEXT = "GRATING"
            elif isinstance(oe, S4Lens):
                TEXT = "LENS"
            elif isinstance(oe, S4CRL):
                TEXT = "COMPOUND REFRACTIVE LENS"
            elif isinstance(oe, S4Transfocator):
                TEXT = "TRANSFOCATOR"
            elif isinstance(oe, S4Screen):
                TEXT = "SCREEN/SLIT/ABSORBER"
            else:
                TEXT = "*** Error: NOT IDENTIFIED OE %s ***" % repr(oe)

            T_SOURCE, T_IMAGE, T_INCIDENCE, T_REFLECTION, ALPHA = coor.get_positions()


            txt +=  '\n'
            txt +=  TEXT+"\n"
            txt +=  '\n'
            txt +=  '  Orientation        %f deg\n'  %(numpy.degrees(ALPHA))
            txt +=  '  Source Plane       %f m  \n'   %(T_SOURCE)
            txt +=  '  Incidence Ang.     %f deg (grazing: %f mrad)\n'  %( numpy.degrees(T_INCIDENCE), (numpy.pi / 2 - T_INCIDENCE) * 1e3)
            txt +=  '  Reflection Ang.    %f deg (grazing: %f mrad)\n'  %( numpy.degrees(T_REFLECTION),(numpy.pi / 2 - T_REFLECTION) * 1e3)
            txt +=  '  Image Plane        %f m  \n'    %(T_IMAGE)
            txt +=  '  O.E. (inter) thickness        %f m  \n' % (oe.interthickness())
            txt += 	BREAK

        txt += "\n\n                          OPTICAL SYSTEM CONFIGURATION\n"
        txt += "                           Laboratory Reference Frame.\n"

        txt += "OPT. Elem #       X =                 Y =                 Z =\n\n"


        dic = self.syspositions()
        txt += "       0    %18.11f     %18.11f     %18.11f\n"%(dic["source"][0],dic["source"][1],dic["source"][2])
        for i in range(n_elements):
            txt += "       %d    %18.11f     %18.11f     %18.11f\n"%(    i+1,dic["mirr"][0,i],dic["mirr"][1,i],dic["mirr"][2,i])
            txt += "          %d'    %18.11f     %18.11f     %18.11f\n"%(i+1,dic["star"][0,i],dic["star"][1,i],dic["star"][2,i])


        txt += TOPLIN
        txt += '********                 E N D                  ***************\n'
        txt += TOPLIN


        return(txt)

    def distances_summary(self,file=''):
        """
        Returns a summary of the real distances, focal distances and orientation angles.

        Returns
        -------
        str
        """
        txt = '  ********  SUMMARY OF DISTANCES ********\n'
        txt += '   ** DISTANCES FOR ALL O.E. [m] **           \n'
        txt += "%12s %12s %14s %14s %14s %14s \n"%('OE','TYPE','p[m]','q[m]','src-oe','src-screen')

        tot = 0.0
        alphatot=0.0
        deflection_H = 0.0
        deflection_V = 0.0
        pihalf = numpy.pi/2
        txt1 = ''
        txt2 = ''
        oeshape = '?'

        n_elements = self.get_beamline_elements_number()

        for i in range(n_elements):
            ble = self.get_beamline_element_at(i)
            coor = ble.get_coordinates()
            oe = ble.get_optical_element()
            T_SOURCE, T_IMAGE, T_INCIDENCE, T_REFLECTION, ALPHA = coor.get_positions()
            is_focusing = 0

            # txt += ' \n'
            # txt += 'Optical Element # %d\n'%(i + 1)

            if isinstance(oe, S4Empty):
                oetype = "EMPTY ELEMENT"
            elif isinstance(oe, S4Mirror):
                oetype = "MIRROR"
                is_focusing = 1
            elif isinstance(oe, S4Crystal):
                oetype = "CRYSTAL"
            elif isinstance(oe, S4Grating):
                oetype = "GRATING"
            elif isinstance(oe, S4Lens):
                oetype = "LENS"
                is_focusing = 1
            elif isinstance(oe, S4CRL):
                oetype = "COMPOUND REFRACTIVE LENS"
                is_focusing = 1
            elif isinstance(oe, S4Transfocator):
                oetype = "TRANSFOCATOR"
                is_focusing = 1
            elif isinstance(oe, S4Screen):
                oetype = "SCR/SLIT/ABS"
            else:
                oetype = "*** Error: NOT IDENTIFIED OE %s ***" % repr(oe)

            #1) Distances summary
            tot = tot + T_SOURCE + oe.interthickness() + T_IMAGE
            totoe = tot - T_IMAGE
            line="%12d %12s %14.4f %14.4f %14.4f %14.4f \n"%(i+1,oetype,T_SOURCE,T_IMAGE,totoe,tot)
            txt1 += line

            # 2) focusing summary
            try:
                sh1 = oe.get_surface_shape_instance()
                oeshape =  sh1.__class__.__name__
            except:
                oeshape = "?????"


            got_focal_distances = 0

            try:
                sc = oe._surface_calculation
                if sc == 0:
                    pp = sh1.get_p_focus()
                    qq = sh1.get_q_focus()
                    got_focal_distances = 1
            except:
                pass

            if is_focusing:
                if got_focal_distances:
                    line = '%10d %20s %10.2f %10.2f %10.2f \n' % (i + 1, oeshape, pp, qq, pp / qq)
                else:
                    line = '%10d %20s %10s %10s %10s \n'%( i+1,oeshape,'?','?','?')

                txt2 += line

            T_INCIDENCE  = T_INCIDENCE
            T_REFLECTION = T_REFLECTION
            ALPHA        = ALPHA

            # 3) total deflection
            alphatot = alphatot + ALPHA
            deflection_H = deflection_H +  numpy.sin(alphatot) *  ( (pihalf-T_INCIDENCE) + (pihalf-T_REFLECTION) )
            deflection_V = deflection_V +  numpy.cos(alphatot) *  ( (pihalf-T_INCIDENCE) + (pihalf-T_REFLECTION) )

        txt += txt1
        txt += '\n'
        txt += '   ** FOCUSING ELEMENTS **           \n'
        # focusing elements
        line = '%10s %20s %10s %10s %10s \n'%('OE','SHAPE','p_foc','q_foc','1/M')
        txt += line
        txt += txt2

        txt += '\n'
        line = 'Sum of Alphas [deg]:         %.3f \n'%(alphatot*180/numpy.pi)
        txt += line
        line = 'Sum of Alphas Mod 180 [deg]: %.3f \n'%( numpy.mod(alphatot*180/numpy.pi,180))
        txt += line
        line = 'Sum of Alphas Mod 360 [deg]: %.3f \n'%( numpy.mod(alphatot*180/numpy.pi,360))
        txt += line
        txt += '\n'

        line = 'Total deflection angle H = %12.6f rad = %9.3f deg\n'%(deflection_H,deflection_H*180/numpy.pi)
        txt += line
        line = 'Total deflection angle V = %12.6f rad = %9.3f deg \n'%(deflection_V,deflection_V*180/numpy.pi)
        txt += line

        return(txt)

if __name__ == "__main__":
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror, S4PlaneMirrorElement
    from shadow4.beamline.optical_elements.refractors.s4_transfocator import S4Transfocator, S4TransfocatorElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from srxraylib.plot.gol import set_qt
    set_qt()

    light_source = SourceGeometrical(name='SourceGeometrical', nrays=10000, seed=5676561)

    m1 = S4PlaneMirror()
    # m2 = S4PlaneMirror()
    m2 = S4Transfocator()


    e1 = S4PlaneMirrorElement(m1, ElementCoordinates())
    # e2 = S4PlaneMirrorElement(m2, ElementCoordinates())
    e2 = S4TransfocatorElement(m2, ElementCoordinates())

    bl = S4Beamline(light_source=light_source) # , beamline_elements_list=[e1,e2])

    # print(bl.info())
    #
    # print(bl.to_python_code())
    #
    # output_beam, output_mirr = bl.run_beamline()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = output_beam.get_rays()
    # plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(Xp,Zp) in microns')

    # print(bl.sourcinfo())

    # print(bl.oeinfo())
    #
    # print(bl.sysinfo())
    #
    # print(bl.to_json())

    # print(bl.syspositions())

    print(bl.distances_summary())

    print(">>>> oe info: ", bl.oeinfo())

    print(">>>> sys info: ", bl.sysinfo())

