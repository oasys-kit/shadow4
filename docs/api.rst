.. currentmodule:: shadow4

===========
Package API
===========
This page lists main classes in this package.


``shadow4`` classes.

main
----
``shadow4`` Base class

.. autosummary::
   :toctree: generated/

   shadow4_object

beamline
--------
``shadow4.beamline`` A beamline is made by a light_source and one or several beamline_element

.. autosummary::
   :toctree: generated/

   shadow4.beamline.s4_beamline

light source
------------
Light source = Electron beam + magnetic structure

* ``shadow4.sources.s4_light_source`` Electron beam + magnetic structure

   .. autosummary::
      :toctree: generated/

      shadow4.sources.s4_light_source

* ``shadow4.sources.s4_light_source_base`` Base light source

   .. autosummary::
      :toctree: generated/

      shadow4.sources.s4_light_source_base

eletron beam
-------------
Electron beam

``shadow4.sources.s4_electron_beam`` Electron beam

.. autosummary::
   :toctree: generated/

   shadow4.sources.s4_electron_beam

magnetic structures
-------------------
Magnetic structures


``shadow4.sources`` Sources

.. autosummary::
   :toctree: generated/

   shadow4.sources
   shadow4.sources.source_geometrical
   shadow4.sources.bending_magnet
   shadow4.sources.wiggler
   shadow4.sources.undulator



beamline element
----------------
``shadow4.beamline.s4_beamline_element`` A beamline element is made by an optical_element and the element_coordinates

.. autosummary::
   :toctree: generated/

   shadow4.beamline.s4_beamline_element


optical elements
----------------

``shadow4.beamline.optical_elements`` Optical elements

* ``shadow4.beamline.optical_elements.absorbers`` absorber, beam_stopper, filter, holed_filter, slit
.. autosummary::
   :toctree: generated/

   shadow4.beamline.optical_elements.absorbers.s4_screen

* ``shadow4.beamline.optical_elements.crystals`` crystals

.. autosummary::
   :toctree: generated/

   shadow4.beamline.optical_elements.crystals
   shadow4.beamline.optical_elements.crystals.s4_crystal
   shadow4.beamline.optical_elements.crystals.s4_plane_crystal
   shadow4.beamline.optical_elements.crystals.s4_sphere_crystal
   shadow4.beamline.optical_elements.crystals.s4_toroid_crystal
   shadow4.beamline.optical_elements.crystals.s4_conic_crystal
   shadow4.beamline.optical_elements.crystals.s4_ellipsoid_crystal
   shadow4.beamline.optical_elements.crystals.s4_hyperboloid_crystal
   shadow4.beamline.optical_elements.crystals.s4_paraboloid_crystal
   shadow4.beamline.optical_elements.crystals.s4_numerical_mesh_crystal
   shadow4.beamline.optical_elements.crystals.s4_additional_numerical_mesh_crystal


* ``shadow4.beamline.optical_elements.gratings`` gratings

.. autosummary::
   :toctree: generated/

   shadow4.beamline.optical_elements.gratings
   shadow4.beamline.optical_elements.gratings.s4_grating
   shadow4.beamline.optical_elements.gratings.s4_plane_grating
   shadow4.beamline.optical_elements.gratings.s4_sphere_grating
   shadow4.beamline.optical_elements.gratings.s4_toroid_grating
   shadow4.beamline.optical_elements.gratings.s4_conic_grating
   shadow4.beamline.optical_elements.gratings.s4_ellipsoid_grating
   shadow4.beamline.optical_elements.gratings.s4_hyperboloid_grating
   shadow4.beamline.optical_elements.gratings.s4_paraboloid_grating
   shadow4.beamline.optical_elements.gratings.s4_numerical_mesh_grating
   shadow4.beamline.optical_elements.gratings.s4_additional_numerical_mesh_grating

* ``shadow4.beamline.optical_elements.mirrors`` mirrors

.. autosummary::
   :toctree: generated/

   shadow4.beamline.optical_elements.mirrors
   shadow4.beamline.optical_elements.mirrors.s4_mirror
   shadow4.beamline.optical_elements.mirrors.s4_plane_mirror
   shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror
   shadow4.beamline.optical_elements.mirrors.s4_toroid_mirror
   shadow4.beamline.optical_elements.mirrors.s4_conic_mirror
   shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror
   shadow4.beamline.optical_elements.mirrors.s4_hyperboloid_mirror
   shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror
   shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror
   shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror

* ``shadow4.beamline.optical_elements.refractors`` interface, lens, crl

.. autosummary::
   :toctree: generated/

   shadow4.beamline.optical_elements.refractors
   shadow4.beamline.optical_elements.refractors.s4_interface
   shadow4.beamline.optical_elements.refractors.s4_conic_interface
   shadow4.beamline.optical_elements.refractors.s4_lens
   shadow4.beamline.optical_elements.refractors.s4_crl

* ``shadow4.beamline.optical_elements.ideal_elements`` empty, ideal_lens, beam_movement

.. autosummary::
   :toctree: generated/

   shadow4.beamline.optical_elements.ideal_elements
   shadow4.beamline.optical_elements.ideal_elements.s4_beam_movement
   shadow4.beamline.optical_elements.ideal_elements.s4_empty
   shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens

optical surfaces
----------------
``shadow4.optical_surfaces`` Internal classes for dealing with optical surfaces

``shadow4.optical_surfaces.s4_optical_surface`` Base class
   .. autosummary::
      :toctree: generated/

      shadow4.optical_surfaces.s4_optical_surface

``shadow4.optical_surfaces.s4_optical_surface.s4_conic`` Conics
   .. autosummary::
      :toctree: generated/

      shadow4.optical_surfaces.s4_optical_surface.s4_conic


``shadow4.optical_surfaces.s4_optical_surface.s4_toroid`` Toroid
   .. autosummary::
      :toctree: generated/

      shadow4.optical_surfaces.s4_optical_surface.s4_toroid

``shadow4.optical_surfaces.s4_optical_surface.s4_mesh`` Mesh
   .. autosummary::
      :toctree: generated/

      shadow4.optical_surfaces.s4_optical_surface.s4_mesh

physical models
----------------
``shadow4.physical_models`` Internal classes for dealing with reflectivity of optical surfaces

``shadow4.physical_models.bragg`` Crystal stuff
   .. autosummary::
      :toctree: generated/

      shadow4.physical_models.bragg

``shadow4.physical_models.mlayer`` Multilayer stuff
   .. autosummary::
      :toctree: generated/

      shadow4.physical_models.mlayer

``shadow4.physical_models.prerefl`` Mirror and lens stuff
   .. autosummary::
      :toctree: generated/

      shadow4.physical_models.prerefl

