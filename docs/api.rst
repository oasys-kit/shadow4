.. currentmodule:: shadow4

===========
Package API
===========
This page lists main classes in this package.


beam
----
``shadow4.beam`` The S4 beam

.. autosummary::
   :toctree: generated/

   shadow4.beam.s4_beam

sources
-------
``shadow4.sources`` The S4 sources

* ``shadow4.sources.source_geometrical`` The S4 geometrical sources

   .. autosummary::
      :toctree: generated/

      shadow4.sources.source_geometrical.source_geometrical
      shadow4.sources.source_geometrical.source_grid_cartesian
      shadow4.sources.source_geometrical.source_grid_polar
      shadow4.sources.source_geometrical.source_gaussian

* ``shadow4.sources.bending_magnet`` The S4 bending magnet

   .. autosummary::
      :toctree: generated/

      shadow4.sources.bending_magnet.s4_bending_magnet
      shadow4.sources.bending_magnet.s4_bending_magnet_light_source

* ``shadow4.sources.wiggler`` The S4 wiggler

   .. autosummary::
      :toctree: generated/

      shadow4.sources.wiggler.s4_wiggler
      shadow4.sources.wiggler.s4_wiggler_light_source

* ``shadow4.sources.undulator`` The S4 undulators

   .. autosummary::
      :toctree: generated/

      shadow4.sources.undulator.s4_undulator_gaussian
      shadow4.sources.undulator.s4_undulator
      shadow4.sources.undulator.s4_undulator_gaussian_light_source
      shadow4.sources.undulator.s4_undulator_light_source

beamline optical elements
-------------------------
``shadow4.beamline.optical_elements`` The S4 optical elements

* ``shadow4.sources.optical_elements.absorbers`` The absorbers, slits, sops, screens

   .. autosummary::
      :toctree: generated/

      shadow4.beamline.optical_elements.absorbers.s4_screen

* ``shadow4.sources.optical_elements.mirrors`` The mirrors

   .. autosummary::
      :toctree: generated/

      shadow4.beamline.optical_elements.mirrors.s4_plane_mirror
      shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror
      shadow4.beamline.optical_elements.mirrors.s4_toroid_mirror
      shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror
      shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror
      shadow4.beamline.optical_elements.mirrors.s4_hyperboloid_mirror
      shadow4.beamline.optical_elements.mirrors.s4_conic_mirror
      shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror
      shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror

* ``shadow4.sources.optical_elements.gratings`` The gratings

   .. autosummary::
      :toctree: generated/

      shadow4.beamline.optical_elements.gratings.s4_plane_grating
      shadow4.beamline.optical_elements.gratings.s4_sphere_grating
      shadow4.beamline.optical_elements.gratings.s4_toroid_grating
      shadow4.beamline.optical_elements.gratings.s4_ellipsoid_grating
      shadow4.beamline.optical_elements.gratings.s4_paraboloid_grating
      shadow4.beamline.optical_elements.gratings.s4_hyperboloid_grating
      shadow4.beamline.optical_elements.gratings.s4_conic_grating
      shadow4.beamline.optical_elements.gratings.s4_numerical_mesh_grating
      shadow4.beamline.optical_elements.gratings.s4_additional_numerical_mesh_grating

* ``shadow4.sources.optical_elements.crystals`` The crystals

   .. autosummary::
      :toctree: generated/

      shadow4.beamline.optical_elements.crystals.s4_plane_crystal
      shadow4.beamline.optical_elements.crystals.s4_sphere_crystal
      shadow4.beamline.optical_elements.crystals.s4_toroid_crystal
      shadow4.beamline.optical_elements.crystals.s4_ellipsoid_crystal
      shadow4.beamline.optical_elements.crystals.s4_paraboloid_crystal
      shadow4.beamline.optical_elements.crystals.s4_hyperboloid_crystal
      shadow4.beamline.optical_elements.crystals.s4_conic_crystal
      shadow4.beamline.optical_elements.crystals.s4_numerical_mesh_crystal
      shadow4.beamline.optical_elements.crystals.s4_additional_numerical_mesh_crystal

* ``shadow4.sources.optical_elements.refractors`` The refractors

   .. autosummary::
      :toctree: generated/

      shadow4.beamline.optical_elements.refractors.s4_lens
      shadow4.beamline.optical_elements.refractors.s4_crl
      shadow4.beamline.optical_elements.refractors.s4_conic_interface

* ``shadow4.sources.optical_elements.ideal_elements`` The ideal elements

   .. autosummary::
      :toctree: generated/

      shadow4.beamline.optical_elements.ideal_elements.s4_beam_movement
      shadow4.beamline.optical_elements.ideal_elements.s4_empty
      shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens