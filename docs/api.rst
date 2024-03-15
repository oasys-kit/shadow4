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

      shadow4.sources.source_geometrical.s4_source_geometrical
      shadow4.sources.source_geometrical.s4_source_grid_cartesian
      shadow4.sources.source_geometrical.s4_source_grid_polar
      shadow4.sources.source_geometrical.s4_source_gaussian

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