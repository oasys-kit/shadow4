This directory contains the optical surfaces in shadow4

Basically there are only three surfaces:

Conic: implements all the conic surfaces (plane, spherical hyperboloid, etc)
Toroid: implement toroidal surfaces
Mesh: implement numeric mesh.

TODO: right now Mesh is an independent optical surface. In shadow3 the numeric mesh
    (slope errors, etc.) are added to a mathematical surface (conic, toroid). It is
    necessary to implement this functionality.


