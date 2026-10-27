# PLAN_SYNED_AI: shadow4 beamlines in syned-json, and their reconstruction in OASYS2

Branch `syned-ai` in shadow4, shadow4tests and OASYS2-SHADOW4 (development paused 2026-10-27 to complete the
polarisation work in branch `polarisers`; `syned-ai` was created from `polarisers`).

## Goal

A shadow4 beamline dumped with `S4Beamline.to_json()` must be:

1. read back by python (`shadow4.tools.json_tools.load_from_json_file`) giving an *identical* beamline
   (same ray tracing, same `to_python_code()`);
2. loaded into an OASYS2 canvas, creating and connecting the corresponding widgets with all their settings
   (and, later, the reverse: export an OASYS2 workspace to syned-json);
3. a common, documented exchange format (syned-json) usable by other tools (e.g. AI agents creating/modifying
   beamlines, other codes reading syned).

## Done (commits in `syned-ai`)

### shadow4

- [x] `S4OpticalElementDecorator.__init_subclass__`: wraps the `__init__` of every optical element class and records
      the arguments of the outermost constructor call in `_s4_init_kwargs`.
- [x] `S4BeamlineElement.to_dictionary()`: adds `"s4_init_kwargs"` (encoded) to the optical element and
      `"movements"` to the beamline element.
- [x] `shadow4/tools/json_tools.py` rewritten:
  - writing: `encode_value` (SynedObject + nested `s4_init_kwargs`, movements, DabaxXraylib via `"s4_dabax"`,
    dict/list, numpy), `add_shadow4_items_to_dictionary` (warns and skips kwargs that cannot be encoded);
  - reading: `load_from_dictionary` rebuilds shadow4 optical elements with `cls(**kwargs)` and overlays simple keys;
    other objects with the syned mechanism; movements restored; old files: `_curved_surface_shape = _surface_shape`;
  - `load_from_json_text/file/url(..., exec_commands=None)`.
- [x] `SourceGeometrical`: keys `"f_color "`, `"f_foher "` (trailing spaces) fixed; `f_phot` added to the keys.
- [x] `s4_beamline.py`: commented self-tests use the shadow4 loader.

Bugs fixed: (1) curved surfaces lost on reload (0 good rays), (2) movements not saved, (3) `__inputs` /
`to_python_code()` wrong after reload, (4) source keys.

### shadow4tests

- [x] `test_library/test_json_roundtrip.py` (4 tests; all fail without the fix):
  - beamline round trip (source in A, coherent; plane mirror xraylib; p/q sphere + movements; sphere with
    elliptical boundary at 30 deg): identical rays, identical code, idempotent json;
  - all optical element classes with a usable default (>45);
  - special elements: CRL, lens, additional numerical mesh, compound with inner movements, DABAX crystal/mirror;
  - files of previous versions (no `s4_init_kwargs`): sphere, ellipsoid, hyperboloid, toroid.

### OASYS2-SHADOW4 (prototype)

- [x] `orangecontrib/shadow4/util/beamline_from_json.py`: `load_beamline`, `build_canvas` (creates nodes and
      "Shadow Data" links), populate functions for `OWGeometrical` and `OWMirror` (plane, sphere/cylinder via
      external radius, boundaries, reflectivity, coordinates, movements).
- [x] Shadow4 menu: "Load Shadow4 beamline from JSON file...".
- [x] Verified headless (`QT_QPA_PLATFORM=offscreen`): 4 nodes, 3 links, rebuilt beamline traces the original
      within 3.5e-11.

## TODO

### A. shadow4 core (json)

- [ ] Decide the merge strategy: merge `syned-ai` into master after `polarisers` is merged (rebase if needed;
      `s4_optical_element_decorators.py` may conflict).
- [ ] Python 3.8 check of the new code (OASYS1 compatibility).
- [ ] Paraboloids in old-format files: the syned `Paraboloid` shape does not hold all the shadow4 parameters
      (cannot be restored exactly). Document, or store the missing parameters in syned.
- [ ] Light sources other than `SourceGeometrical`: check round trip of `S4Undulator(Gaussian)`, `S4Wiggler`,
      `S4BendingMagnet`, `SourceGridCartesian/Polar`, electron beam (`S4ElectronBeam`) — keys, derived state,
      `to_python_code()`. Add them to the test.
- [ ] `S4Beamline.to_json()` of beamlines with `S4Beam`-dependent elements (input beam is not written: OK) and
      with external files (prerefl, bragg, mlayer, mesh): decide whether to embed the data or keep file names.
- [ ] Unify with syned: propose upstream (syned) a hook for `s4_init_kwargs`/movements, so that
      `syned.util.json_tools.load_from_json_file` delegates to the package loader (today syned's generic loader
      ignores the new keys).
- [ ] syned: `get_exec_commands_for_package(syned)` imports `syned/__test/test.py`, which prints an electron beam
      at each load (exclude `__test` modules).
- [ ] Versioning of the json format (e.g. `"s4_json_version": 1`) for future changes.
- [ ] Pre-existing (not json) bugs found: `S4CRL.to_python_code()` and `S4Lens.to_python_code()` fail with
      `prerefl_file=None` (the default); default objects of `S4Compound`, `S4AdditionalNumericalMesh*`, `S4Grating`
      cannot generate code.

### B. OASYS2-SHADOW4 importer (json -> canvas)

- [ ] Sources: `OWBendingMagnet`, `OWUndulator`, `OWUndulatorGaussian`, `OWWiggler`, `OWGrid` (+ electron beam
      settings in `ow_synchrotron_source.py`).
- [ ] Mirrors: all surface shapes (ellipsoid, hyperboloid, paraboloid, toroid, conic, numerical mesh, additional
      numerical mesh), internal (p/q) vs external definitions, convexity, cylinder orientation.
- [ ] Crystals (`OWCrystal`, `OWMosaicCrystal`, `OWPhasePlateCrystal`), multilayers, gratings (VLS).
- [ ] Refractors: lens, CRL, transfocator, refractive interface, ideal lens, ideal FZP.
- [ ] Screen/slits, empty element, beam movements.
- [ ] Compound elements: decide mapping (several widgets or not supported).
- [ ] Layout of the canvas (positions, multiple rows), titles from the element names.
- [ ] Report unsupported elements clearly (today: warnings list).
- [ ] Test (headless) for every widget: populate -> rebuild -> trace identical to the original beamline
      (scratchpad prototype `canvas_test.py` to be converted into a test in the OASYS2-SHADOW4 repository).
- [ ] Optionally, run the beamline after loading.

### C. OASYS2 canvas -> json (reverse)

- [ ] Menu action "Save Shadow4 beamline to JSON file...": collect the chain of shadow4 widgets, build the
      `S4Beamline` from the widgets (`get_lightsource()`, `get_beamline_element_instance()`...), `to_json()`.
- [ ] Round trip test: canvas -> json -> canvas.

### D. AI / external usage

- [ ] Document the json format (keys, `s4_init_kwargs`, `movements`, `s4_dabax`) in the shadow4 docs, with an
      example file.
- [ ] Example: create/modify a beamline json by script (or by an AI agent) and load it in OASYS2.
- [ ] Validation tool: load a json, report unknown classes/keys, trace a test beam.
