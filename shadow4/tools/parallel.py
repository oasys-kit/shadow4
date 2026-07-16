import importlib
import multiprocessing
import os
import pathlib
import sys
import time

import joblib
from joblib import Parallel, delayed


UNSUPPORTED_PARALLEL_SOURCE_MESSAGE = (
    "S4LightSourceFromBeamlines, S4LightSourceFromFile, SourceGridCartesian and SourceGridPolar"
    " are deterministic sources and do not support parallel repetitions. "
)


def run_generated_trace(module_name, module_dir, seed, nrays):
    module_dir = os.path.abspath(module_dir)

    if module_dir not in sys.path:
        sys.path.insert(0, module_dir)

    importlib.invalidate_caches()
    module = importlib.import_module(module_name)
    return module.trace_beamline(seed=seed, nrays=nrays)


def load_runner_module(module_path):
    module_path = pathlib.Path(module_path).resolve()
    module_dir = str(module_path.parent)
    module_name = module_path.stem

    if module_dir not in sys.path:
        sys.path.insert(0, module_dir)

    importlib.invalidate_caches()
    if module_name in sys.modules:
        del sys.modules[module_name]

    return importlib.import_module(module_name)


def get_parallel_runner_prototype(beamline):
    from shadow4.sources.s4_light_source_from_beamlines import S4LightSourceFromBeamlines
    from shadow4.sources.source_geometrical.source_grid_cartesian import SourceGridCartesian
    from shadow4.sources.source_geometrical.source_grid_polar import SourceGridPolar
    from shadow4.sources.s4_light_source_from_file import S4LightSourceFromFile

    light_source = beamline.get_light_source()
    prototype_beamline = beamline

    if isinstance(light_source, (SourceGridCartesian, SourceGridPolar, 
                                 S4LightSourceFromFile, S4LightSourceFromBeamlines)):
        raise ValueError(UNSUPPORTED_PARALLEL_SOURCE_MESSAGE)

    return prototype_beamline


def cpu_info_text():
    lines = [
        "CPU availability:",
        "    os.cpu_count(): %s" % os.cpu_count(),
        "    multiprocessing.cpu_count(): %s" % multiprocessing.cpu_count(),
        "    joblib.cpu_count(): %s" % joblib.cpu_count(),
    ]

    try:
        lines.append("\n    CPU affinity: %s" % len(os.sched_getaffinity(0)))
    except AttributeError:
        lines.append("    CPU affinity: not available on this OS")

    return "\n".join(lines)


def print_cpu_info():
    text = cpu_info_text()
    print(text)
    return text


def runner_module_path(script_file):
    if script_file is None:
        return "<interactive Shadow4 script>"
    return os.path.abspath(script_file)


def runner_module_dir(script_file):
    if script_file is None:
        return None
    return os.path.dirname(os.path.abspath(script_file))


def runner_module_name(script_file):
    if script_file is None:
        return None
    return os.path.splitext(os.path.basename(script_file))[0]


def seed_for_iteration(base_seed, iteration):
    if base_seed == 0:
        return 0
    return int(base_seed + iteration * 2)


def concatenate_beams(beam_list, footprint_list, seed_list, verbose=True):
    ntimes = len(seed_list)
    for i in range(ntimes):
        beam_list[i].clean_lost_rays()
        if footprint_list[i] is not None:
            footprint_list[i].clean_lost_rays()

        if i == 0:
            beam_acc = beam_list[i].duplicate()
            footprint_acc = None if footprint_list[i] is None else footprint_list[i].duplicate()
        else:
            beam_acc.append_beam(beam_list[i])
            if footprint_acc is not None and footprint_list[i] is not None:
                footprint_acc.append_beam(footprint_list[i])

        if verbose:
            print("Iteration %d: seed=%d, rays=%d" % (i, seed_list[i], beam_list[i].N))

    return beam_acc, footprint_acc


def build_accumulated_beamline(trace_beamline, seed_list, number_of_rays):
    from shadow4.beamline.s4_beamline import S4Beamline
    from shadow4.sources.s4_light_source_from_beamlines import S4LightSourceFromBeamlines

    light_source_acc = S4LightSourceFromBeamlines(name="Accumulate Parallel Run")

    for seed in seed_list:
        _, _, _, beamline = trace_beamline(
            seed=seed,
            nrays=number_of_rays,
            dry_run=True,
            return_beamline=True,
        )
        light_source_acc.append_beamline(
            beamline,
            id="beamline seed: %d" % seed,
            weight=1.0,
        )

    beamline_acc = S4Beamline()
    beamline_acc.set_light_source(light_source_acc)
    return beamline_acc


def run_parallel_from_generated_script(
    trace_beamline,
    script_file,
    number_of_repetitions,
    number_of_rays,
    n_jobs,
    base_seed,
):
    t_total = time.perf_counter()

    print_cpu_info()
    print("")
    print("Number of repetitions:", number_of_repetitions)
    print("Number of rays:", number_of_rays)
    print("Base seed:", base_seed)
    if n_jobs == -1:
        n_jobs = joblib.cpu_count()
    print("Number of cores:", n_jobs)

    print("")
    print("Runner module:", runner_module_path(script_file))
    print("")

    seed_list = [seed_for_iteration(base_seed, i) for i in range(number_of_repetitions)]

    t_parallel = time.perf_counter()
    module_dir = runner_module_dir(script_file)
    module_name = runner_module_name(script_file)
    if module_dir is None:
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(trace_beamline)(seed=seed, nrays=number_of_rays)
            for seed in seed_list
        )
    else:
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(run_generated_trace)(module_name, module_dir, seed, number_of_rays)
            for seed in seed_list
        )
    parallel_elapsed = time.perf_counter() - t_parallel

    seed_list = [result[0] for result in results]
    beam_list = [result[1] for result in results]
    footprint_list = [result[2] for result in results]

    t_concatenate = time.perf_counter()
    beam_acc, footprint_acc = concatenate_beams(
        beam_list,
        footprint_list,
        seed_list,
        verbose=True,
    )
    beamline_acc = build_accumulated_beamline(trace_beamline, seed_list, number_of_rays)
    concatenate_elapsed = time.perf_counter() - t_concatenate

    print("")
    print("Parallel elapsed: %.3f s" % parallel_elapsed)
    print("Concatenation elapsed: %.3f s" % concatenate_elapsed)
    print("Total elapsed: %.3f s" % (time.perf_counter() - t_total))
    print("Accumulated rays:", beam_acc.N)

    return seed_list, beamline_acc, beam_acc, footprint_acc
