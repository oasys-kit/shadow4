"""
Read/write shadow4 objects (e.g., S4Beamline) in json files.

The json files are written by the SYNED mechanism (SynedObject.to_json(), e.g., S4Beamline.to_json()), with two
additions for shadow4:

* the arguments passed to the constructor of every shadow4 optical element (recorded in its attribute
  _s4_init_kwargs, see S4OpticalElementDecorator) are written in the key "s4_init_kwargs" of the optical element;
* the movements of the beamline elements are written in the key "movements" of the beamline element.

When reading, the shadow4 optical elements are rebuilt by calling their constructor with these arguments, so that all
the quantities derived in the constructor (e.g., the optical surface used for the ray tracing, or the inputs used by
to_python_code) are correctly restored; the other objects are read with the SYNED mechanism (default constructor and
setting the values of the keys). Files written by previous versions (without "s4_init_kwargs") are also read: the
optical surface of curved elements is then set from the SYNED surface shape stored in the file.
"""
import inspect

import json_tricks
from urllib.request import urlopen

from syned.syned_object import SynedObject
from syned.util.json_tools import get_exec_commands_for_package
import syned
import shadow4

from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

#
# writing
#
def movements_to_dictionary(movements):
    """Returns a json-compatible dictionary with an S4BeamlineElementMovements."""
    return {"CLASS_NAME": "S4BeamlineElementMovements",
            "f_move"    : movements.f_move,
            "offset_x"  : movements.offset_x,
            "offset_y"  : movements.offset_y,
            "offset_z"  : movements.offset_z,
            "rotation_x": movements.rotation_x,
            "rotation_y": movements.rotation_y,
            "rotation_z": movements.rotation_z,
            }

def _dabax_to_dictionary(dabax):
    params = {}
    for name in list(inspect.signature(type(dabax).__init__).parameters)[1:]:
        if name == "verbose": continue
        getter = getattr(dabax, "get_" + name, None)
        if getter is not None: params[name] = getter()
    return {"CLASS_NAME": type(dabax).__name__, "s4_dabax": params}

def encode_value(value):
    """
    Returns a json-compatible version of a value (e.g., an argument of the constructor of a shadow4 optical element).
    Raises TypeError if the value cannot be written.
    """
    if value is None or isinstance(value, (bool, int, float, complex, str)):
        return value
    if isinstance(value, SynedObject):
        dictionary = value.to_dictionary()
        init_kwargs = getattr(value, "_s4_init_kwargs", None)
        if init_kwargs is not None:
            dictionary["s4_init_kwargs"] = encode_value(init_kwargs)
        return dictionary
    if isinstance(value, S4BeamlineElementMovements):
        return movements_to_dictionary(value)
    if type(value).__module__.startswith("dabax"):
        return _dabax_to_dictionary(value)
    if isinstance(value, dict):
        return {key: encode_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode_value(item) for item in value]
    if type(value).__module__ == "numpy":   # arrays and numpy scalars (written by json_tricks)
        return value
    raise TypeError("cannot write in json an object of type %s" % type(value).__name__)

def add_shadow4_items_to_dictionary(dictionary, beamline_element):
    """
    Adds to the dictionary of a beamline element (from SynedObject.to_dictionary) the constructor arguments of its
    optical element ("s4_init_kwargs") and its movements ("movements").
    """
    optical_element = beamline_element.get_optical_element()
    init_kwargs = getattr(optical_element, "_s4_init_kwargs", None)
    if init_kwargs is not None and isinstance(dictionary.get("optical_element", None), dict):
        try:
            dictionary["optical_element"]["s4_init_kwargs"] = encode_value(init_kwargs)
        except TypeError as e:
            print("Warning: the constructor arguments of %s are not written in json (%s)" % (
                  type(optical_element).__name__, str(e)))
    movements = beamline_element.get_movements()
    if movements is not None:
        dictionary["movements"] = movements_to_dictionary(movements)
    return dictionary

#
# reading
#
def _namespace(exec_commands=None):
    ns = {}
    commands = [get_exec_commands_for_package(syned), get_exec_commands_for_package(shadow4)]
    if exec_commands is not None: commands.append(exec_commands)
    for command in commands:
        exec(command, ns)
    try:
        exec("from dabax.dabax_xraylib import DabaxXraylib", ns)
    except Exception:
        pass
    return ns

def _class(class_name, ns):
    try:
        return eval(class_name, ns)
    except Exception:
        raise RuntimeError("Error evaluating: %s ** you could use the exec_commands keyword to import it **" % class_name)

def _decode(value, ns):
    if isinstance(value, dict):
        if "CLASS_NAME" in value: return load_from_dictionary(value, ns)
        return {key: _decode(item, ns) for key, item in value.items()}
    if isinstance(value, list):
        return [_decode(item, ns) for item in value]
    return value

def load_from_dictionary(jsn, ns):
    """Returns the object (shadow4 or SYNED) stored in a dictionary (from a json file)."""
    class_name = jsn["CLASS_NAME"]
    cls = _class(class_name, ns)

    if "s4_dabax" in jsn:
        return cls(**jsn["s4_dabax"])
    if class_name == "S4BeamlineElementMovements":
        return S4BeamlineElementMovements(**{key: value for key, value in jsn.items() if key != "CLASS_NAME"})

    init_kwargs = jsn.get("s4_init_kwargs", None)
    if isinstance(init_kwargs, dict):
        # shadow4 optical element: rebuilt with its constructor; then the values of simple keys are set (they may
        # have been modified after the construction)
        obj = cls(**_decode(init_kwargs, ns))
        if obj.keys() is not None:
            for key in obj.keys():
                if key in jsn and not isinstance(jsn[key], (dict, list)):
                    obj.set_value_from_key_name(key, jsn[key])
        return obj

    # SYNED mechanism: default constructor and values of the keys
    obj = cls()
    if obj.keys() is not None:
        for key in obj.keys():
            if key in jsn:
                obj.set_value_from_key_name(key, _decode(jsn[key], ns))
    if "movements" in jsn and hasattr(obj, "set_movements"):
        obj.set_movements(_decode(jsn["movements"], ns))
    # files of previous versions: the optical surface of curved shadow4 elements is the SYNED surface shape
    if hasattr(obj, "_curved_surface_shape") and getattr(obj, "_surface_shape", None) is not None:
        obj._curved_surface_shape = obj._surface_shape
    return obj

def load_from_json_text(text, exec_commands=None):
    """
    Load a syned/shadow4 object from a json string.

    Parameters
    ----------
    text : str
        The json-encoded object as a string.
    exec_commands : str, optional
        Additional import commands (e.g. for classes outside of syned and shadow4).

    Returns
    -------
    instance of SynedObject
    """
    return load_from_dictionary(json_tricks.loads(text), _namespace(exec_commands))

def load_from_json_file(file_name, exec_commands=None):
    """
    Load a syned/shadow4 object from a json file.

    Parameters
    ----------
    file_name : str
        Path to the json file.
    exec_commands : str, optional
        Additional import commands (e.g. for classes outside of syned and shadow4).

    Returns
    -------
    instance of SynedObject
    """
    with open(file_name) as f:
        text = f.read()
    return load_from_json_text(text, exec_commands=exec_commands)

def load_from_json_url(file_url, exec_commands=None):
    """
    Load a syned/shadow4 object from a remote json file.

    Parameters
    ----------
    file_url : str
        URL of the remote json file.
    exec_commands : str, optional
        Additional import commands (e.g. for classes outside of syned and shadow4).

    Returns
    -------
    instance of SynedObject
    """
    u = urlopen(file_url)
    return load_from_json_text(u.read().decode(encoding='UTF-8'), exec_commands=exec_commands)


if __name__ == "__main__":
    file_url = "/home/srio/Oasys2/tmp0.json"
    syned_obj = load_from_json_file(file_url)
    print(syned_obj)
    print(syned_obj.info())
