"""
Adapt to Shadow4 the SYNED functions to read/write syned objects in json files.

Functions here wrap the syned equivalents and automatically inject all shadow4
class imports as exec_commands, so callers do not need to supply them manually.
"""

from syned.util.json_tools import load_from_json_file as SY_load_from_json_file
from syned.util.json_tools import load_from_json_text as SY_load_from_json_text
from syned.util.json_tools import get_exec_commands_for_package

import shadow4

from urllib.request import urlopen

def load_from_json_file(file_name, exec_commands=None):
    """
    Load a syned/shadow4 object from a json file.

    All shadow4 class imports are added automatically via
    :func:`get_shadow4_exec_commands`, so no manual exec_commands are needed
    for typical shadow4 json files.

    Parameters
    ----------
    file_name : str
        Path to the json file.
    exec_commands : str, optional
        Additional import commands to prepend before the shadow4 ones.

    Returns
    -------
    instance of SynedObject
    """
    if exec_commands is None:
        exec_commands = get_exec_commands_for_package(shadow4)
    else:
        exec_commands = exec_commands + "\n" + get_exec_commands_for_package(shadow4)
    return SY_load_from_json_file(file_name, exec_commands=exec_commands)

def load_from_json_text(text, exec_commands=None):
    """
    Load a syned/shadow4 object from a json string.

    All shadow4 class imports are added automatically via
    :func:`get_shadow4_exec_commands`, so no manual exec_commands are needed
    for typical shadow4 json strings.

    Parameters
    ----------
    text : str
        The json-encoded object as a string.
    exec_commands : str, optional
        Additional import commands to prepend before the shadow4 ones.

    Returns
    -------
    instance of SynedObject
    """
    if exec_commands is None:
        exec_commands = get_exec_commands_for_package(shadow4)
    else:
        exec_commands = exec_commands + "\n" + get_exec_commands_for_package(shadow4)
    return SY_load_from_json_text(text, exec_commands=exec_commands)

def load_from_json_url(file_url, exec_commands=None):
    """
    Load a syned/shadow4 object from a remote json file.

    All shadow4 class imports are added automatically via
    :func:`get_shadow4_exec_commands`, so no manual exec_commands are needed
    for typical shadow4 json files.

    Parameters
    ----------
    file_url : str
        URL of the remote json file.
    exec_commands : str, optional
        Additional import commands to prepend before the shadow4 ones.

    Returns
    -------
    instance of SynedObject
    """
    if exec_commands is None:
        exec_commands = get_exec_commands_for_package(shadow4)
    else:
        exec_commands = exec_commands + "\n" + get_exec_commands_for_package(shadow4)

    u = urlopen(file_url)
    ur = u.read()
    url = ur.decode(encoding='UTF-8')
    return load_from_json_text(url, exec_commands=exec_commands)


if __name__ == "__main__":
    file_url = "/home/srio/Oasys2/tmp0.json"
    syned_obj = load_from_json_file(file_url)
    print(syned_obj)
    print(syned_obj.info())
    # syned_obj.run_beamline()
