import logging

def is_verbose():
    return logging.root.level <= logging.INFO

def is_debug():
    return logging.root.level <= logging.DEBUG

def set_verbose(status=1):
    if status:
        logging.getLogger().setLevel(logging.INFO)
    else:
        logging.getLogger().setLevel(logging.WARNING)

def set_debug(status):
    if status:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARNING)

def printlog(*args):
    logging.info(*args)

if __name__ == "__main__":
    from shadow4.tools.logger import is_verbose, is_debug

    print("is_verbose:", is_verbose())
    logging.info("This is not seen")
    set_verbose(1)
    print("is_verbose:", is_verbose())
    logging.info("This is \n seen")
    printlog("via pringlog...")
    print("is_debug; ", is_debug())

