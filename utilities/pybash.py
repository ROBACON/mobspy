#!/usr/bin/env python3

#####################################################
#                                                   #
#   Run bash commands through python                #
#                                                   #
#####################################################

import subprocess, traceback
import os, sys


# TODO: clean out and err output
# also in subprocess.run() should pass diff default file for stdout and stderr


def run(command, logfile=None, cwd=None):
    if logfile is None:
        err = try_subprocess(command, sys.stderr, cwd=cwd)
    else:
        with open(logfile, 'a') as file:
            # TODO is there a cleaner way of handling sys.stderr vs custom file
            err = try_subprocess(command, file, cwd=cwd)
    return err


def try_subprocess(command, file, cwd=None):
    err = False

    try:
        subprocess.run(
            command,
            stdin=sys.stdin,
            stdout=file,
            stderr=file,
            check=True,
            cwd=cwd,
        )

    except subprocess.CalledProcessError as e:
        err = True
        print(
            "\n^ ERROR ABOVE ^\nThe above error occured during execution of the subprocess, the wrapper stacktrace is:",
            file=file,
        )
        traceback.print_tb(e.__traceback__, file=file)
        print("Error:", e, file=file)
    finally:
        pass
    return err
