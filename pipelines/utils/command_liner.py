import os
import sys

def cmd_exit(cmd_string, task):
    return exit_if_error(os.system(cmd_string), task)

def exit_if_error(code, task):
    if code == 0:
        return
    elif code == 256:
        print("No action was performed")
    else:
        print("Error occurred while performing task: ", task, code)
        sys.exit(code)
