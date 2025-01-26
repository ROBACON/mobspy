from copy import deepcopy
import inspect
import sys
import traceback
import mobspy.modules.mobspy_expressions as me

global_simlog_level = 3

# For debugging
_PRINT_FULL_EXCEPTION_LOG = 0

def error(message, stack_index=-1, full_exception_log = False):

    if stack_index > -1:
        # Get information about the code line where the error occurred
        stack = inspect.stack()
        code_line = stack[stack_index].code_context[0][:-1]
        line_number = stack[stack_index].lineno
        message = f'At: {code_line} \n' + f'Line number: {line_number} \n' + message

    # Print the full traceback if the flag is set
    if _PRINT_FULL_EXCEPTION_LOG or full_exception_log:
        traceback_details = ''.join(traceback.format_stack())
        message = f'Full Exception Log:\n{traceback_details}\n{message}'

    if global_simlog_level >= 0:
        print('\033[91m' + 'ERROR: ' + message + '\033[0m', file=sys.stderr)

    sys.exit(1)


def debug(message=''):
    global global_simlog_level
    if global_simlog_level >= 2:
        message_copy = deepcopy(message)
        print(message_copy.replace('_dot_', '.'), file=sys.stderr)
    pass


def warning(message):
    global global_simlog_level
    if global_simlog_level >= 1:
        message_copy = deepcopy(message)
        message.replace('_dot_', '.')
        print('\033[91m' + 'WARNING: ' + message_copy + '\033[0m', file=sys.stderr)


