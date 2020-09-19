#
# Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File:    simple.py
#
# Purpose: Demonstrates a very simple example using MOSEK by
# reading a problem file, solving the problem and
# writing the solution to a file.
#

import mosek
import sys

def streamprinter(msg):
    sys.stdout.write (msg)
    sys.stdout.flush ()

if len(sys.argv) <= 1:
    print ("Missing argument, syntax is:")
    print ("  simple inputfile [ solutionfile ]")
else:
    # Create the mosek environment. 
    env  = mosek.Env ()

    # Create a task object linked with the environment env.
    # We create it with 0 variables and 0 constraints initially, 
    # since we do not know the size of the problem.
    task = env.Task (0, 0)
    task.set_Stream (mosek.streamtype.log, streamprinter)

    # We assume that a problem file was given as the first command
    # line argument (received in `argv')
    task.readdata (sys.argv[1])

    # Solve the problem
    task.optimize ()

    # Print a summary of the solution
    task.solutionsummary (mosek.streamtype.log)

    # If an output file was specified, write a solution
    if len(sys.argv) >= 3:
        # We define the output format to be OPF, and tell MOSEK to
        # leave out parameters and problem data from the output file.
        task.putintparam (mosek.iparam.write_data_format,    mosek.dataformat.op)
        task.putintparam (mosek.iparam.opf_write_solutions,  mosek.onoffkey.on)
        task.putintparam (mosek.iparam.opf_write_hints,      mosek.onoffkey.off)
        task.putintparam (mosek.iparam.opf_write_parameters, mosek.onoffkey.off)
        task.putintparam (mosek.iparam.opf_write_problem,    mosek.onoffkey.off)

        task.writedata (sys.argv[2])

