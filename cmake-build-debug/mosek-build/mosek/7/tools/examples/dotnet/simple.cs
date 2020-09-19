/*
  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  File:    simple.cs

  Purpose: Demonstrates a very simple example using MOSEK by
  reading a problem file, solving the problem and
  writing the solution to a file.
*/

using System;

class msgclass : mosek.Stream 
{
  public override void streamCB (string msg)
  {
    Console.Write ("{0}", msg);
  }
}

public class simple
{  
  public static void Main (string[] args)
  {
    if (args.Length == 0)
    {
      Console.WriteLine ("Missing argument, syntax is:");
      Console.WriteLine ("  simple inputfile [ solutionfile ]");
    }
    else
    {
      using (mosek.Env env = new mosek.Env())
      { 
        using (mosek.Task task = new mosek.Task(env))
        {
          task.set_Stream (mosek.streamtype.log, new msgclass ());

          // We assume that a problem file was given as the first command
          // line argument (received in `args')
          task.readdata (args[0]);

          // Solve the problem
          task.optimize ();

          // Print a summary of the solution
          task.solutionsummary (mosek.streamtype.log);
          
          // If an output file was specified, write a solution
          if (args.Length >= 2)
          {
            // We define the output format to be OPF, and tell MOSEK to
            // leave out parameters and problem data from the output file.
            task.putintparam (mosek.iparam.write_data_format,    mosek.dataformat.op);
            task.putintparam (mosek.iparam.opf_write_solutions,  mosek.onoffkey.on);
            task.putintparam (mosek.iparam.opf_write_hints,      mosek.onoffkey.off);
            task.putintparam (mosek.iparam.opf_write_parameters, mosek.onoffkey.off);
            task.putintparam (mosek.iparam.opf_write_problem,    mosek.onoffkey.off);
            
            task.writedata(args[1]);
          }
        }
      }
    }
  }
}
