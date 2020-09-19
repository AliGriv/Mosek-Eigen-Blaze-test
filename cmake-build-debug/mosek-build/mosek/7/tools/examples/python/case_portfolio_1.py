"""
  File : case_portfolio_1.py

  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  Description :  Implements a basic portfolio optimization model.
"""

import mosek

try:
    from numpy import zeros
except ImportError:
    from mosek.array import zeros

def streamprinter(text):
    print("%s" % text),
    
if __name__ == '__main__':    

  n     = 3
  gamma = 0.05
  mu    = [0.1073,  0.0737,  0.0627]
  GT    = [[0.1667,  0.0232,  0.0013],
           [0.0000,  0.1033, -0.0022],
           [0.0000,  0.0000,  0.0338]]
  x0    = [0.0, 0.0, 0.0]
  w     = 1.0

  inf   = 0.0 # This value has no significance 

  with mosek.Env() as env:
      with env.Task(0,0) as task:
          task.set_Stream(mosek.streamtype.log,streamprinter)

          rtemp = w
          for j in range(0,n):
              rtemp += x0[j]

          # Constraints. 
          task.appendcons(1+n)                 
          task.putconbound(0,mosek.boundkey.fx,rtemp,rtemp)
          task.putconname(0,"budget")

          task.putconboundlist(range(1+0,1+n),n*[mosek.boundkey.fx],n*[0.0],n*[0.0])
          for j in range(1,1+n) :
              task.putconname(j,"GT[%d]" % j)

          # Variables.
          task.appendvars(1+2*n)

          # Offset of variables into the API variable.      
          offsetx = 0   
          offsets = n   
          offsett = n+1 

          # x variables. 
          task.putclist(range(offsetx+0,offsetx+n),mu)
          task.putaijlist(n*[0],range(offsetx+0,offsetx+n),n*[1.0])
          for j in range(0,n):
              task.putaijlist(n*[1+j],range(offsetx+0,offsetx+n),GT[j])

          task.putvarboundlist(range(offsetx+0,offsetx+n),n*[mosek.boundkey.lo],n*[0.0],n*[inf])
          for j in range(0,n):
              task.putvarname(offsetx+j,"x[%d]" % (1+j))

          # s variable.
          task.putvarbound(offsets+0,mosek.boundkey.fx,gamma,gamma)
          task.putvarname(offsets+0,"s")

          # t variables. 
          task.putaijlist(range(1,n+1),range(offsett+0,offsett+n),n*[-1.0])
          task.putvarboundlist(range(offsett+0,offsett+n),n*[mosek.boundkey.fr],n*[-inf],n*[inf])
          for j in range(0,n):
              task.putvarname(offsett+j,"t[%d]" % (1+j))
     
          task.appendcone(mosek.conetype.quad,0.0,[offsets] + range(offsett,offsett+n))
          task.putconename(0,"stddev")

          task.putobjsense(mosek.objsense.maximize)
          
          # Turn all log output off. 
          task.putintparam(mosek.iparam.log,1)

          # Dump the problem to a human readable OPF file.  
          #task.writedata("dump.opf")
        
          task.optimize()

          # Display the solution summary for quick inspection of results. 
          task.solutionsummary(mosek.streamtype.msg)

          expret = 0.0
          x      = zeros(n,float)
          task.getxxslice(mosek.soltype.itr,offsetx+0,offsetx+n,x)
          for j in range(0,n):
              expret += mu[j]*x[j]

          stddev = zeros(1,float)     
          task.getxxslice(mosek.soltype.itr,offsets+0,offsets+1,stddev)

          print("\nExpected return %e for gamma %e\n" % (expret,stddev[0])) 	

