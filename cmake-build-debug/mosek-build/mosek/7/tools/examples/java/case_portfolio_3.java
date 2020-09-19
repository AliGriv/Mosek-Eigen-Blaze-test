package com.mosek.example;

/*
  File : case_portfolio_3.java

  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  Description :  Implements a basic portfolio optimization model.
*/

public class case_portfolio_3
{
    static final int n = 3;

    public static void main (String[] args)
    {
        // Since the value infinity is never used, we define
        // 'infinity' symbolic purposes only
        double infinity = 0;
        double gamma = 0.05;
        double[]   mu = {0.1073,  0.0737,  0.0627};
        double[][] GT={
            {0.1667,  0.0232,  0.0013},
            {0.0000,  0.1033, -0.0022},
            {0.0000,  0.0000,  0.0338}
        };
        double[] x0 = {0.0, 0.0, 0.0};
        double   w = 1.0;
        double[] m={0.01, 0.01, 0.01}; 

        int offsetx = 0;    
        int offsets = offsetx + n;    
        int offsett = offsets + 1;
        int offsetc = offsett + n;
        int offsetv = offsetc + n;
        int offsetz = offsetv + n;
        int offsetf = offsetz + n;
        int offsetg = offsetf + 3*n;

        int numvar = offsetg + 3*n;

        int offset_con_budget = 0;
        int offset_con_gx_t = offset_con_budget+1;
        int offset_con_abs1 = offset_con_gx_t+n;
        int offset_con_abs2 = offset_con_abs1+n;
        int offset_con_f=  offset_con_abs2+n; 
        int offset_con_g=  offset_con_f + 3*n; 

        int numcon = 1 + 3*n + 2*3*n;

        mosek.Env env = null;
        mosek.Task task = null;
        
        try
            {
                // Make mosek environment. 
                env  = new mosek.Env ();
                // Create a task object linked with the environment env.
                task = new mosek.Task (env, 0, 0);

                // Directs the log task stream to the user specified
                // method task_msg_obj.stream
                task.set_Stream(
                                mosek.Env.streamtype.log,
                                new mosek.Stream() 
                                { public void stream(String msg) { System.out.print(msg); }});

                //Set up constraint bounds, names and variable coefficients
                task.appendcons(numcon);
                for( int i=0; i<n; ++i)
                    {
                        w+=x0[i];
                        task.putconbound(offset_con_gx_t+i, mosek.Env.boundkey.fx, 0.,0.);
                        task.putconname(offset_con_gx_t+i,"GT["+(i+1)+"]");

                        task.putconbound(offset_con_abs1+i, mosek.Env.boundkey.lo, -x0[i],infinity);
                        task.putconname(offset_con_abs1+i,"zabs1["+(i+1)+"]");

                        task.putconbound(offset_con_abs2+i, mosek.Env.boundkey.lo, x0[i],infinity);
                        task.putconname(offset_con_abs2+i,"zabs2["+(i+1)+"]");

                        for(int j=0;j<3;++j)
                            {
                                task.putconbound(offset_con_f+3*i+j, mosek.Env.boundkey.fx, 0.,0.);
                                task.putconname(offset_con_f+3*i+j,"f["+(i+1)+","+(j+1)+"]");

                                task.putconbound(offset_con_g+3*i+j, mosek.Env.boundkey.fx, 0.,0.);
                                task.putconname(offset_con_g+3*i+j,"g["+(i+1)+","+(j+1)+"]");
                            }

                        task.putconbound(offset_con_g+3*i+1, mosek.Env.boundkey.fx, -1./8,-1./8.);

                    }
                // e x = w + e x0
                task.putconbound(offset_con_budget, mosek.Env.boundkey.fx,w,w);
                task.putconname(offset_con_budget,"budget");

                //Variables.
                task.appendvars(numvar);
    
                //the objective function coefficients
                int[] xindx={offsetx+0,offsetx+1,offsetx+2};
                task.putclist(xindx,mu);

                double[] one_m_one={1.0,-1.0};
                double[] one_one={1.0,1.0};

                //set up variable bounds and names
                for( int i=0; i<n; ++i)
                    {
                        task.putvarbound(offsetx+i, mosek.Env.boundkey.lo,0.,infinity);
                        task.putvarbound(offsett+i, mosek.Env.boundkey.fr,infinity,infinity);
                        task.putvarbound(offsetc+i, mosek.Env.boundkey.fr,infinity,infinity);
                        task.putvarbound(offsetz+i, mosek.Env.boundkey.fr,infinity,infinity);
                        task.putvarbound(offsetv+i, mosek.Env.boundkey.fr,infinity,infinity);
                        for(int j=0;j<3;++j)
                            {
                                task.putvarbound(offsetf+j+i*3, mosek.Env.boundkey.fr,infinity,infinity);
                                task.putvarbound(offsetg+j+i*3, mosek.Env.boundkey.fr,infinity,infinity);
                            }
                        task.putvarname(offsetx+i,"x["+(i+1)+"]");
                        task.putvarname(offsett+i,"t["+(i+1)+"]");
                        task.putvarname(offsetc+i,"c["+(i+1)+"]");
                        task.putvarname(offsetz+i,"z["+(i+1)+"]");
                        task.putvarname(offsetv+i,"v["+(i+1)+"]");
                        for(int j=0;j<3;++j)
                            {
                                task.putvarname(offsetf+j+i*3, "f["+(i+1)+","+(j+1)+"]");
                                task.putvarname(offsetg+j+i*3, "g["+(i+1)+","+(j+1)+"]");
                            }

                        for( int j=i; j<n;++j)
                            task.putaij(offset_con_gx_t+i,j, GT[i][j]); 

                        task.putaij(offset_con_gx_t+i, offsett+i,-1.0);

                        task.putaij(offset_con_budget,offsetx+i, 1.0); 
                        task.putaij(offset_con_budget,offsetc+i, m[i]); 

                        // z_j - x_j >= -x0_j
                        int[] indx1={offsetz+i,offsetx+i};
                        task.putarow(offset_con_abs1+i, indx1, one_m_one);
                        // z_j + x_j >= +x0_j
                        int[] indx2={offsetz+i, offsetx+i};
                        task.putarow(offset_con_abs2+i,indx2,one_one);

                        int[] indxf1={ offsetv+i, offsetf+i*3};
                        task.putarow(offset_con_f+3*i, indxf1,one_m_one);
                        int[] indxf2={offsetc+i, offsetf+i*3+1};
                        task.putarow(offset_con_f+1+3*i,indxf2,one_m_one);
                        int[] indxf3={offsetz+i, offsetf+i*3+2};
                        task.putarow(offset_con_f+2+3*i,indxf3,one_m_one);

                        int[] indxg1={offsetz+i, offsetg+i*3};
                        task.putarow(offset_con_g+3*i,indxg1,one_m_one);

                        task.putaij(offset_con_g+3*i+1,offsetg+i*3+1,-1.);

                        int[] indxg3={offsetv+i, offsetg+i*3+2};
                        task.putarow(offset_con_g+3*i+2,indxg3,one_m_one);
                    }
                task.putvarbound(offsets, mosek.Env.boundkey.fx,gamma,gamma);
                task.putvarname(offsets,"s");

                //Cones.
                int conecount=0;

                int[] csub = {offsets,offsett+0,offsett+1,offsett+2};
                task.appendcone(mosek.Env.conetype.quad, 0.0,  csub);
                task.putconename(conecount,"stddev");
                ++conecount;

                for(int j=0; j<n; ++j,++conecount)
                    {
                        int[] coneindx={offsetf + j*3 , offsetf +j*3 +1,  offsetf +j*3 +2};
                        task.appendcone(mosek.Env.conetype.rquad, 0.0, coneindx);
                        task.putconename(conecount,"f["+(j+1)+"]");
                    }

                for(int j=0; j<n; ++j, ++conecount)
                    {
                        int[] coneindx={offsetg + j*3 , offsetg +j*3 +1,  offsetg +j*3 +2};
                        task.appendcone(mosek.Env.conetype.rquad, 0.0, coneindx);
                        task.putconename(conecount,"g["+(j+1)+"]");
                    }
                /* A maximization problem */ 
                task.putobjsense(mosek.Env.objsense.maximize);

                /* Solve the problem */
                try
                    {
                        //Turn all log output off.  
                        //task.putintparam(mosek.Env.iparam.log,0);  

                        //task.writedata("dump.opf");

                        task.optimize();

                        task.solutionsummary(mosek.Env.streamtype.log);

                        double expret=0.0,stddev=0.0;
                        double[] xx = new double[numvar];

                        task.getxx(mosek.Env.soltype.itr,xx);

                        for(int j=0; j<n; ++j)  
                            expret += mu[j]*xx[j+offsetx];

                        System.out.printf("eglected return %e for gamma %e\n\n",expret, xx[offsets]); 
    
                    }
                catch (mosek.Warning mw)
                    {
                        System.out.println (" Mosek warning:");
                        System.out.println (mw.toString ());
                    }
            }
        catch( mosek.Exception e)
            {
                System.out.println ("An error/warning was encountered"); 
                System.out.println (e.toString ());
                throw e;
            }
        finally 
            { 
                if (task != null) task.dispose (); 
                if (env  != null)  env.dispose (); 
            } 
    }
}