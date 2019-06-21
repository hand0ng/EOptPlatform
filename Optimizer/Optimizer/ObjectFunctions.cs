using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace OptimizerSet
{
    public class ObjectFunctions
    {
        public double ObjectFunctionWithCon(double[] X, out double delta)
        {
            //double Fun = MathNet.Numerics.TestFunctions.Ackley(X);
            ////double Fun = MathNet.Numerics.TestFunctions.Rastrigin(X);
            //delta = 0.0;
            //return Fun;

            return Graham(X, out delta);
        }

        public double[] MultiObjectFunction(int Obj_Dim, double[] X, out double delta)
        {
            DTLZ_Set DTLZ = new DTLZ_Set();
            ZDT_Set ZDT = new ZDT_Set();
            MOEAD_DE_Set F = new MOEAD_DE_Set();
            WFG_Set WFG = new WFG_Set();
            //return DTLZ.DTLZ6(Obj_Dim, X, out delta);
            //return ZDT.ZDT4(Obj_Dim, X, out delta);
            //return F.F6(Obj_Dim, X, out delta);
            return WFG.WFG1(Obj_Dim, X, out delta);
        }

        private double Graham(double[] X, out double Delta)
        {
            int nv = X.Length / 2;
            double f = 0.0;
            double f1 = 0.0;
            double[] r = new double[X.Length / 2];
            double[] theata = new double[X.Length / 2];
            for (int i = 0; i < X.Length / 2; i++)
            {
            	r[i] = X[i];
                theata[i] = X[i + X.Length / 2];
            }
            for (int i = 0; i < X.Length / 2 - 1; i++)
            {
            	f1 = f1 + (r[i + 1] * r[i] * Math.Sin(theata[i + 1] - theata[i]));
            }
            f = 1.0 / 0.5 / f1;
            Delta = nonlcon(X);
            return f;
        }

        private double nonlcon(double[] X)
        {
            int nv = X.Length / 2;
            double tmp_c = 0.0;
            double delta = 0.0;
            double[] r = new double[X.Length / 2];
            double[] rsqure = new double[X.Length / 2];
            double[] theata = new double[X.Length / 2];
            for (int i = 0; i < X.Length / 2; i++)
            {
            	r[i] = X[i];
                rsqure[i] = r[i] * r[i];
                theata[i] = X[i + X.Length / 2];
            }
            for (int i = 0; i < nv - 1; i++)
            {
            	for (int j = i + 1; j < nv; j++)
            	{
            		tmp_c = rsqure[i] + rsqure[j] - 2 * r[i] * r[j] * Math.Cos(theata[i] - theata[j]) - 1;
                    if (tmp_c > 0.0)
                    {
                        delta = delta + tmp_c;
                    }
            	}
            }
            for (int i = 0; i < nv - 1; i++)
            {
            	tmp_c = theata[i] - theata[i + 1];
                if (tmp_c > 0.0)
                {
                    delta = delta + tmp_c;
                }
            }
            return delta;
        }
    }
}
