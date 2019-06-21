using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace OptimizerSet
{
    public class ZDT_Set
    {
        public double[] ZDT1(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            double h = 0.0;
            delta = 0.0;
            Fun[0] = X[0];
            g = 1.0 + 9.0 * (X.Sum() - X[0]) / (X_Dim - 1.0);
            h = 1.0 - Math.Sqrt(Fun[0] / g);
            Fun[1] = g * h;
            return Fun;
        }

        public double[] ZDT2(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            double h = 0.0;
            delta = 0.0;
            Fun[0] = X[0];
            g = 1.0 + 9.0 * (X.Sum() - X[0]) / (X_Dim - 1.0);
            h = 1.0 - Math.Pow(Fun[0] / g, 2.0);
            Fun[1] = g * h;
            return Fun;
        }

        public double[] ZDT3(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            double h = 0.0;
            delta = 0.0;
            Fun[0] = X[0];
            g = 1.0 + 9.0 * (X.Sum() - X[0]) / (X_Dim - 1.0);
            h = 1.0 - Math.Sqrt(Fun[0] / g) - (Fun[0] / g) * Math.Sin(10.0 * Math.PI * Fun[0]);
            Fun[1] = g * h;
            return Fun;
        }

        public double[] ZDT4(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            double h = 0.0;
            delta = 0.0;
            Fun[0] = X[0];
            g = 1.0 + 10.0 * (X_Dim - 1.0);
            for (int i = 1; i < X_Dim; i++)
            {
                g = g + (X[i] * X[i] - 10.0 * Math.Cos(4.0 * Math.PI * X[i]));
            }
            h = 1.0 - Math.Sqrt(Fun[0] / g);
            Fun[1] = g * h;
            return Fun;
        }

        public double[] ZDT6(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            double h = 0.0;
            delta = 0.0;
            Fun[0] = 1.0 - Math.Exp(-4.0 * X[0]) * Math.Pow(Math.Sin(6.0 * Math.PI * X[0]), 6.0);
            g = 1.0 + 9.0 * Math.Pow((X.Sum() - X[0]) / (X_Dim - 1.0), 0.25);
            h = 1.0 - Math.Pow(Fun[0] / g, 2.0);
            Fun[1] = g * h;
            return Fun;
        }
    }
}
