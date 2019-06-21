using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace OptimizerSet
{
    public class DTLZ_Set
    {
        public double[] DTLZ1(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + Math.Pow(X[i] - 0.5, 2.0) - Math.Cos(20.0 * Math.PI * (X[i] - 0.5));
            }
            g = 100.0 * (X_Dim - Obj_Dim + 1.0 + g);
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = 0.5 * (1.0 + g);
                for (int j = 0; j < Obj_Dim - 1 - i; j++)
                {
                    Fun[i] = Fun[i] * X[j];
                }
                if (i > 0)
                {
                    Fun[i] = Fun[i] * (1.0 - X[Obj_Dim - 1 - i]);
                }
            }
            delta = 0.0;
            return Fun;
        }

        public double[] DTLZ2(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + Math.Pow(X[i] - 0.5, 2.0);
            }
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = (1.0 + g);
                for (int j = 0; j < Obj_Dim - 1 - i; j++)
                {
                    Fun[i] = Fun[i] * Math.Cos(X[j] * Math.PI / 2.0);
                }
                if (i > 0)
                {
                    Fun[i] = Fun[i] * Math.Sin(X[Obj_Dim - 1 - i] * Math.PI / 2.0);
                }
            }
            delta = 0.0;
            return Fun;
        }

        public double[] DTLZ3(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + Math.Pow(X[i] - 0.5, 2.0) - Math.Cos(20.0 * Math.PI * (X[i] - 0.5));
            }
            g = 100.0 * (X_Dim - Obj_Dim + 1.0 + g);
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = (1.0 + g);
                for (int j = 0; j < Obj_Dim - 1 - i; j++)
                {
                    Fun[i] = Fun[i] * Math.Cos(X[j] * Math.PI / 2.0);
                }
                if (i > 0)
                {
                    Fun[i] = Fun[i] * Math.Sin(X[Obj_Dim - 1 - i] * Math.PI / 2.0);
                }
            }
            delta = 0.0;
            return Fun;
        }

        public double[] DTLZ4(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + Math.Pow(X[i] - 0.5, 2.0);
            }
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = (1.0 + g);
                for (int j = 0; j < Obj_Dim - 1 - i; j++)
                {
                    Fun[i] = Fun[i] * Math.Cos(Math.Pow(X[j], 100.0) * Math.PI / 2.0);
                }
                if (i > 0)
                {
                    Fun[i] = Fun[i] * Math.Sin(Math.Pow(X[Obj_Dim - 1 - i], 100.0) * Math.PI / 2.0);
                }
            }
            delta = 0.0;
            return Fun;
        }

        public double[] DTLZ5(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + Math.Pow(X[i] - 0.5, 2.0);
            }
            for (int i = 0; i < Obj_Dim - 2; i++)
            {
                X[i + 1] = (1.0 + 2.0 * g * X[i + 1]) / (2.0 * (1.0 + g));
            }
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = (1.0 + g);
                for (int j = 0; j < Obj_Dim - 1 - i; j++)
                {
                    double theta = Math.PI / (4.0 * (1.0 + g)) * (1.0 + 2.0 * g * X[j]);
                    Fun[i] = Fun[i] * Math.Cos(X[j] * Math.PI / 2.0);
                }
                if (i > 0)
                {
                    double theta = Math.PI / (4.0 * (1.0 + g)) * (1.0 + 2.0 * g * X[Obj_Dim - 1 - i]);
                    Fun[i] = Fun[i] * Math.Sin(X[Obj_Dim - 1 - i] * Math.PI / 2.0);
                }
            }
            delta = 0.0;
            return Fun;
        }

        public double[] DTLZ6(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + Math.Pow(X[i], 0.1);
            }
            for (int i = 0; i < Obj_Dim - 2; i++)
            {
                X[i + 1] = (1.0 + 2.0 * g * X[i + 1]) / (2.0 * (1.0 + g));
            }
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = (1.0 + g);
                for (int j = 0; j < Obj_Dim - 1 - i; j++)
                {
                    double theta = Math.PI / (4.0 * (1.0 + g)) * (1.0 + 2.0 * g * X[j]);
                    Fun[i] = Fun[i] * Math.Cos(X[j] * Math.PI / 2.0);
                }
                if (i > 0)
                {
                    double theta = Math.PI / (4.0 * (1.0 + g)) * (1.0 + 2.0 * g * X[Obj_Dim - 1 - i]);
                    Fun[i] = Fun[i] * Math.Sin(X[Obj_Dim - 1 - i] * Math.PI / 2.0);
                }
            }
            delta = 0.0;
            return Fun;
        }

        public double[] DTLZ7(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            double g = 0.0;
            double h = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                g = g + X[i];
            }
            g = 1.0 + 9.0 * g / (X_Dim - Obj_Dim + 1.0);
            for (int i = 0; i < Obj_Dim - 1; i++)
            {
                Fun[i] = X[i];
                h = h + Fun[i] / (1.0 + g) * (1.0 + Math.Sin(3.0 * Math.PI * Fun[i]));
            }
            Fun[Obj_Dim - 1] = (1.0 + g) * (Obj_Dim - h);
            delta = 0.0;
            return Fun;
        }
    }
}
