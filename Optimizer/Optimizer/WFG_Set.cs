using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace OptimizerSet
{
    public class WFG_Set
    {
        public double[] WFG1(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            int K = Obj_Dim - 1;
            int L = X_Dim - K;
            double[] S = new double[Obj_Dim];
            double[] A = new double[Obj_Dim - 1];
            double[] z01 = new double[X_Dim];
            double[] t1 = new double[X_Dim];
            double[] t2 = new double[X_Dim];
            double[] t3 = new double[X_Dim];
            double[] t4 = new double[Obj_Dim];
            double[] x = new double[Obj_Dim];
            double[] h = new double[Obj_Dim];
            for (int i = 0; i < S.Length; i++)
            {
                S[i] = 2.0 + 2.0 * i;
            }
            for (int i = 0; i < A.Length; i++)
            {
                A[i] = 1.0;
            }
            for (int i = 0; i < X_Dim; i++)
            {
                z01[i] = X[i] / (2.0 + 2.0 * i);
                if (i < K)
                {
                    t1[i] = z01[i];
                    t2[i] = t1[i];
                    t3[i] = Math.Pow(t2[i], 0.02);
                }
                else
                {
                    t1[i] = s_linear(z01[i], 0.35);
                    t2[i] = b_flat(t1[i], 0.8, 0.75, 0.85);
                    t3[i] = Math.Pow(t2[i], 0.02);
                }
            }
            for (int i = 0; i < Obj_Dim - 1; i++)
            {
                t4[i] = t3[i];
            }
            double tmp = 0.0;
            double sum = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                tmp = tmp + t3[i] * 2.0 * (i + 1.0);
                sum = sum + 2.0 * (i + 1.0);
            }
            t4[Obj_Dim - 1] = tmp / sum;
            for (int i = 0; i < Obj_Dim - 1; i++)
            {
                x[i] = Math.Max(t4[Obj_Dim - 1], A[i]) * (t4[i] - 0.5) + 0.5;
            }
            x[Obj_Dim - 1] = t4[Obj_Dim - 1];
            h = convex(x);
            h[Obj_Dim - 1] = mixed(x);
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = x[Obj_Dim - 1] + S[i] * h[i];
            }
            delta = 0.0;
            return Fun;
        }

        public double[] WFG2(int Obj_Dim, double[] X, out double delta)
        {
            double[] Fun = new double[Obj_Dim];
            int X_Dim = X.Length;
            int K = Obj_Dim - 1;
            int L = X_Dim - K;
            double[] S = new double[Obj_Dim];
            double[] A = new double[Obj_Dim - 1];
            double[] z01 = new double[X_Dim];
            double[] t1 = new double[X_Dim];
            double[] t2 = new double[K + L / 2];
            double[] t3 = new double[X_Dim];
            double[] t4 = new double[Obj_Dim];
            double[] x = new double[Obj_Dim];
            double[] h = new double[Obj_Dim];
            for (int i = 0; i < S.Length; i++)
            {
                S[i] = 2.0 + 2.0 * i;
            }
            for (int i = 0; i < A.Length; i++)
            {
                A[i] = 1.0;
            }
            for (int i = 0; i < X_Dim; i++)
            {
                z01[i] = X[i] / (2.0 + 2.0 * i);
                if (i < K)
                {
                    t1[i] = z01[i];
                    t3[i] = Math.Pow(t2[i], 0.02);
                }
                else
                {
                    t1[i] = s_linear(z01[i], 0.35);
                    t3[i] = Math.Pow(t2[i], 0.02);
                }
            }
            for (int i = 0; i < K + L / 2; i++)
            {
                if (i < K)
                {
                    t2[i] = t1[i];
                }
                else
                {
                    double[] y = new double[K + 2 * (i - K) - (K + 2 * (i - K) - 1)];
                    t2[i] = b_flat(t1[i], 0.8, 0.75, 0.85);
                }
            }
            for (int i = 0; i < Obj_Dim - 1; i++)
            {
                t4[i] = t3[i];
            }
            double tmp = 0.0;
            double sum = 0.0;
            for (int i = Obj_Dim - 1; i < X_Dim; i++)
            {
                tmp = tmp + t3[i] * 2.0 * (i + 1.0);
                sum = sum + 2.0 * (i + 1.0);
            }
            t4[Obj_Dim - 1] = tmp / sum;
            for (int i = 0; i < Obj_Dim - 1; i++)
            {
                x[i] = Math.Max(t4[Obj_Dim - 1], A[i]) * (t4[i] - 0.5) + 0.5;
            }
            x[Obj_Dim - 1] = t4[Obj_Dim - 1];
            h = convex(x);
            h[Obj_Dim - 1] = mixed(x);
            for (int i = 0; i < Obj_Dim; i++)
            {
                Fun[i] = x[Obj_Dim - 1] + S[i] * h[i];
            }
            delta = 0.0;
            return Fun;
        }

        private double s_linear(double y, double A)
        {
            return Math.Abs(y - A) / Math.Abs(Math.Floor(A - y) + A);
        }

        private double b_flat(double y, double A, double B, double C)
        {
            return A + Math.Min(0.0, Math.Floor(y - B)) * A * (B - y) / B - Math.Min(0.0, Math.Floor(C - y)) * (1.0 - A) * (y - C) / (1.0 - C);
        }

        private double[] convex(double[] x)
        {
            int M = x.Length;
            double[] h = new double[M];
            for (int i = 0; i < M - 1; i++)
            {
                h[i] = 1.0;
                for (int j = 0; j < M - 1 - i; j++)
                {
                    h[i] = h[i] * (1.0 - Math.Cos(x[j] * Math.PI / 2.0));
                }
                if (i > 0)
                {
                    h[i] = h[i] * (1.0 - Math.Sin(x[M -1 - i] * Math.PI / 2.0));
                }
            }
            h[M - 1] = 1 - Math.Sin(x[0] * Math.PI / 2.0);
            return h;
        }

        private double mixed(double[] x)
        {
            return 1.0 - x[0] - Math.Cos(10.0 * Math.PI * x[0] + Math.PI / 2.0) / 10.0 / Math.PI;
        }

        private double r_nonsep(double[] y, int A)
        {
            double Output = 0.0;
            double tmp = 0.0;
            for (int j = 0; j < y.Length; j++)
            {
                tmp = 0.0;
            	for (int k = 0; k <= A - 2; k++)
            	{
                    tmp = tmp + Math.Abs(y[j]-y[(j+k)%y.Length]);
            	}
                Output = Output + y[j] + tmp;
            }
            Output = Output / (y.Length / A) / Math.Ceiling(A / 2.0) / (1.0 + 2.0 * A - 2.0 * Math.Ceiling(A / 2.0));
            return Output;
        }
    }
}
