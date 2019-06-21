using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace OptimizerSet
{
    public class MOEAD_DE_Set
    {
        public double[] F2(int Obj_Dim, double[] X, out double delta)
        {
            double[] fun = new double[2];
            int m = X.Length;
            int count1 = 0;
            int count2 = 0;
            for (int i = 2; i <= m; i++)
            {
                if (i % 2 == 0)
                    count2++;
                else
                    count1++;
            }
            int[] J1 = new int[count1];
            int[] J2 = new int[count2];
            delta = 0.0;
            double[] X1_tmp = new double[count1];
            double[] X2_tmp = new double[count2];

            for (int i = 2, j = 0, k = 0; i <= m; i++)
            {
                if (i % 2 == 0)
                {
                    J2[j] = i;
                    j++;
                }
                else
                {
                    J1[k] = i;
                    k++;
                }
            }

            for (int i = 0; i < count1; i++)
                X1_tmp[i] = Math.Pow(X[J1[i] - 1] - Math.Sin(6.0 * Math.PI * X[0] + J1[i] * Math.PI / m), 2);
            for (int i = 0; i < count2; i++)
                X2_tmp[i] = Math.Pow(X[J2[i] - 1] - Math.Sin(6.0 * Math.PI * X[0] + J2[i] * Math.PI / m), 2);
            fun[0] = X[0] + (2.0 / (double)count1) * X1_tmp.Sum();
            fun[1] = 1 - Math.Sqrt(X[0]) + 2.0 / (double)count2 * X2_tmp.Sum();

            return fun;
        }

        public double[] F5(int Obj_Dim, double[] X, out double delta)
        {
            double[] fun = new double[2];
            int m = X.Length;
            delta = 0.0;
            int count1 = 0;
            int count2 = 0;
            for (int i = 2; i <= m; i++)
            {
                if (i % 2 == 0)
                    count2++;
                else
                    count1++;
            }
            int[] J1 = new int[count1];
            int[] J2 = new int[count2];

            for (int i = 2, j = 0, k = 0; i <= m; i++)
            {
                if (i % 2 == 0)
                {
                    J2[j] = i;
                    j++;
                }
                else
                {
                    J1[k] = i;
                    k++;
                }
            }

            double[] Y1_tmp = new double[count1];
            double[] Y2_tmp = new double[count2];
            for (int i = 0; i < count1; i++)
                Y1_tmp[i] = Math.Pow(X[J1[i] - 1] - (0.3 * Math.Pow(X[0], 2) * Math.Cos(24.0 * Math.PI * X[0] + 4.0 * J1[i] * Math.PI / m) + 0.6 * X[0]) * Math.Cos(6.0 * Math.PI * X[0] + J1[i] * Math.PI / m), 2);

            for (int i = 0; i < count2; i++)
                Y2_tmp[i] = Math.Pow(X[J2[i] - 1] - (0.3 * Math.Pow(X[0], 2) * Math.Cos(24.0 * Math.PI * X[0] + 4.0 * J2[i] * Math.PI / m) + 0.6 * X[0]) * Math.Sin(6.0 * Math.PI * X[0] + J2[i] * Math.PI / m), 2);

            fun[0] = X[0] + 2.0 / (double)count1 * Y1_tmp.Sum();
            fun[1] = 1 - Math.Sqrt(X[0]) + 2.0 / (double)count2 * Y2_tmp.Sum();
            return fun;
        }

        public double[] F8(int Obj_Dim, double[] X, out double delta)
        {
            double[] fun = new double[2];
            delta = 0.0;
            int m = X.Length;
            int count1 = 0;
            int count2 = 0;
            for (int i = 2; i <= m; i++)
            {
                if (i % 2 == 0)
                    count2++;
                else
                    count1++;
            }
            int[] J1 = new int[count1];
            int[] J2 = new int[count2];

            for (int i = 2, j = 0, k = 0; i <= m; i++)
            {
                if (i % 2 == 0)
                {
                    J2[j] = i;
                    j++;
                }
                else
                {
                    J1[k] = i;
                    k++;
                }
            }

            double[] Y1_tmp = new double[count1];
            double[] Y2_tmp = new double[count2];
            double[] Y1_tmp1 = new double[count1];
            double[] Y2_tmp1 = new double[count2];
            for (int i = 0; i < count1; i++)
            {
                Y1_tmp[i] = Math.Pow(X[J1[i] - 1] - Math.Pow(X[0], 0.5 * (1.0 + 3 * (J1[i] - 2) / (m - 2))), 2);
                Y1_tmp1[i] = Math.Cos(20.0 * J1[i] * Math.PI / Math.Sqrt(J1[i]));
            }
            for (int i = 0; i < count2; i++)
            {
                Y2_tmp[i] = Math.Pow(X[J2[i] - 1] - Math.Pow(X[0], 0.5 * (1.0 + 3 * (J2[i] - 2) / (m - 2))), 2);
                Y2_tmp1[i] = Math.Cos(20.0 * J2[i] * Math.PI / Math.Sqrt(J2[i]));
            }
            double y1 = 1.0;
            double y2 = 1.0;
            foreach (double tmp in Y1_tmp1)
                y1 = y1 * tmp;
            foreach (double tmp in Y2_tmp1)
                y2 = y2 * tmp;

            fun[0] = X[0] + 2.0 / (double)count1 * (4.0 * Y1_tmp.Sum() - 2.0 * y1 + 2);
            fun[1] = 1 - Math.Sqrt(X[0]) + 2.0 / (double)count2 * (4.0 * Y2_tmp.Sum() - 2.0 * y2 + 2);

            return fun;
        }

        public double[] F6(int Obj_Dim, double[] X, out double delta)
        {
            double[] fun = new double[3];
            delta = 0.0;
            int m = X.Length;
            int count1 = 0, count2 = 0, count3 = 0;
            for (int i = 3; i <= m; i++)
            {
                if ((i - 1) % 3 == 0)
                    count1++;
                else if ((i - 2) % 3 == 0)
                    count2++;
                else
                    count3++;
            }
            int[] J1 = new int[count1];
            int[] J2 = new int[count2];
            int[] J3 = new int[count3];

            for (int i = 3, j = 0, k = 0, p = 0; i <= m; i++)
            {
                if ((i - 1) % 3 == 0)
                {
                    J1[j] = i;
                    j++;
                }
                if ((i - 2) % 3 == 0)
                {
                    J2[k] = i;
                    k++;
                }
                if (i % 3 == 0)
                {
                    J3[p] = i;
                    p++;
                }
            }
            double[] F1_tmp = new double[count1];
            double[] F2_tmp = new double[count2];
            double[] F3_tmp = new double[count3];

            for (int i = 0; i < count1; i++)
                F1_tmp[i] = Math.Pow(X[J1[i] - 1] - 2.0 * X[1] * Math.Sin(2.0 * Math.PI * X[0] + J1[i] * Math.PI / (double)m), 2);
            for (int i = 0; i < count2; i++)
                F2_tmp[i] = Math.Pow(X[J2[i] - 1] - 2.0 * X[1] * Math.Sin(2.0 * Math.PI * X[0] + J2[i] * Math.PI / (double)m), 2);
            for (int i = 0; i < count3; i++)
                F3_tmp[i] = Math.Pow(X[J3[i] - 1] - 2.0 * X[1] * Math.Sin(2.0 * Math.PI * X[0] + J3[i] * Math.PI / (double)m), 2);

            fun[0] = Math.Cos(0.5 * X[0] * Math.PI) * Math.Cos(0.5 * X[1] * Math.PI) + 2.0 / (double)count1 * F1_tmp.Sum();
            fun[1] = Math.Cos(0.5 * X[0] * Math.PI) * Math.Sin(0.5 * X[1] * Math.PI) + 2.0 / (double)count2 * F2_tmp.Sum();
            fun[2] = Math.Sin(0.5 * X[0] * Math.PI) + 2.0 / (double)count3 * F3_tmp.Sum();

            return fun;
        }

        private double[] F9(int Obj_Dim, double[] X, out double delta)
        {
            double[] fun = new double[3];
            delta = 0.0;
            int m = X.Length;
            int count1 = 0, count2 = 0, count3 = 0;
            double kir = 0.1;

            for (int i = 3; i <= m; i++)
            {
                if ((i - 1) % 3 == 0)
                    count1++;
                else if ((i - 2) % 3 == 0)
                    count2++;
                else
                    count3++;
            }
            int[] J1 = new int[count1];
            int[] J2 = new int[count2];
            int[] J3 = new int[count3];

            for (int i = 3, j = 0, k = 0, p = 0; i <= m; i++)
            {
                if ((i - 1) % 3 == 0)
                {
                    J1[j] = i;
                    j++;
                }
                if ((i - 2) % 3 == 0)
                {
                    J2[k] = i;
                    k++;
                }
                if (i % 3 == 0)
                {
                    J3[p] = i;
                    p++;
                }
            }

            double[] F1_tmp = new double[count1];
            double[] F2_tmp = new double[count2];
            double[] F3_tmp = new double[count3];

            for (int i = 0; i < count1; i++)
                F1_tmp[i] = Math.Pow(X[J1[i] - 1] - 2.0 * X[1] * Math.Sin(2.0 * Math.PI * X[0] + J1[i] * Math.PI / (double)m), 2);
            for (int i = 0; i < count2; i++)
                F2_tmp[i] = Math.Pow(X[J2[i] - 1] - 2.0 * X[1] * Math.Sin(2.0 * Math.PI * X[0] + J2[i] * Math.PI / (double)m), 2);
            for (int i = 0; i < count3; i++)
                F3_tmp[i] = Math.Pow(X[J3[i] - 1] - 2.0 * X[1] * Math.Sin(2.0 * Math.PI * X[0] + J3[i] * Math.PI / (double)m), 2);

            fun[0] = 0.5 * (Math.Max(0.0, (1.0 + kir) * (1 - 4.0 * Math.Pow(2.0 * X[0] - 1, 2))) + 2.0 * X[0]) * X[1] + 2.0 / (double)count1 * F1_tmp.Sum();
            fun[1] = 0.5 * (Math.Max(0.0, (1.0 + kir) * (1 - 4.0 * Math.Pow(2.0 * X[0] - 1, 2))) - 2.0 * X[0] + 2.0) * X[1] + 2.0 / (double)count2 * F2_tmp.Sum();
            fun[2] = 1 - X[1] + 2.0 / (double)count3 * F3_tmp.Sum();

            return fun;
        }
    }
}
