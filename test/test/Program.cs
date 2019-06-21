using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Accord.Math;
using MathNet.Numerics;
//using MathNet.Numerics.Statistics;
using MathNet.Numerics.LinearAlgebra;

namespace ConsoleApplication1
{
    class Program
    {
        private static Random rd = new Random(); // 随机种子
        private static int N_Pop = 100; // 种群规模
        private static int X_Dim = 0; // 变量维度
        private static int Obj_Dim = 0; // 目标维度
        private static double ProC = 1.0;
        private static double ProM = 1.0;
        private static double DisC = 20.0;
        private static double DisM = 20.0;
        private static double[,] Z; // 参考点
        private static double[] Zmin;
        private static double[,] X; // 种群
        private static double[,] X_Offspring; // 子代种群
        private static double[,] X_R; // 合并种群
        private static double[] X_Init; // 优化变量初始值
        private static double[] X_Lb; // 优化变量下限
        private static double[] X_Ub; // 优化变量上限
        private static double[,] Fun; // 各个个体的适应度
        private static double[,] Fun_Offspring; // 子代个体的适应度
        private static double[,] Fun_R; // 合并个体的适应度
        private static double[] Delta; // 各个个体的约束违反度
        private static double[] Delta_Offspring; // 各个个体的约束违反度
        private static double[] Delta_R; // 合并个体的约束违反度
        private static int[] Index; // 个体编号
        private static int[] FrontValue; // 非支配排序层编号
        private static int MaxFront = 0;
        private static OptimizerSet.General G = new OptimizerSet.General();



        static void Main(string[] args)
        {

            //double[] f = new double[] { 4300.000000, 3900.000000, 5000.000000, 6000.000000, 0.000000, 0.000000, 100000.000000 };
            //double[,] A = new double[,] { {0.935000,0.885000,1.030000,1.080000,0.940000,0.000000,-1.000000},
            //                              {0.300000,0.020000,0.009000,0.000000,0.010000,0.000000,-1.000000},
            //                            {0.200000,0.024000,0.300000,0.000000,0.000000,0.000000,-1.000000},
            //                            {0.005500,0.000000,0.010000,0.000000,0.000000,0.000000,-1.000000},
            //                            {0.000000,0.000000,0.000000,0.183000,0.000000,0.000000,-1.000000},
            //                            {0.430000,0.400000,0.830000,0.300000,0.000000,0.000000,-1.000000},
            //                            {0.820000,0.600000,1.150000,0.500000,0.000000,0.000000,-1.000000},
            //                            {1.700000,1.500000,1.540000,1.400000,0.000000,0.000000,-1.000000},
            //                            {2.070000,1.800000,1.970000,1.700000,0.000000,0.000000,-1.000000},
            //                            {7.360000,6.680000,8.000000,7.270000,0.000000,0.000000,-1.000000},
            //                            {0.180000,0.000000,0.000000,0.000000,0.000000,0.000000,-1.000000},
            //                            {0.680000,0.825000,0.440000,0.440000,0.000000,0.000000,-1.000000},
            //                            {-0.935000,-0.885000,-1.030000,-1.080000,-0.940000,-0.000000,-1.000000},
            //                            {-0.300000,-0.020000,-0.009000,-0.000000,-0.010000,-0.000000,-1.000000},
            //                            {-0.200000,-0.024000,-0.300000,-0.000000,-0.000000,-0.000000,-1.000000},
            //                            {-0.005500,-0.000000,-0.010000,-0.000000,-0.000000,-0.000000,-1.000000},
            //                            {-0.000000,-0.000000,-0.000000,-0.183000,-0.000000,-0.000000,-1.000000},
            //                            {-0.430000,-0.400000,-0.830000,-0.300000,-0.000000,-0.000000,-1.000000},
            //                            {-0.820000,-0.600000,-1.150000,-0.500000,-0.000000,-0.000000,-1.000000},
            //                            {-1.700000,-1.500000,-1.540000,-1.400000,-0.000000,-0.000000,-1.000000},
            //                            {-2.070000,-1.800000,-1.970000,-1.700000,-0.000000,-0.000000,-1.000000},
            //                            {-7.360000,-6.680000,-8.000000,-7.270000,-0.000000,-0.000000,-1.000000},
            //                            {-0.180000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-1.000000},
            //                            {-0.680000,-0.825000,-0.440000,-0.440000,-0.000000,-0.000000,-1.000000} };
            //double[] b = new double[] { -1.815042, 5.180390, 19.855153, 0.400111, 0.505404, 18.210585, 33.352646, 27.837326, 5.054039, 23.505292, -1.163231, 26.192758, 2.817827, 21.560557, 24.713092, 0.714095, 2.502953, 59.783844, 100.352089, 183.861838, 223.358217, 37.776045, 12.305292, 18.375487 };
            //double[,] Aeq = new double[,] { {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,0.000000},
            //                                {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
            //                                {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
            //                                {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
            //                                {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
            //                                {0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,0.000000},
            //                                {0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000} };
            //double[] beq = new double[] { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 };
            //double[] UB = new double[] { 4.000000, 4.000000, 4.000000, 4.000000, 0.000000, 0.000000, 10000.000000 };
            //double[] LB = new double[] { -4.000000, 0.000000, -4.000000, -4.000000, 0.000000, 0.000000, 0.000000 };
            //double Fun = 0.0;
            //int Flag = 0;
            //int[] i = new int[1];
            //OptimizerSet.Linprog lin = new OptimizerSet.Linprog();
            //double[] X = lin.StartOpt(f, LB, UB, out Fun, out Flag, A, b, Aeq, beq, null, null);
            //double f1 = Fun + 1;

            //savedata();

            test1 t1 = new test1();
            test2 t2 = new test2();
            test3 t3 = new test3();
            t3.test();

            //ServiceController scSQL = new ServiceController();
            //scSQL.MachineName = "localhost";
            //scSQL.ServiceName = "MSSQLSERVER";
            //scSQL.Start();
        }

        private static double[,] ReferencePoint(int N, int M, out int N_New)
        {
            int N1 = 0, N2 = 0, H1 = 1;
            while (MathNet.Numerics.Combinatorics.Combinations(H1 + M, M - 1) <= N)
            {
                H1 = H1 + 1;
            }
            N1 = (int)MathNet.Numerics.Combinatorics.Combinations(H1 + M - 1, M - 1);
            double[,] Z = new double[N1, M];
            double[,] Z1 = new double[N1, M];
            double[] H1_Vector = new double[H1 + M - 1];
            for (int i = 0; i < H1_Vector.Length; i++)
            {
                H1_Vector[i] = i + 1;
            }
            List<double[]> W_tmp1 = H1_Vector.Combinations(M - 1).ToList();
            for (int i = 0; i < N1; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    if (j == 0)
                    {
                        //Z[i, j] = (W_tmp1[i][j] - 1 - j) / H1;
                        Z1[i, j] = (W_tmp1[i][j] - 1) / H1;
                    }
                    else if (j == M - 1)
                    {
                        //Z[i, j] = (H1 - (W_tmp1[i][j - 1] - j)) / H1;
                        Z1[i, j] = (H1 - W_tmp1[i][j - 1] + j) / H1;
                    }
                    else
                    {
                        //Z[i, j] = (W_tmp1[i][j] - 1 - j - (W_tmp1[i][j - 1] - j)) / H1;
                        Z1[i, j] = (W_tmp1[i][j] - W_tmp1[i][j - 1] - 1) / H1;
                    }
                }
                for (int j = 0; j < N1; j++)
                {
                    for (int k = 0; k < M; k++)
                    {
                        Z[j, k] = Z1[j, k];
                    }
                }
            }
            if (H1 < M)
            {
                int H2 = 0;
                while ((MathNet.Numerics.Combinatorics.Combinations(H1 + M - 1, M - 1) + MathNet.Numerics.Combinatorics.Combinations(H2 + M, M - 1)) <= N)
                {
                    H2 = H2 + 1;
                }
                if (H2 > 0)
                {
                    N2 = (int)MathNet.Numerics.Combinatorics.Combinations(H2 + M - 1, M - 1);
                    double[] H2_Vector = new double[H2 + M - 1];
                    for (int i = 0; i < H2_Vector.Length; i++)
                    {
                        H2_Vector[i] = i + 1;
                    }
                    List<double[]> W_tmp2 = H2_Vector.Combinations(M - 1).ToList();
                    double[,] Z2 = new double[N2, M];
                    for (int i = 0; i < N2; i++)
                    {
                        for (int j = 0; j < M; j++)
                        {
                            if (j == 0)
                            {
                                //Z2[i, j] = (W_tmp2[i][j] - 1 - j) / H2;
                                Z2[i, j] = (W_tmp2[i][j] - 1) / H2;
                            }
                            else if (j == M - 1)
                            {
                                //Z2[i, j] = (H1 - (W_tmp2[i][j - 1] - j)) / H2;
                                Z2[i, j] = (H1 - W_tmp2[i][j - 1] + j) / H2;
                            }
                            else
                            {
                                //Z2[i, j] = (W_tmp2[i][j] - 1 - j - (W_tmp2[i][j - 1] - j)) / H2;
                                Z2[i, j] = (W_tmp2[i][j] - W_tmp2[i][j - 1] - 1) / H2;
                            }
                        }
                    }
                    Z = new double[N1 + N2, M];
                    for (int i = 0; i < N1; i++)
                    {
                        for (int j = 0; j < M; j++)
                        {
                            Z[i, j] = Z1[i, j];
                        }
                    }
                    for (int i = 0; i < N2; i++)
                    {
                        for (int j = 0; j < M; j++)
                        {
                            Z[i + N1, j] = Z2[i, j];
                        }
                    }
                }
            }
            N_New = N1 + N2;
            return Z;
        }

        private static double[] UpdateZmin(double[] Z, double[,] Fun, double[] Delta)
        {
            for (int i = 0; i < Z.Length; i++)
            {
                for (int j = 0; j < Fun.GetLength(0); j++)
                {
                    if (Z[i] > Fun[j, i] && Delta[j] == 0.0)
                    {
                        Z[i] = Fun[j, i];
                    }
                }
            }
            return Z;
        }
        /// <summary>
        /// 非支配排序
        /// An Efficient Approach to Nondominated Sorting for Evolutionary Multiobjective Optimization, TEVC, 2015
        /// </summary>
        private static void NDSort()
        {
            MaxFront = 0;
            int N = X_R.GetLength(0);
            int sum_s = 0;
            int[] Sorted = new int[N];
            BubbleSort();
            while (sum_s < N)
            {
                sum_s = 0;
                int[] ThisFront = new int[N];
                for (int i = 0; i < N; i++)
                {
                    if (Sorted[i] == 0)
                    {
                        bool NonDominate = true;
                        for (int j = 0; j < N; j++)
                        {
                            if (ThisFront[j] == 1)
                            {
                                NonDominate = false;
                                for (int j2 = 1; j2 < Obj_Dim; j2++)
                                {
                                    if ((Delta_R[i] == Delta_R[j] && Fun_R[i, j2] <= Fun_R[j, j2]) || Delta_R[i] < Delta_R[j])
                                    {
                                        NonDominate = true;
                                        break;
                                    }
                                }
                                if (!NonDominate)
                                {
                                    break;
                                }
                            }
                        }
                        if (NonDominate)
                        {
                            ThisFront[i] = 1;
                            Sorted[i] = 1;
                        }
                    }
                }
                for (int i = 0; i < N; i++)
                {
                    if (ThisFront[i] == 1)
                    {
                        FrontValue[i] = MaxFront;
                        ThisFront[i] = ThisFront[i] + 1;
                    }
                }
                MaxFront = MaxFront + 1;
                for (int i = 0; i < Sorted.Length; i++)
                {
                    sum_s = sum_s + Sorted[i];
                }
            }
        }

        /// <summary>
        /// 冒泡排序
        /// </summary>
        private static void BubbleSort()
        {
            bool flag;
            int N = Fun_R.GetLength(0);
            for (int i = N - 1; i > 0; i--)
            {
                flag = true;
                for (int j = 0; j < i; j++)
                {
                    flag = SortRow(j, j + 1, 0, flag);
                }
                if (flag) break;
            }
        }

        private static void BubbleSort1()
        {
            bool flag;
            int N = Z.GetLength(0);
            for (int i = N - 1; i > 0; i--)
            {
                flag = true;
                for (int j = 0; j < i; j++)
                {
                    flag = SortRow1(j, j + 1, 0, flag);
                }
                if (flag) break;
            }
        }

        /// <summary>
        /// 冒泡排序
        /// </summary>
        private static bool SortRow(int row1, int row2, int cols, bool flag)
        {
            double X_tmp = 0.0;
            double Fun_tmp = 0.0;
            double Delta_tmp = 0.0;
            if ((Delta_R[row1] == Delta_R[row2] && Fun_R[row1, cols] > Fun_R[row2, cols]) || Delta_R[row1] > Delta_R[row2])
            {
                for (int k = 0; k < X_Dim; k++)
                {
                    X_tmp = X_R[row1, k];
                    X_R[row1, k] = X_R[row2, k];
                    X_R[row2, k] = X_tmp;
                }
                for (int k = 0; k < Obj_Dim; k++)
                {
                    Fun_tmp = Fun_R[row1, k];
                    Fun_R[row1, k] = Fun_R[row2, k];
                    Fun_R[row2, k] = Fun_tmp;
                }
                Delta_tmp = Delta_R[row1];
                Delta_R[row1] = Delta_R[row2];
                Delta_R[row2] = Delta_tmp;
                flag = false;
            }
            else if ((Delta_R[row1] == Delta_R[row2]) && (Fun_R[row1, cols] == Fun_R[row2, cols]) && (cols < Obj_Dim - 1))
            {
                flag = SortRow(row1, row2, cols + 1, flag);
            }
            return flag;
        }

        private static bool SortRow1(int row1, int row2, int cols, bool flag)
        {
            double Z_tmp = 0.0;
            if (Z[row1, cols] > Z[row2, cols])
            {
                for (int k = 0; k < Obj_Dim; k++)
                {
                    Z_tmp = Z[row1, k];
                    Z[row1, k] = Z[row2, k];
                    Z[row2, k] = Z_tmp;
                }
                flag = false;
            }
            else if (Z[row1, cols] == Z[row2, cols])
            {
                flag = SortRow1(row1, row2, cols + 1, flag);
            }
            return flag;
        }

        /// <summary>
        /// 环境选择
        /// </summary>
        private static void EnvironmentalSelection()
        {
            bool Last_Stop = false;
            int Front_Last = 0;
            int N_LastPop = 0;
            int N_Pop_Next = 0;
            int[] N_Pop_EachFront = new int[MaxFront - 1];
            int[] FrontValue_tmp;
            double[,] X_tmp;
            double[,] Fun_tmp;
            double[] Delta_tmp;
            // Get the last front
            for (int i = 0; i < MaxFront - 1; i++)
            {
                for (int j = 0; j < 2 * N_Pop; j++)
                {
                    if (FrontValue[j] == i)
                    {
                        N_Pop_EachFront[i] = N_Pop_EachFront[i] + 1;
                    }
                }
                if (!Last_Stop)
                {
                    N_Pop_Next = N_Pop_Next + N_Pop_EachFront[i];
                    if (N_Pop_Next > N_Pop)
                    {
                        N_LastPop = N_Pop_EachFront[i];
                        Front_Last = i;
                        Last_Stop = true;
                    }
                }
            }
            FrontValue_tmp = new int[N_Pop_Next];
            X_tmp = new double[N_Pop_Next, X_Dim];
            Fun_tmp = new double[N_Pop_Next, Obj_Dim];
            Delta_tmp = new double[N_Pop_Next];
            N_Pop_Next = 0;
            for (int i = 0; i <= Front_Last; i++)
            {
                for (int j = 0; j < 2 * N_Pop; j++)
                {
                    if (FrontValue[j] == i)
                    {
                        for (int k = 0; k < X_Dim; k++)
                        {
                            X_tmp[N_Pop_Next, k] = X_R[j, k];
                        }
                        for (int k = 0; k < Obj_Dim; k++)
                        {
                            Fun_tmp[N_Pop_Next, k] = Fun_R[j, k];
                        }
                        Delta_tmp[N_Pop_Next] = Delta_R[j];
                        FrontValue_tmp[N_Pop_Next] = i;
                        N_Pop_Next = N_Pop_Next + 1;
                    }
                }
            }
            if ((N_Pop_Next - N_LastPop) < N_Pop)
            {
                var M = Matrix<double>.Build;
                bool Intercept_isNan = false;
                int K = N_Pop - (N_Pop_Next - N_LastPop);
                int[] ExtremeIndex = new int[Obj_Dim];
                int[] AssociateIndex = new int[N_Pop_Next];
                int[] AssociateCount = new int[N_Pop];
                int[] Pop_Last_Front_Index = new int[N_LastPop];
                double[] Hyperplane = new double[Obj_Dim];
                double[] Intercept = new double[Obj_Dim];
                double[,] w = new double[Obj_Dim, Obj_Dim];
                double[,] Fun_tmp_Norm = new double[N_Pop_Next, Obj_Dim];
                double[,] Distance_to_Z = new double[N_Pop_Next, N_Pop];
                double[,] X_K = new double[K, X_Dim];
                double[,] Fun_K = new double[K, Obj_Dim];
                double[] Delta_K = new double[K];
                Matrix<double> H = M.Dense(Obj_Dim, Obj_Dim);
                for (int i = 0; i < Obj_Dim; i++)
                {
                    for (int j = 0; j < Obj_Dim; j++)
                    {
                        w[i, j] = 1.0E-6;
                        if (i == j)
                        {
                            w[i, j] = w[i, j] + 1.0;
                        }
                    }
                }
                // Translate objectives
                for (int i = 0; i < N_Pop_Next; i++)
                {
                    for (int j = 0; j < Obj_Dim; j++)
                    {
                        Fun_tmp_Norm[i, j] = Fun_tmp[i, j] - Zmin[j];
                    }
                }



                // Detect the extreme points
                for (int i = 0; i < Obj_Dim; i++)
                {
                    double[] Fun_E = new double[N_Pop_Next];
                    for (int j = 0; j < N_Pop_Next; j++)
                    {
                        double[] Fun_E_tmp = new double[Obj_Dim];
                        for (int k = 0; k < Obj_Dim; k++)
                        {
                            Fun_E_tmp[k] = Fun_tmp_Norm[j, k] * w[i, k];
                        }
                        Fun_E[j] = Fun_E_tmp.Max();
                    }
                    ExtremeIndex[i] = G.minIndex(Fun_E);
                }
                // Calculate the intercepts of the hyperplane constructed by the extreme
                for (int i = 0; i < Obj_Dim; i++)
                {
                    for (int j = 0; j < Obj_Dim; j++)
                    {
                        H[i, j] = Fun_tmp_Norm[ExtremeIndex[i], j];
                    }
                }
                Hyperplane = (H.Inverse() * M.Dense(Obj_Dim, 1, 1)).ToColumnWiseArray();
                for (int i = 0; i < Obj_Dim; i++)
                {
                    if (Hyperplane[i] != 0.0)
                    {
                        Intercept[i] = 1 / Hyperplane[i];
                    }
                    else
                    {
                        Intercept_isNan = true;
                        break;
                    }
                }
                if (Intercept_isNan)
                {
                    for (int i = 0; i < Obj_Dim; i++)
                    {
                        Intercept[i] = M.DenseOfArray(Fun_tmp_Norm).Column(i).Max();
                    }
                }
                // Normalization
                for (int i = 0; i < N_Pop_Next; i++)
                {
                    for (int j = 0; j < Obj_Dim; j++)
                    {
                        Fun_tmp_Norm[i, j] = Fun_tmp_Norm[i, j] / Intercept[j];
                    }
                }

                // Associate each solution with one reference point
                // Calculate the distance of each solution to each reference vector
                // Associate each solution with its nearest reference point
                int Index_Count = 0;
                for (int i = 0; i < N_Pop_Next; i++)
                {
                    double[] cos = new double[N_Pop];
                    double[] distance = new double[N_Pop];
                    for (int j = 0; j < N_Pop; j++)
                    {
                        double[] Ftmp = new double[Obj_Dim];
                        double[] Ztmp = new double[Obj_Dim];
                        for (int k = 0; k < Obj_Dim; k++)
                        {
                            Ftmp[k] = Fun_tmp_Norm[i, k];
                            Ztmp[k] = Z[j, k];
                        }
                        cos[j] = 1 - MathNet.Numerics.Distance.Cosine(Ftmp, Ztmp);
                        distance[j] = Math.Sqrt(Ftmp.SquareEuclidean()) * Math.Sqrt(1 - cos[j] * cos[j]);
                        Distance_to_Z[i, j] = distance[j];
                    }
                    AssociateIndex[i] = G.minIndex(distance);
                    //AssociateCount[AssociateIndex[i]] = AssociateCount[AssociateIndex[i]] + 1;
                    // Calculate the number of associated solutions except for the last front of each reference point
                    if (FrontValue_tmp[i] != Front_Last)
                    {
                        AssociateCount[AssociateIndex[i]] = AssociateCount[AssociateIndex[i]] + 1;
                    }
                    else
                    {
                        Pop_Last_Front_Index[Index_Count] = i;
                        Index_Count++;
                    }
                }
                //List<string[]> comeout = new List<string[]>();
                //for (int k = 0; k < Distance_to_Z.GetLength(0); k++)
                //{
                //    string[] Fun1 = new string[1];
                //    for (int j = 0; j < Fun1.Length; j++)
                //    {
                //        Fun1[j] = Distance_to_Z[k, j].ToString();
                //    }
                //    comeout.Add(Fun1);
                //}
                //string comeoutPath = "C:\\Users\\Administrator\\Desktop\\222.csv";
                //G.WriteCSV(comeoutPath, false, comeout);

                //double[] min_D = new double[N_Pop_Next];
                //for (int i = 0; i < N_Pop_Next; i++)
                //{
                //    double[] Distance_to_Z_tmp = new double[N_Pop];
                //    for (int j = 0; j < N_Pop; j++)
                //    {
                //        Distance_to_Z_tmp[j] = Distance_to_Z[i, j];
                //    }
                //    min_D[i] = Distance_to_Z_tmp.Min();
                //}
                //double min_D1 = min_D[0];

                //List<string[]> comeout = new List<string[]>();
                //for (int k = 0; k < AssociateIndex.GetLength(0); k++)
                //{
                //    string[] Fun1 = new string[1];
                //    Fun1[0] = AssociateIndex[k].ToString();
                //    comeout.Add(Fun1);
                //}
                //string comeoutPath = "C:\\Users\\Administrator\\Desktop\\222.csv";
                //G.WriteCSV(comeoutPath, false, comeout);

                //List<string[]> comeout = new List<string[]>();
                //for (int k = 0; k < Fun_tmp_Norm.GetLength(0); k++)
                //{
                //    string[] Fun1 = new string[3];
                //    for (int i = 0; i < Obj_Dim; i++)
                //    {
                //        Fun1[i] = Fun_tmp_Norm[k, i].ToString();
                //    }
                //    Fun1[2] = AssociateIndex[k].ToString();
                //    comeout.Add(Fun1);
                //}
                //string comeoutPath = "C:\\Users\\Administrator\\Desktop\\222.csv";
                //G.WriteCSV(comeoutPath, false, comeout);

                // Niching
                // Select K solutions one by one
                int Sel = 0;
                int[] Chosen = new int[N_Pop_Next];
                while (Sel < K)
                {
                    // Select the least crowded reference point
                    ArrayList J_Min = G.minIndexArray(AssociateCount);
                    ArrayList Candi_Index = new ArrayList();
                    ArrayList Candi_Distance = new ArrayList();
                    int J_Chose = 0;
                    if (J_Min.Count > 0)
                    {
                        J_Chose = G.rndInt_uni(0, J_Min.Count - 1, rd);
                    }
                    int J_Chose_Min = (int)J_Min[J_Chose];
                    for (int j = 0; j < N_LastPop; j++)
                    {
                        if (AssociateIndex[Pop_Last_Front_Index[j]] == J_Chose_Min && Chosen[Pop_Last_Front_Index[j]] == 0)
                        {
                            Candi_Index.Add(Pop_Last_Front_Index[j]);
                            Candi_Distance.Add(Distance_to_Z[Pop_Last_Front_Index[j], (int)J_Min[J_Chose]]);
                        }
                    }

                    if (Candi_Index.Count > 0)
                    {
                        int candi;
                        if (AssociateCount[(int)J_Min[J_Chose]] == 0)
                        {
                            double[] candi_d = (double[])Candi_Distance.ToArray(typeof(double));
                            candi = G.minIndex(candi_d);
                            for (int j = 0; j < X_Dim; j++)
                            {
                                X_K[Sel, j] = X_tmp[(int)Candi_Index[candi], j];
                            }
                            for (int j = 0; j < Obj_Dim; j++)
                            {
                                Fun_K[Sel, j] = Fun_tmp[(int)Candi_Index[candi], j];
                            }
                            Delta_K[Sel] = Delta_tmp[(int)Candi_Index[candi]];
                        }
                        else
                        {
                            candi = G.rndInt_uni(0, Candi_Index.Count - 1, rd);
                            for (int j = 0; j < X_Dim; j++)
                            {
                                X_K[Sel, j] = X_tmp[(int)Candi_Index[candi], j];
                            }
                            for (int j = 0; j < Obj_Dim; j++)
                            {
                                Fun_K[Sel, j] = Fun_tmp[(int)Candi_Index[candi], j];
                            }
                            Delta_K[Sel] = Delta_tmp[(int)Candi_Index[candi]];
                        }
                        Sel = Sel + 1;
                        Chosen[(int)Candi_Index[candi]] = 1;
                        AssociateCount[J_Chose] = AssociateCount[J_Chose] + 1;
                    }
                    else
                    {
                        // Exclude the reference point
                        AssociateCount[(int)J_Min[J_Chose]] = int.MaxValue;
                    }
                }
                for (int i = 0; i < N_Pop; i++)
                {
                    if (i < (N_Pop_Next - N_LastPop))
                    {
                        for (int j = 0; j < X_Dim; j++)
                        {
                            X[i, j] = X_tmp[i, j];
                        }
                        for (int j = 0; j < Obj_Dim; j++)
                        {
                            Fun[i, j] = Fun_tmp[i, j];
                        }
                        Delta[i] = Delta_tmp[i];
                    }
                    else
                    {
                        for (int j = 0; j < X_Dim; j++)
                        {
                            X[i, j] = X_K[i - (N_Pop_Next - N_LastPop), j];
                        }
                        for (int j = 0; j < Obj_Dim; j++)
                        {
                            Fun[i, j] = Fun_K[i - (N_Pop_Next - N_LastPop), j];
                        }
                        Delta[i] = Delta_K[i - (N_Pop_Next - N_LastPop)];
                    }
                }
            }
            else
            {
                for (int i = 0; i < N_Pop; i++)
                {
                    for (int j = 0; j < X_Dim; j++)
                    {
                        X[i, j] = X_tmp[i, j];
                    }
                    for (int j = 0; j < Obj_Dim; j++)
                    {
                        Fun[i, j] = Fun_tmp[i, j];
                    }
                    Delta[i] = Delta_tmp[i];
                }
            }
        }

        private static void Variation(int[] Index_Parent)
        {
            
            double beta;
            double mu;
            double k;
            double[] X_tmp = new double[X_Dim];
            double[] Fun_tmp = new double[Obj_Dim];
            double[,] X_Offspring_tmp = new double[Index_Parent.Length, X_Dim];
            double[,] p1 = new double[N_Pop/2, X_Dim];
            double[,] p2 = new double[N_Pop/2, X_Dim];
            for (int i = 0; i < N_Pop/2; i++)
            {
                for (int j = 0; j < X_Dim; j++)
            	{
                    p1[i, j] = X[Index_Parent[i], j];
                    p2[i, j] = X[Index_Parent[i + N_Pop/2], j];
            	}
            }
            // SBX crossover
            //for (int i = 0; i < Index_Parent.Length; i = i + 2)
            for (int i = 0; i < Index_Parent.Length/2; i++)
            {
                for (int j = 0; j < X_Dim; j++)
                {
                    //mu = G.rnd_uni(rd);
                    mu = 0.7;
                    if (mu <= 0.5)
                    {
                        beta = Math.Pow(2.0 * mu, (1.0 / (DisC + 1.0)));
                    }
                    else
                    {
                        beta = Math.Pow(2.0 - 2.0 * mu, (-1.0 / (DisC + 1.0)));
                    }
                    //beta = beta * Math.Pow(-1, G.rndInt_uni(0, 1, rd));
                    beta = beta * Math.Pow(-1, 1);
                    //if (G.rnd_uni(rd) > ProC)
                    //{
                    //    beta = 1.0;
                    //}
                    //X_Offspring_tmp[i, j] = (X[Index_Parent[i], j] + X[Index_Parent[i + 1], j]) / 2 + beta * (X[Index_Parent[i], j] - X[Index_Parent[i + 1], j]) / 2;
                    //X_Offspring_tmp[i + 1, j] = (X[Index_Parent[i], j] + X[Index_Parent[i + 1], j]) / 2 - beta * (X[Index_Parent[i], j] - X[Index_Parent[i + 1], j]) / 2;
                    X_Offspring_tmp[i, j] = (X[Index_Parent[i], j] + X[Index_Parent[i + N_Pop / 2], j]) / 2 + beta * (X[Index_Parent[i], j] - X[Index_Parent[i + N_Pop / 2], j]) / 2;
                    X_Offspring_tmp[i + Index_Parent.Length / 2, j] = (X[Index_Parent[i], j] + X[Index_Parent[i + N_Pop / 2], j]) / 2 - beta * (X[Index_Parent[i], j] - X[Index_Parent[i + N_Pop / 2], j]) / 2;
                }
            }
            
            // 多项式变异
            for (int i = 0; i < X_Offspring.GetLength(0); i++)
            {
                for (int j = 0; j < X_Dim; j++)
                {
                    k = G.rnd_uni(rd);
                    //k = 1;
                    mu = G.rnd_uni(rd);
                    //mu = 0.4;
                    if (k <= ProM && mu < 0.5)
                    {
                        X_Offspring_tmp[i, j] = X_Offspring_tmp[i, j] + (X_Ub[j] - X_Lb[j]) * (Math.Pow(2.0 * mu + (1.0 - 2.0 * mu) * Math.Pow(1 - (X_Offspring_tmp[i, j] - X_Lb[j]) / (X_Ub[j] - X_Lb[j]), DisM + 1.0), 1.0 / (DisM + 1.0)) - 1.0);
                        X_Offspring[i, j] = G.BoundConstraint_2(X_Offspring_tmp[i, j], X_Lb[j], X_Ub[j]);

                    }
                    else if (k <= ProM && mu >= 0.5)
                    {
                        X_Offspring_tmp[i, j] = X_Offspring_tmp[i, j] + (X_Ub[j] - X_Lb[j]) * (1 - Math.Pow(2.0 * (1.0 - mu) + 2.0 * (mu - 0.5) * Math.Pow(1 - (X_Ub[j] - X_Offspring_tmp[i, j]) / (X_Ub[j] - X_Lb[j]), DisM + 1.0), 1.0 / (DisM + 1.0)));
                        X_Offspring[i, j] = G.BoundConstraint_2(X_Offspring_tmp[i, j], X_Lb[j], X_Ub[j]);
                    }
                    else
                    {
                        X_Offspring[i, j] = G.BoundConstraint_2(X_Offspring_tmp[i, j], X_Lb[j], X_Ub[j]);
                    }
                    X_tmp[j] = X_Offspring[i, j];
                }
                //Fun_tmp = ObjectFunction(X_tmp, out Delta_Offspring[i]); // 得到种群适应度和约束违反度
                //for (int j = 0; j < Obj_Dim; j++)
                //{
                //    Fun_Offspring[i, j] = Fun_tmp[j];
                //}
            }
        }
    }
}