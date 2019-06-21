using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Accord.Math;
using InterfaceDeclare;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace OptimizerSet
{
    /// <summary>
    /// An Evolutionary Many-Objective Optimization Algorithm Using Reference-point Based Non-dominated Sorting Approach, Part I: Solving Problems with Box Constraints, TEVC, 2014
    /// </summary>
    public class NSGAIII
    {
        #region 私有变量
        private Random rd; // 随机种子
        //private Random rd1; // 随机种子
        //private Random rd2; // 随机种子
        //private Random rd3; // 随机种子
        //private Random rd4; // 随机种子
        //private Random rd5; // 随机种子
        private int Max_FEs = 30000; // 最大评价次数
        private int N_Pop = 100; // 种群规模
        private int X_Dim = 0; // 变量维度
        private int Obj_Dim = 0; // 目标维度
        private double ProC = 1.0;
        private double ProM = 1.0;
        private double DisC = 30.0;
        private double DisM = 20.0;  
        private double[,] Z; // 参考点
        private double[] Zmin;
        private double[] Zmax;
        private double[,] X; // 种群
        private double[,] X_Offspring; // 子代种群
        private double[,] X_R; // 合并种群
        private double[] X_Init; // 优化变量初始值
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[,] Fun; // 各个个体的适应度
        private double[,] Fun_Offspring; // 子代个体的适应度
        private double[,] Fun_R; // 合并个体的适应度
        private double[] Delta; // 各个个体的约束违反度
        private double[] Delta_Offspring; // 各个个体的约束违反度
        private double[] Delta_R; // 合并个体的约束违反度
        private int[] FrontValue; // 非支配排序层编号
        private int MaxFront = 0;
        private General G;
        private IModel[] model;//模型接口
        private IMultiFitness fitness;//适应度函数接口
        #endregion

        /// <summary>
        /// 设置适应度函数
        /// </summary>
        /// <param name="fitval">适应度函数接口</param>
        /// <param name="modelval">模型接口</param>
        public void SetFitness(IMultiFitness fitval, IModel[] modelval)
        {
            fitness = fitval;
            model = modelval;
        }

        /// <summary>
        /// 开始优化
        /// </summary>
        /// <param name="ConfigFile">算法参数配置文件路径</param>
        /// <param name="X_Init">优化变量初值</param>
        /// <param name="X_MaxStep">优化变量最大步长</param>
        /// <param name="X_Lb">优化变量下限</param>
        /// <param name="X_Ub">优化变量上限</param>
        /// <returns>使适应度函数最小的变量值</returns>
        /// <returns>最优目标函数</returns>
        public double[,] StartOpt(string ConfigFile, double[] X_Init, double[] X_MaxStep, double[] X_Lb, double[] X_Ub, out double[,] Fun)
        {
            if (X_Init.Length != X_MaxStep.Length || X_Init.Length != X_Lb.Length || X_Init.Length != X_Ub.Length)
            {
                throw new Exception("Variable number are not set correctly!");
            }
            this.G = new General();
            List<string[]> Config = new List<string[]>();
            try
            {
                Config = G.ReadCSV(ConfigFile);
                this.Max_FEs = int.Parse(Config[0][1]); // 最大评价次数
                this.N_Pop = int.Parse(Config[1][1]); // 种群规模
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.Obj_Dim = fitness.GetObjectiveNumber();
            //this.Obj_Dim = 3;
            this.Z = ReferencePoint(N_Pop, Obj_Dim, out N_Pop); // 生成参考点
            this.X_Dim = X_Init.Length;
            this.ProM = 1.0 / X_Dim;
            this.Zmin = new double[Obj_Dim];
            this.Zmax = new double[Obj_Dim];
            this.X = new double[N_Pop, X_Dim];
            this.X_Offspring = new double[N_Pop, X_Dim];
            this.X_R = new double[N_Pop * 2, X_Dim];
            this.X_Init = new double[X_Dim];
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.Fun = new double[N_Pop, Obj_Dim]; // 各个个体的适应度
            this.Fun_Offspring = new double[N_Pop, Obj_Dim];
            this.Fun_R = new double[N_Pop * 2, Obj_Dim];
            this.Delta = new double[N_Pop]; // 各个个体的约束违反度
            this.Delta_Offspring = new double[N_Pop];
            this.Delta_R = new double[N_Pop * 2];
            this.FrontValue = new int[N_Pop * 2];

            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Init[i] = X_Init[i];
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
            }
            for (int i = 0; i < Obj_Dim; i++)
            {
                Zmin[i] = double.MaxValue;
            }
            this.rd = new Random();
            //this.rd1 = new Random();
            //this.rd2 = new Random();
            //this.rd3 = new Random();
            //this.rd4 = new Random();
            //this.rd5 = new Random();
            Opt(); // 优化
            Fun = new double[N_Pop, Obj_Dim];
            for (int i = 0; i < N_Pop; i++)
            {
                for (int j = 0; j < Obj_Dim; j++)
            	{
                    Fun[i, j] = this.Fun[i, j];
            	}
            }
            return this.X;
        }

        /// <summary>
        /// NSGAIII
        /// </summary>
        private void Opt()
        {
            int FEs = 0; // 当前评价次数
            int[] Index_Parent; // 父代个体编号
            // 上下限内初始化种群
            Initial();

            //更新Zmin
            Zmin = UpdateZmin(Zmin, Fun, Delta);
            Zmax = UpdateZmax(Fun, Delta);

            // 开始寻优
            // 种群停止进化，结束优化流程
            while (FEs < Max_FEs)
            {
                if (N_Pop % 2 != 0)
                {
                    // 二元锦标赛选择N_Pop + 1个父代个体编号
                    Index_Parent = BinaryTournamentSelection(N_Pop + 1);
                }
                else
                {
                    // 二元锦标赛选择N_Pop个父代个体编号
                    Index_Parent = BinaryTournamentSelection(N_Pop);
                }
                
                // SBX与多项式变异生成子代
                Variation(Index_Parent);

                //更新Zmin
                Zmin = UpdateZmin(Zmin, Fun_Offspring, Delta_Offspring);

                // 合并父代子代
                for (int i = 0; i < N_Pop; i++)
                {
                    for (int j = 0; j < X_Dim; j++)
                    {
                        X_R[i, j] = X[i, j];
                        X_R[i + N_Pop, j] = X_Offspring[i, j];
                    }
                    for (int j = 0; j < Obj_Dim; j++)
                    {
                        Fun_R[i, j] = Fun[i, j];
                        Fun_R[i + N_Pop, j] = Fun_Offspring[i, j];
                    }
                    Delta_R[i] = Delta[i];
                    Delta_R[i + N_Pop] = Delta_Offspring[i];
                }

                // 非支配排序
                NDSort();

                // 环境选择
                EnvironmentalSelection();
                Zmax = UpdateZmax(Fun, Delta);
                //if (FEs%10500 == 0)
                //{
                //    //rd = new Random();
                //    System.Threading.Thread.Sleep(10);
                //}
                FEs = FEs + N_Pop;
                //Console.WriteLine("NSGAIII " + FEs);
            }
        }

        /// <summary>
        /// 生成参考点
        /// </summary>
        /// <param name="N">个体数</param>
        /// <param name="M">目标数</param>
        /// <param name="N_New">新个体数</param>
        /// <returns>参考点</returns>
        private double[,] ReferencePoint(int N, int M, out int N_New)
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
            for (int i = 0; i < N1 + N2; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    if (Z[i, j] <= 1e-6)
                    {
                        Z[i, j] = 1e-6;
                    }
                }
            }
            N_New = N1 + N2;
            return Z;
        }

        /// <summary>
        /// 二元锦标赛选择
        /// </summary>
        /// <param name="N">选择个体数</param>
        /// <returns>选择结果</returns>
        private int[] BinaryTournamentSelection(int N)
        {
            int[] p = new int[N];
            for (int i = 0; i < N_Pop; i++)
            {
                double[] Fun_p1 = new double[Obj_Dim];
                double[] Fun_p2 = new double[Obj_Dim];
                double Delta_p1 = 0.0;
                double Delta_p2 = 0.0;
                int p1 = G.rndInt_uni(0, N_Pop - 1, rd);
                int p2 = G.rndInt_uni(0, N_Pop - 1, rd);
                while (p2 == p1) p2 = G.rndInt_uni(0, N_Pop - 1, rd);
                for (int j = 0; j < Obj_Dim; j++)
                {
                    Fun_p1[j] = Fun[p1, j];
                    Fun_p2[j] = Fun[p2, j];
                }
                Delta_p1 = Delta[p1];
                Delta_p2 = Delta[p2];
                //if (isDominate(Fun_p1, Fun_p2, Delta_p1, Delta_p2))
                if ((Delta_p1 == Delta_p2 && Fun_p1.Sum() < Fun_p2.Sum()) || (Delta_p1 < Delta_p2))
            	{
                    p[i] = p1;
            	}
                //else if (isDominate(Fun_p2, Fun_p1, Delta_p2, Delta_p1))
                else if ((Delta_p1 == Delta_p2 && Fun_p1.Sum() > Fun_p2.Sum()) || (Delta_p1 > Delta_p2))
                {
                    p[i] = p2;
                }
                else
                {
                    if (G.rnd_uni(rd) > 0.5)
                    {
                        p[i] = p1;
                    }
                    else
                    {
                        p[i] = p2;
                    }
                }
            }
            if (p.Length > N_Pop)
            {
                p[N_Pop] = p[0];
            }
            return p;
        }

        /// <summary>
        /// a是否支配b
        /// </summary>
        /// <param name="a">个体a的目标函数值</param>
        /// <param name="b">个体b的目标函数值</param>
        /// <returns>支配结果</returns>
        private bool isDominate(double[] f1, double[] f2, double d1, double d2)
        {
            int flag = 0;
            bool dominate = true;
            if (d1 > d2)
            {
                dominate = false;
            }
            else if (d1 == d2)
            {
                for (int i = 0; i < f1.Length; i++)
                {
                    if (f1[i] > f2[i])
                    {
                        dominate = false;
                        break;
                    }
                    if (f1[i] == f2[i])
                    {
                        flag = flag + 1;
                    }
                }
                if (flag == f1.Length)
                {
                    dominate = false;
                }
            }
            return dominate;
        }

        /// <summary>
        /// SBX与多项式变异
        /// </summary>
        /// <param name="Index_Parent">父代个体编号</param>
        private void Variation(int[] Index_Parent)
        {
            object tmp;
            double x = 0.0;
            tmp = x;
            double beta;
            double mu;
            double k;
            double[] X_tmp = new double[X_Dim];
            double[] Fun_tmp = new double[Obj_Dim];
            double[,] X_Offspring_tmp = new double[Index_Parent.Length, X_Dim];
            // SBX crossover
            for (int i = 0; i < Index_Parent.Length / 2; i++)
            {
            	for (int j = 0; j < X_Dim; j++)
            	{
                    mu = G.rnd_uni(rd);
                    if (mu <= 0.5)
                    {
                        beta = Math.Pow(2.0 * mu, (1.0 / (DisC + 1.0)));
                    }
                    else
                    {
                        beta = Math.Pow(2.0 - 2.0 * mu, (-1.0 / (DisC + 1.0)));
                    }
                    beta = beta * Math.Pow(-1, G.rndInt_uni(0, 1, rd));
                    if (G.rnd_uni(rd) > ProC)
                    {
                        beta = 1.0;
                    }
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
                    mu = G.rnd_uni(rd);
                    if (k <= ProM && mu < 0.5)
                    {
                        X_Offspring_tmp[i, j] = X_Offspring_tmp[i, j] + (X_Ub[j] - X_Lb[j]) * (Math.Pow(2.0 * mu + (1.0 - 2.0 * mu) * Math.Pow(1 - (X_Offspring_tmp[i, j] - X_Lb[j]) / (X_Ub[j] - X_Lb[j]), DisM + 1.0), 1.0 / (DisM + 1.0)) - 1.0);
                        X_Offspring[i, j] = G.BoundConstraint(X_Offspring_tmp[i, j], X_Lb[j], X_Ub[j]);

                    }
                    else if (k <= ProM && mu >= 0.5)
                    {
                        X_Offspring_tmp[i, j] = X_Offspring_tmp[i, j] + (X_Ub[j] - X_Lb[j]) * (1 - Math.Pow(2.0 * (1.0 - mu) + 2.0 * (mu - 0.5) * Math.Pow(1 - (X_Ub[j] - X_Offspring_tmp[i, j]) / (X_Ub[j] - X_Lb[j]), DisM + 1.0), 1.0 / (DisM + 1.0)));
                        X_Offspring[i, j] = G.BoundConstraint(X_Offspring_tmp[i, j], X_Lb[j], X_Ub[j]);
                    }
                    else
                    {
                        X_Offspring[i, j] = G.BoundConstraint(X_Offspring_tmp[i, j], X_Lb[j], X_Ub[j]);
                    }
                    X_tmp[j] = X_Offspring[i, j];
                }
                Fun_tmp = ObjectFunction(X_tmp, out Delta_Offspring[i]); // 得到种群适应度和约束违反度
                for (int j = 0; j < Obj_Dim; j++)
                {
                    Fun_Offspring[i, j] = Fun_tmp[j];
                }
            }
        }

        /// <summary>
        /// 非支配排序
        /// An Efficient Approach to Nondominated Sorting for Evolutionary Multiobjective Optimization, TEVC, 2015
        /// </summary>
        private void NDSort()
        {
            MaxFront = 0;
            int N = X_R.GetLength(0);
            int[] Sorted = new int[N];
            BubbleSort();
            while (Sorted.Sum() < N)
            {
                MaxFront = MaxFront + 1;
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
            }
        }

        /// <summary>
        /// 冒泡排序
        /// </summary>
        private void BubbleSort()
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

        /// <summary>
        /// 替换
        /// </summary>
        private bool SortRow(int row1, int row2, int cols, bool flag)
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

        /// <summary>
        /// 环境选择
        /// </summary>
        private void EnvironmentalSelection()
        {
            bool Last_Stop = false;
            int Front_Last = 0;
            int N_LastPop = 0;
            int N_Pop_Next = 0;
            int[] N_Pop_EachFront = new int[MaxFront];
            int[] FrontValue_tmp;
            double[,] X_tmp; 
            double[,] Fun_tmp;
            double[] Delta_tmp;
            // Get the last front
            for (int i = 1; i <= MaxFront; i++)
            {
                for (int j = 0; j < 2 * N_Pop; j++)
                {
                    if (FrontValue[j] == i)
                    {
                        N_Pop_EachFront[i - 1] = N_Pop_EachFront[i - 1] + 1;
                    } 
                }
                if (!Last_Stop)
                {
                    N_Pop_Next = N_Pop_Next + N_Pop_EachFront[i - 1];
                    if (N_Pop_Next > N_Pop)
                    {
                        N_LastPop = N_Pop_EachFront[i - 1];
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
            for (int i = 1; i <= Front_Last; i++)
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
                var V = Vector<double>.Build;
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
                double[,] Distance_to_Z = new double[N_Pop_Next,N_Pop];
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
                            Fun_E_tmp[k] = Fun_tmp_Norm[j, k] / w[i, k];
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
                        //distance[j] = V.DenseOfArray(Ftmp).L2Norm() * Math.Sqrt(1 - cos[j] * cos[j]);
                        Distance_to_Z[i, j] = distance[j];
                	}
                    AssociateIndex[i] = G.minIndex(distance);
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
                    for (int j = 0; j < N_LastPop; j++)
                    {
                        if (AssociateIndex[Pop_Last_Front_Index[j]] == (int)J_Min[J_Chose] && Chosen[Pop_Last_Front_Index[j]] == 0)
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
                        AssociateCount[(int)J_Min[J_Chose]] = AssociateCount[(int)J_Min[J_Chose]] + 1;
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

        /// <summary>
        /// 随机初始化种群
        /// </summary>
        private void Initial()
        {
            double[] X_tmp = new double[X_Dim];
            double[] Fun_tmp = new double[Obj_Dim];
            for (int i = 0; i < N_Pop; i++)
            {
                for (int k = 0; k < X_Dim; k++)
                {
                    if (i == 0)
                    {
                        X[i, k] = X_Init[k];
                    }
                    else
                    {
                        X[i, k] = G.rnd_uni(rd) * (this.X_Ub[k] - this.X_Lb[k]) + X_Lb[k];
                    }
                    X_tmp[k] = X[i, k];
                }
                Fun_tmp = ObjectFunction(X_tmp, out Delta[i]); // 得到种群适应度和约束违反度
                for (int k = 0; k < Obj_Dim; k++)
                {
                    Fun[i, k] = Fun_tmp[k];
                }
            }
        }

        /// <summary>
        /// 更新Zmin
        /// </summary>
        /// <param name="Z">Zmin</param>
        /// <param name="Fun">目标函数值</param>
        /// <returns>Zmin</returns>
        private double[] UpdateZmin(double[] Z, double[,] Fun, double[] Delta)
        {
            for (int i = 0; i < Z.Length; i++)
            {
                for (int j = 0; j < Fun.GetLength(0); j++)
                {
                    if (Z[i] > Fun[j,i] && Delta[j] == 0.0)
                    {
                        Z[i] = Fun[j, i];
                    }
                }
            }
            return Z;
        }

        private double[] UpdateZmax(double[,] Fun, double[] Delta)
        {
            double[] Z = new double[Fun.GetLength(1)];
            for (int i = 0; i < Z.Length; i++)
            {
                for (int j = 0; j < Fun.GetLength(0); j++)
                {
                    if (Z[i] < Fun[j, i] && Delta[j] == 0.0)
                    {
                        Z[i] = Fun[j, i];
                    }
                }
            }
            return Z;
        }

        /// <summary>
        /// 目标函数
        /// </summary>
        /// <param name="X">目标函数输入——优化变量</param>
        /// <param name="model">模型接口</param>
        /// <param name="Delta">违反度，输出参数，未违反约束时为0</param>
        /// <returns>适应度函数值</returns>
        private double[] ObjectFunction(double[] X, out double Delta)
        {
            return fitness.Fitness(X, model, out Delta);
            //ObjectFunctions O = new ObjectFunctions();
            //return O.MultiObjectFunction(Obj_Dim, X, out Delta);
            //double[] fun;
            //fun = O.MultiObjectFunction(Obj_Dim, X, out Delta);
            //return fun;
        }
    }
}