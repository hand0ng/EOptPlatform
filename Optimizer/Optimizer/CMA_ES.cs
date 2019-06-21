using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using InterfaceDeclare;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace OptimizerSet
{
    /// <summary>
    /// Completely Derandomized Self-Adaptation in Evolution Strategies, Evolutionary Computation, 2001
    /// <summary>
    public class CMA_ES
    {
        #region 私有变量
        private Random rd; // 随机种子
        private int Max_FEs = 30000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private int X_Dim = 0; // 变量维度
        private double[] X;
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] X_Opt; // 最优解

        private Matrix<double> xmeanw; // Start point
        private Matrix<double> xold;
        private Matrix<double> zmeanw;
        private bool flginiphase = true; // Initial phase
        private double sigma = 0.25; // Initial step size
        private double Min_sigma = 1.0e-15; // Minimal step size
        private double Max_sigma = 0.0; // Maximum step size      
        // Parameter setting: selection
        private int lambda = 0;
        private int mu = 0;
        private Matrix<double> arweights;
        // Parameter setting: adaptation
        private double cc = 0.0;
        private double ccov = 0.0;
        private double cs = 0.0;
        private double damp = 0.0;
        // Initialize dynamic strategy parameters and constants
        private Matrix<double> B;
        private Matrix<double> D;
        private Matrix<double> BD;
        private Matrix<double> C;
        private Matrix<double> CC;
        private Matrix<double> pc;
        private Matrix<double> ps;
        private double cw = 0.0;
        private double chiN = 0.0;
        // Lambda offspring
        private Matrix<double> arz;
        private Matrix<double> arx;
        private Matrix<double> arz_mu;
        private Matrix<double> arx_mu;
        private double[] Fun; // 各个个体的适应度
        private double[] Delta; // 各个个体的约束违反度
        private int[] Index; // 各个个体的序号
        private double Best_Fun = 0.0; // 最优个体适应度
        private double Best_Delta = 0.0; // 最优个体约束违反度
        private double Last_Best_Fun = 0.0; // 上次最优个体适应度
        private double TolFun = 1; // 两代最优解之差
        private General G;
        private IModel[] model;//模型接口
        private IFitness fitness;//适应度函数接口
        #endregion

        /// <summary>
        /// 设置适应度函数
        /// </summary>
        /// <param name="fitval">适应度函数接口</param>
        /// <param name="modelval">模型接口</param>
        public void SetFitness(IFitness fitval, IModel[] modelval)
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
        public double[] StartOpt(string ConfigFile, double[] X_Init, double[] X_MaxStep, double[] X_Lb, double[] X_Ub, out double Fun, out double Delta)
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
                this.Min_TolFun = double.Parse(Config[1][1]); // 结束条件
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            var M = Matrix<double>.Build; // 生成矩阵
            //var V = Vector<double>.Build; // 生成向量
            this.X_Dim = X_Init.Length;
            this.X = new double[X_Dim];
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            double[] X_Diff = new double[X_Dim];
            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
                X_Diff[i] = this.X_Ub[i] - this.X_Lb[i];
            }
            this.X_Opt = new double[X_Dim];

            this.xmeanw = M.Dense(X_Dim, 1);
            this.xold = M.Dense(X_Dim, 1);
            this.zmeanw = M.Dense(X_Dim, 1);
            this.flginiphase = true;
            this.sigma = 0.25; // Initial step size
            this.Min_sigma = 1.0e-15; // Minimal step size
            this.Max_sigma = X_Diff.Max() / Math.Sqrt((double)X_Dim);
            for (int i = 0; i < X_Dim; i++)
            {
                this.xmeanw[i, 0] = X_Init[i];
            }

            this.lambda = 4 + (int)Math.Floor(3.0 * Math.Log((double)X_Dim));
            this.mu = (int)Math.Floor((double)lambda / 2.0);
            this.arweights = M.Dense(mu, 1);
            for (int i = 0; i < mu; i++)
            {
                this.arweights[i, 0] = Math.Log((lambda + 1.0) / 2.0) - Math.Log(i + 1.0);
            }

            this.cc = 4.0 / (X_Dim + 4.0);
            this.ccov = 2.0 / Math.Pow((X_Dim + Math.Pow(2.0, 0.5)), 2.0);
            this.cs = 4.0 / (X_Dim + 4.0);
            this.damp = (1.0 - Math.Min(0.7, X_Dim * lambda / Max_FEs)) / cs + 1.0;

            this.B = M.DenseIdentity(X_Dim);
            this.D = M.DenseIdentity(X_Dim);
            this.BD = M.Dense(X_Dim, X_Dim);
            this.C = M.Dense(X_Dim, X_Dim);
            this.CC = M.Dense(X_Dim, X_Dim);
            this.BD = B * D;
            this.C = BD * BD.Transpose();
            this.pc = M.Dense(X_Dim, 1);
            this.ps = M.Dense(X_Dim, 1);
            this.cw = arweights.RowSums().Sum() / arweights.L2Norm();
            this.chiN = Math.Pow(X_Dim, 0.5) * (1.0 - 1.0 / (4.0 * X_Dim) + 1.0 / (21.0 * Math.Pow(X_Dim, 2.0)));

            this.arz = M.Dense(X_Dim, lambda);
            this.arx = M.Dense(X_Dim, lambda);
            this.arz_mu = M.Dense(X_Dim, mu);
            this.arx_mu = M.Dense(X_Dim, mu);

            this.Fun = new double[lambda]; // 各个个体的适应度
            this.Delta = new double[lambda]; // 各个个体的约束违反度
            this.Index = new int[lambda]; // 各个个体的序号

            this.rd = new Random();
            Opt(); // 优化
            Fun = Best_Fun;
            Delta = Best_Delta;
            return X_Opt;
        }

        /// <summary>
        /// CMA_ES
        /// </summary>
        private void Opt()
        {
            var M = Matrix<double>.Build; // 生成矩阵
            int FEs = 0; // 当前评价次数
            // 开始寻优
            while (FEs < Max_FEs)
            {
                // Generate and evaluate lambda offspring
                for (int i = 0; i < X_Dim; i++)
                {
                    for (int k = 0; k < lambda; k++)
                    {
                        arz[i, k] = G.rnd_normal(0.0, 1.0, rd); // 正态分布
                    }
                }
                arx = xmeanw * M.Dense(1, lambda, 1) + sigma * (BD * arz);
                // Handle the elements of the variable which violate the boundary
                for (int i = 0; i < X_Dim; i++)
                {
                    for (int k = 0; k < lambda; k++)
                    {
                        arx[i, k] = G.BoundConstraint_CMA_ES(arx[i, k], X_Lb[i], X_Ub[i]);
                    }
                }
                for (int k = 0; k < lambda; k++)
                {
                    for (int i = 0; i < X_Dim; i++)
                    {
                        X[i] = arx[i, k];
                    }
                    Fun[k] = ObjectFunction(X, out Delta[k]);
                    Index[k] = k;
                }

                // Sort by fitness and compute weighted mean in xmeanw
                BubbleSort();
                for (int i = 0; i < X_Dim; i++)
                {
                    xold[i, 0] = xmeanw[i, 0];
                }
                for (int k = 0; k < mu; k++)
                {
                    for (int i = 0; i < X_Dim; i++)
                    {
                        arx_mu[i, k] = arx[i, Index[k]];
                        arz_mu[i, k] = arz[i, Index[k]];
                    }
                }
                xmeanw = arx_mu * arweights / arweights.RowSums().Sum();
                zmeanw = arz_mu * arweights / arweights.RowSums().Sum();

                // Adapt covariance matrix
                pc = (1 - cc) * pc + (Math.Sqrt(cc * (2 - cc)) * cw / sigma) * (xmeanw - xold);
                // Do not adapt in the initial phase
                if (!flginiphase)
                {
                    C = (1 - ccov) * C + ccov * pc * pc.Transpose();
                }
                // Adapt sigma
                ps = (1 - cs) * ps + (Math.Sqrt(cs * (2 - cs)) * cw) * (B * zmeanw);
                sigma = sigma * Math.Exp((ps.L2Norm() - chiN) / chiN / damp);

                //Update B and D from C
                if (((FEs / lambda) % (1 / ccov / X_Dim / 5)) < 1)
                {
                    for (int i = 0; i < X_Dim; i++)
                    {
                        for (int k = 0; k < X_Dim; k++)
                        {
                            if (i == k)
                            {
                                CC[i, k] = 0.0;
                            }
                            else
                            {
                                CC[i, k] = C[i, k];
                            }
                        }
                    }
                    C = C.UpperTriangle() + CC.UpperTriangle().Transpose(); // Enforce symmetry
                    try
                    {
                        D = C.Evd().D;
                        B = C.Evd().EigenVectors;
                    }
                    catch (System.Exception ex)
                    {
                        break;
                    }
                    // Limit condition of C to 1e14 + 1
                    if (D.Diagonal().Max() > 1.0e14 * D.Diagonal().Min())
                    {
                        double tmp1 = D.Diagonal().Max() / 1.0e14 - D.Diagonal().Min();
                        C = C + tmp1 * M.DenseIdentity(X_Dim);
                        D = D + tmp1 * M.DenseIdentity(X_Dim);
                    }
                    double[] tmp2 = new double[X_Dim];
                    for (int i = 0; i < X_Dim; i++)
                    {
                        tmp2[i] = Math.Sqrt(D.Diagonal().ToArray()[i]);
                    }
                    D = M.DenseOfDiagonalArray(tmp2); // D contains standard deviations now
                    BD = B * D; // For speed up only
                }

                // Adjust minimal step size
                if ((sigma * D.Diagonal().Min() < Min_sigma)
                    || (Fun[0] == Fun[Math.Min(mu, lambda - 1)])
                    || (xmeanw == xmeanw + 0.2 * sigma * BD.SubMatrix(0, X_Dim, (int)Math.Floor((double)(FEs / lambda) % X_Dim), 1)))
                {
                    sigma = 1.4 * sigma;
                }
                if (sigma > Max_sigma)
                {
                    sigma = Max_sigma;
                }

                // Test for end of initial phase
                if (flginiphase && ((FEs / lambda) > (2 / cs)))
                {
                    // Step size is not much too small
                    if (((ps.L2Norm() - chiN) / chiN) < 0.05)
                    {
                        flginiphase = false;
                    }
                }

                // 种群最优解
                if (FEs > 0)
                {
                    TolFun = Math.Abs((Last_Best_Fun - Best_Fun) / Last_Best_Fun); // 最优适应度相对变化
                    if ((Best_Delta == Delta[0] && Fun[0] < Best_Fun) || Delta[0] < Best_Delta)
                    {
                        for (int j = 0; j < X_Dim; j++)
                        {
                            X_Opt[j] = arx[j, Index[0]];
                        }
                        Best_Fun = Fun[0];
                        Best_Delta = Delta[0];
                    }
                }
                else
                {
                    for (int j = 0; j < X_Dim; j++)
                    {
                        X_Opt[j] = arx[j, Index[0]];
                    }
                    Best_Fun = Fun[0];
                    Best_Delta = Delta[0];
                }
                Last_Best_Fun = Best_Fun;
                FEs = FEs + lambda;
                //if (TolFun < Min_TolFun)
                //{
                //    break;
                //}
                Console.WriteLine("CMA_ES " + FEs + " " + Best_Fun.ToString("f2") + " " + Best_Delta.ToString("f2"));
            }
        }

        /// <summary>
        /// 冒泡排序操作
        /// </summary>
        private void BubbleSort()
        {
            bool flag;
            double Fun_tmp = 0.0;
            double Delta_tmp = 0.0;
            int Index_tmp = 0;
            for (int i = lambda - 1; i > 0; i--)
            {
                flag = true;
                for (int j = 0; j < i; j++)
                {
                    if ((Delta[j] == Delta[j + 1] && Fun[j] > Fun[j + 1]) || Delta[j] > Delta[j + 1])
                    {
                        Fun_tmp = Fun[j];
                        Fun[j] = Fun[j + 1];
                        Fun[j + 1] = Fun_tmp;
                        Delta_tmp = Delta[j];
                        Delta[j] = Delta[j + 1];
                        Delta[j + 1] = Delta_tmp;
                        Index_tmp = Index[j];
                        Index[j] = Index[j + 1];
                        Index[j + 1] = Index_tmp;
                        flag = false;
                    }
                }
                if (flag) break;
            }
        }

        /// <summary>
        /// 目标函数
        /// </summary>
        /// <param name="X">目标函数输入——优化变量</param>
        /// <param name="model">模型接口</param>
        /// <param name="Delta">违反度，输出参数，未违反约束时为0</param>
        /// <returns>适应度函数值</returns>
        private double ObjectFunction(double[] X, out double Delta)
        {
            return fitness.Fitness(X, model, out Delta);
            //ObjectFunctions O = new ObjectFunctions();
            //return O.ObjectFunctionWithCon(X, out Delta);
        }
    }
}