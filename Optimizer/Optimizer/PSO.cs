using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using InterfaceDeclare;

namespace OptimizerSet
{
    public class PSO
    {
        #region 私有变量
        private Random rd; // 随机种子
        private int Max_FEs = 30000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private int N_Pop = 100; // 种群规模
        private double W = 0.0; // 惯性权重
        private double Wmax = 0.9; // 惯性权重Max
        private double Wmin = 0.4; // 惯性权重Min
        private double C1 = 1.4; // 加速因子1
        private double C2 = 1.4; // 加速因子2
        private int X_Dim = 0; // 变量维度
        private double[,] X; // 种群
        private double[,] V; // 个体速度
        private double[] X_Init; // 优化变量初始值
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] Vmax; // 最大速度
        private double[] Vmin; // 最小速度
        private double[] Fun; // 各个个体的适应度
        private double[] Delta; // 各个个体的约束违反度
        private double[] gBest; // 全局最优个体
        private double gBest_Fun = 0.0; // 全局最优个体目标函数值
        private double gBest_Delta = 0.0; // 全局最优个体约束违反度
        private double[,] pBest; // 个体历史最优个体
        private double[] pBest_Fun; // 个体历史最优个体目标函数值
        private double[] pBest_Delta; // 个体历史最优个体约束违反度
        private double[] X_Opt; // 最优解
        private double[] res;  // 返回值
        private double Last_gBest_Fun = 0.0; // 上次最优个体适应度
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
                this.N_Pop = int.Parse(Config[2][1]); // 种群规模
                this.Wmax = double.Parse(Config[3][1]); // 惯性权重Max
                this.Wmin = double.Parse(Config[4][1]); // 惯性权重Min
                this.C1 = double.Parse(Config[5][1]); // 加速因子1
                this.C2 = double.Parse(Config[6][1]); // 加速因子2
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.X_Dim = X_Init.Length;
            this.X = new double[N_Pop, X_Dim];
            this.V = new double[N_Pop, X_Dim];
            this.X_Init = new double[X_Dim];
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.Vmax = new double[X_Dim];
            this.Vmin = new double[X_Dim];
            this.Fun = new double[N_Pop]; // 各个个体的适应度
            this.Delta = new double[N_Pop]; // 各个个体的约束违反度
            this.gBest = new double[X_Dim];
            this.pBest = new double[N_Pop, X_Dim];
            this.pBest_Fun = new double[N_Pop];
            this.pBest_Delta = new double[N_Pop];
            this.X_Opt = new double[X_Dim];
            this.res = new double[X_Dim + 1];
            this.rd = new Random();
            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Init[i] = X_Init[i];
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
                this.Vmax[i] = (this.X_Ub[i] - this.X_Lb[i]) * 0.5;
                this.Vmin[i] = this.Vmax[i] * -1.0;
                for (int k = 0; k < N_Pop; k++)
                {
                    this.V[k, i] = this.Vmin[i] + G.rnd_uni(rd) * (this.Vmax[i] - this.Vmin[i]);
                }
            }
            Opt(); // 优化
            Fun = gBest_Fun;
            Delta = gBest_Delta;
            return X_Opt;
        }

        /// <summary>this</summary>
        private void Opt()
        {
            double[] X_tmp = new double[X_Dim];
            int FEs = 0; // 当前评价次数
            // 上下限内初始化种群
            Initial();
            // 开始寻优
            // 种群停止进化，结束优化流程
            while (FEs < Max_FEs)
            {
                W = Wmax - FEs / Max_FEs * (Wmax - Wmin);
                // 选择个体
                for (int i = 0; i < N_Pop; i++)
                {
                    // 粒子i速度与位置更新
                    for (int k = 0; k < X_Dim; k++)
                    {
                        V[i, k] = W * V[i, k] + C1 * (pBest[i, k] - X[i, k]) + C2 * (gBest[k] - X[i, k]);
                        V[i, k] = G.BoundConstraint(V[i, k], Vmin[k], Vmax[k]);
                        X[i, k] = X[i, k] + V[i, k];
                        X[i, k] = G.BoundConstraint(X[i, k], X_Lb[k], X_Ub[k]);
                        X_tmp[k] = X[i, k];
                    }
                    Fun[i] = ObjectFunction(X_tmp, out Delta[i]);
                    // 更新pBest
                    if ((Delta[i] == gBest_Delta && Fun[i] < pBest_Fun[i]) || Delta[i] < pBest_Delta[i])
                    {
                        for (int k = 0; k < X_Dim; k++)
                        {
                            pBest[i, k] = X[i, k];
                        }
                        pBest_Fun[i] = Fun[i];
                        pBest_Delta[i] = Delta[i];
                    }
                    // 更新gBest
                    if ((Delta[i] == gBest_Delta && Fun[i] < gBest_Fun) || Delta[i] < gBest_Delta)
                    {
                        for (int k = 0; k < X_Dim; k++)
                        {
                            gBest[k] = X[i, k];
                        }
                        gBest_Fun = Fun[i];
                        gBest_Delta = Delta[i];
                    }
                }
                if (FEs > 0)
                {
                    TolFun = Math.Abs((Last_gBest_Fun - gBest_Fun) / Last_gBest_Fun); // 最优适应度变化
                }
                Last_gBest_Fun = gBest_Fun;
                FEs = FEs + N_Pop;
                //if (TolFun < Min_TolFun)
                //{
                //    break;
                //}
                Console.WriteLine("PSO " + FEs + " " + gBest_Fun.ToString("f2") + " " + gBest_Delta.ToString("f2"));
            }
            // 得到最优解
            for (int i = 0; i < X_Dim; i++)
            {
                X_Opt[i] = gBest[i];
            }
        }

        /// <summary>随机初始化种群</summary>
        private void Initial()
        {
            int gBentIndex = 0;
            double[] X_tmp = new double[X_Dim];
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
                        X[i, k] = G.rnd_uni(rd) * (X_Ub[k] - X_Lb[k]) + X_Lb[k];
                    }
                    pBest[i, k] = X[i, k]; // 个体历史最优个体
                    X_tmp[k] = X[i, k];
                }
                Fun[i] = ObjectFunction(X_tmp, out Delta[i]);
                pBest_Fun[i] = Fun[i]; // 个体历史最优个体目标函数值
                pBest_Delta[i] = Delta[i]; // 个体历史最优个体约束违反度
            }
            gBest_Fun = Fun[0];
            gBest_Delta = Delta[0];
            for (int i = 1; i < N_Pop; i++)
            {
                if (Fun[i] < gBest_Fun && Delta[i] <= gBest_Delta)
                {
                    gBest_Fun = Fun[i];
                    gBest_Delta = Delta[i];
                    gBentIndex = i;
                }
            }
            for (int i = 0; i < X_Dim; i++)
            {
                gBest[i] = X[gBentIndex, i];
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