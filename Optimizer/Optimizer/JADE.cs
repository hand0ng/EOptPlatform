using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Accord.Statistics;
using Accord.Statistics.Distributions.Univariate;
using InterfaceDeclare;

namespace OptimizerSet
{
    /// <summary>
    /// JADE: Adaptive Differential Evolution with Optional External Archive, TEVC, 2009
    /// </summary>
    public class JADE
    {
        #region 私有变量
        private Random rd; // 随机种子
        private int Max_FEs = 30000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private int N_Pop = 100; // 种群规模
        private double c = 0.1;
        private double p = 0.05;
        private double mu_F = 0.5; // 缩放比例因子
        private double mu_CR = 0.5; // 交叉概率因子
        private double F = 0.0; // 缩放比例因子
        private double CR = 0.0; // 交叉概率因子
        private List<double[]> Archive; // 变量存档
        private List<double> S_F; // F存档
        private List<double> S_CR; // CR存档
        private int X_Dim = 0; // 变量维度
        private int P_Pop = 0;
        private double[,] X; // 种群
        private double[,] X_P; // 种群
        private double[] X_Init; // 优化变量初始值
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] Fun; // 各个个体的适应度
        private double[] Delta; // 各个个体的约束违反度
        private double[] XBest; // 用于变异的最优X
        private double[] t; // 目标个体
        private double[] v; // 测试个体
        private double[] X_Opt; // 最优解
        private int Best_Index = 0; // 最优个体编号
        private double Fun_tmp = 0.0; // 当前测试个体的适应度
        private double Delta_tmp = 0.0; // 当前测试个体的约束违反度
        private double Last_Best_Fun = 0.0; // 上次最优个体适应度
        private double Best_Fun = 0.0; // 最优个体适应度
        private double Best_Delta = 0.0; // 最优个体约束违反度
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
                this.c = double.Parse(Config[3][1]);
                this.p = double.Parse(Config[4][1]);
                this.mu_F = double.Parse(Config[5][1]); // 缩放比例因子
                this.mu_CR = double.Parse(Config[6][1]); // 交叉概率因子
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.Archive = new List<double[]>(); // 变量存档
            this.S_F = new List<double>(); // F存档
            this.S_CR = new List<double>(); // CR存档
            this.X_Dim = X_Init.Length;
            this.P_Pop = (int)Math.Ceiling(this.p * this.N_Pop);
            this.X = new double[N_Pop, X_Dim];
            this.X_P = new double[P_Pop, X_Dim];
            this.X_Init = new double[X_Dim];
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.Fun = new double[N_Pop]; // 各个个体的适应度
            this.Delta = new double[N_Pop]; // 各个个体的约束违反度
            this.XBest = new double[X_Dim];
            this.v = new double[X_Dim];
            this.t = new double[X_Dim];
            this.X_Opt = new double[X_Dim];
            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Init[i] = X_Init[i];
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
            }
            this.rd = new Random();
            Opt(); // 优化
            Fun = Best_Fun;
            Delta = Best_Delta;
            return X_Opt;
        }

        /// <summary>
        /// JADE
        /// </summary>
        private void Opt()
        {
            int r1, r2;
            int FEs = 0; // 当前评价次数
            // 上下限内初始化种群
            Initial();
            // 开始寻优
            // 种群停止进化，结束优化流程
            while (FEs < Max_FEs)
            {
                // 选择父代个体
                for (int i = 0; i < N_Pop; i++)
                {
                    // 随机选择另外两个解执行变异操作生成临时解，且i≠r1≠r2，随机个体索引
                    r1 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r2 = G.rndInt_uni(0, N_Pop + Archive.Count() - 1, rd);
                    while (i == r1) r1 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r2 == r1 || r2 == i) r2 = G.rndInt_uni(0, N_Pop + Archive.Count() - 1, rd);

                    // 计算F和CR
                    F = CauchyDistribution.Random(mu_F, 0.1); // 柯西分布
                    F = Math.Min(F, 1.0);
                    //F = 0.5;
                    while (F < 0.0)
                    {
                        F = CauchyDistribution.Random(mu_F, 0.1); // 柯西分布
                        F = Math.Min(F, 1.0);
                    }
                    CR = G.rnd_normal(mu_CR, 0.1, rd); // 正态分布
                    CR = Math.Min(1, Math.Max(0.0, CR));

                    // 从X_P中得到Xbest
                    GetXbest();

                    // 变异操作_目标个体
                    Mutation(i, r1, r2);

                    // 交叉操作_测试个体
                    Cross(i);

                    // 选择操作
                    Select(i);

                    // 更新mu_F,mu_CR,Archive
                    Update();
                }

                // 种群最优解
                GlobalOpt();
                if (FEs > 0)
                {
                    TolFun = Math.Abs((Last_Best_Fun - Best_Fun) / Last_Best_Fun); // 最优适应度相对变化
                }
                Last_Best_Fun = Best_Fun;
                FEs = FEs + N_Pop;
                //if (TolFun < Min_TolFun)
                //{
                //    break;
                //}
                Console.WriteLine("JADE " + FEs + " " + Best_Fun.ToString("f2") + " " + Best_Delta.ToString("f2"));
            }
        }

        /// <summary>
        /// 获取XBest
        /// </summary>
        private void GetXbest()
        {
            bool flag;
            double Fun_tmp = 0.0;
            double Delta_tmp = 0.0;
            double[] Fun_P = new double[N_Pop];
            double[] Delta_P = new double[N_Pop];
            int[] Index = new int[N_Pop];
            int Index_tmp = 0;
            int XBest_Index;
            for (int i = 0; i < N_Pop; i++)
            {
            	Fun_P[i] = Fun[i];
                Delta_P[i] = Delta[i];
                Index[i] = i;
            }
            for (int i = N_Pop - 1; i > 0; i--)
            {
                flag = true;
                for (int j = 0; j < i; j++)
                {
                    if ((Delta_P[j] == Delta_P[j + 1] && Fun_P[j] > Fun_P[j + 1]) || Delta_P[j] > Delta_P[j + 1])
                    {
                        Fun_tmp = Fun_P[j];
                        Fun_P[j] = Fun_P[j + 1];
                        Fun_P[j + 1] = Fun_tmp;
                        Delta_tmp = Delta_P[j];
                        Delta_P[j] = Delta_P[j + 1];
                        Delta_P[j + 1] = Delta_tmp;
                        Index_tmp = Index[j];
                        Index[j] = Index[j + 1];
                        Index[j + 1] = Index_tmp;
                        flag = false;
                    }
                }
                if (flag) break;
            }
            for (int i = 0; i < P_Pop; i++)
            {
                for (int j = 0; j < X_Dim; j++)
                {
                	X_P[i, j] = X[Index[i], j];
                }
            }
            XBest_Index = G.rndInt_uni(0, P_Pop - 1, rd);
            for (int i = 0; i < X_Dim; i++)
            {
                XBest[i] = X_P[XBest_Index, i];
            }
        }

        /// <summary>
        /// 变异操作,DE/current-to-pbest/1
        /// </summary>
        private void Mutation(int Index, int r1, int r2)
        {
            if (r2 > N_Pop - 1)
            {
                for (int i = 0; i < X_Dim; i++)
                {
                    t[i] = X[Index, i] + F * (XBest[i] - X[Index, i]) + F * (X[r1, i] - Archive[r2 - N_Pop][i]);
                    t[i] = G.BoundConstraint_JADE(X[Index, i], t[i], X_Lb[i], X_Ub[i]);
                }
            }
            else
            {
                for (int i = 0; i < X_Dim; i++)
                {
                    t[i] = X[Index, i] + F * (XBest[i] - X[Index, i]) + F * (X[r1, i] - X[r2, i]);
                    t[i] = G.BoundConstraint_JADE(X[Index, i], t[i], X_Lb[i], X_Ub[i]);
                }
            }
        }

        /// <summary>
        /// 二项式交叉操作
        /// </summary>
        private void Cross(int Index)
        {
            int D = G.rndInt_uni(0, X_Dim - 1, rd);
            for (int i = 0; i < X_Dim; i++)
            {
                if (G.rnd_uni(rd) < CR || i == D)
                {
                    v[i] = t[i];
                }
                else
                {
                    v[i] = X[Index, i];
                }
            }
        }

        /// <summary>
        /// 选择操作
        /// </summary>
        private void Select(int Index)
        {
            // 计算测试个体的适应度与约束违反度
            Fun_tmp = ObjectFunction(v, out Delta_tmp);
            // 选择适应度与约束违反度小的个体
            if ((Delta_tmp == Delta[Index] && Fun_tmp < Fun[Index]) || Delta_tmp < Delta[Index])
            {
                double[] X_tmp = new double[X_Dim];
                double F_tmp = F;
                double CR_tmp = CR;
                for (int j = 0; j < X_Dim; j++)
                {
                    X_tmp[j] = X[Index, j];
                    X[Index, j] = v[j];
                }
                Fun[Index] = Fun_tmp;
                Delta[Index] = Delta_tmp;
                S_F.Add(F_tmp);
                S_CR.Add(CR_tmp);
                Archive.Add(X_tmp);
            }
        }

        /// <summary>
        /// 更新mu_F,mu_CR,Archive
        /// </summary>
        private void Update()
        {
            // 更新mu_F,mu_CR
            double F_SqureSum = 0.0;
            double F_Sum = 0.0;
            double F_LehmerMean = 0.0;
            double CR_Mean = 0.0;
            for (int i = 0; i < S_F.Count; i++)
            {
                F_SqureSum = F_SqureSum + Math.Pow(S_F[i], 2.0);
                F_Sum = F_Sum + S_F[i];
            }
            if (F_Sum > 0.0)
            {
                F_LehmerMean = F_SqureSum / F_Sum;
            }
            mu_F = (1 - c) * mu_F + c * F_LehmerMean;
            if (S_CR.Count > 0)
            {
                CR_Mean = S_CR.Average();
            }
            mu_CR = (1 - c) * mu_CR + c * CR_Mean;

            if (Archive.Count == 0)
            {
                return;
            }
            // 删除重复数据
            //List<double[]> Archive_Tmp = new List<double[]>();
            //foreach (double[] item in Archive)
            //{
            //    if (!Archive_Tmp.Contains(item))
            //    {
            //        Archive_Tmp.Add(item);
            //    }
            //}
            //Archive = new List<double[]>(Archive_Tmp);
            //Archive = Archive.Distinct(new ComparerClass()).ToList();
            int[] Index_Repeat = new int[Archive.Count];
            for (int i = 0; i < Archive.Count; i++)
            {
                if (Index_Repeat[i] == 0)
                {
                    for (int j = i + 1; j < Archive.Count; j++)
                    {
                        if (Index_Repeat[j] == 0)
                        {
                            Index_Repeat[j] = 1;
                            for (int k = 0; k < Archive[0].Length; k++)
                            {
                                if (Archive[i][k] != Archive[j][k])
                                {
                                    Index_Repeat[j] = 0;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if (Index_Repeat.Sum()>0)
            {
                for (int i = Archive.Count - 1; i >= 0; i--)
                {
                    if (Index_Repeat[i] == 1)
                    {
                        Archive.RemoveAt(i);
                    }
                }
            }

            // 删除多余数据
            while (Archive.Count > N_Pop)
            {
                int Index = G.rndInt_uni(0, Archive.Count - 1, rd);
                Archive.RemoveAt(Index);
            }
        }

        /// <summary>
        /// 种群最优解
        /// </summary>
        private void GlobalOpt()
        {
            bool update = false;
            // 得到最优适应度 
            for (int i = 0; i < N_Pop; i++)
            {
                if ((Delta[i] == Best_Delta && Fun[i] < Best_Fun) || Delta[i] < Best_Delta)
                {
                    Best_Fun = Fun[i];
                    Best_Delta = Delta[i];
                    Best_Index = i;
                    update = true;
                }
            }
            // 得到最优解
            if (update)
            {
                for (int i = 0; i < X_Dim; i++)
                {
                    X_Opt[i] = X[Best_Index, i];
                }
            }
        }

        /// <summary>
        /// 随机初始化种群
        /// </summary>
        private void Initial()
        {
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
                        X[i, k] = G.rnd_uni(rd) * (this.X_Ub[k] - this.X_Lb[k]) + X_Lb[k];
                    }
                    X_tmp[k] = X[i, k];
                }
                Fun[i] = ObjectFunction(X_tmp, out Delta[i]); // 得到种群适应度和约束违反度
            }
            Best_Fun = Fun[0];
            Best_Delta = Delta[0];
            Best_Index = 0;
            // 得到最优适应度 
            for (int i = 1; i < N_Pop; i++)
            {
                if ((Delta[i] == Best_Delta && Fun[i] < Best_Fun) || Delta[i] < Best_Delta)
                {
                    Best_Fun = Fun[i];
                    Best_Delta = Delta[i];
                    Best_Index = i;
                }
            }
            // 得到最优解
            for (int i = 0; i < X_Dim; i++)
            {
                X_Opt[i] = X[Best_Index, i];
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

    class ComparerClass : System.Collections.Generic.IEqualityComparer<double[]>
    {
        public bool Equals(double[] a, double[] b)
        {
            if (a == null && b == null)
                return true;
            else if (a != null && b != null)
            {
                if (a.Length != b.Length)
                    return false;
                a = a.OrderBy(t => t).ToArray();
                b = b.OrderBy(t => t).ToArray();
                for (int i = 0; i < a.Length; i++)
                    if (a[i] != b[i])
                        return false;
                return true;
            }
            else return false;
        }
        public int GetHashCode(double[] ary)
        {
            return base.GetHashCode();
        }
    }
}