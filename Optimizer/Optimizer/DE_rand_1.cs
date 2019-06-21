using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using InterfaceDeclare;

namespace OptimizerSet
{
    public class DE_rand_1
    {
        #region 私有变量
        private Random rd; // 随机种子
        private int Max_FEs = 30000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private int N_Pop = 100; // 种群规模
        private double F = 0.7; // 缩放比例因子
        private double CR = 0.5; // 交叉概率因子
        private int X_Dim = 0; // 变量维度
        private double[,] X; // 种群
        private double[] X_Init; // 优化变量初始值
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] Fun; // 各个个体的适应度
        private double[] Delta; // 各个个体的约束违反度
        private double[] t; // 目标个体
        private double[] v; // 测试个体
        private double[] X_Opt; // 最优解
        private int Mutation_Sel = 1; // DE变异算子选择，DE/rand/1，DE/best/1，DE/rand/2，DE/best/2
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
                this.F = double.Parse(Config[3][1]); // 缩放比例因子
                this.CR = double.Parse(Config[4][1]); // 交叉概率因子
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.X_Dim = X_Init.Length;
            this.X = new double[N_Pop, X_Dim];
            this.X_Init = new double[X_Dim];
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.Fun = new double[N_Pop]; // 各个个体的适应度
            this.Delta = new double[N_Pop]; // 各个个体的约束违反度
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
        /// DE
        /// </summary>
        private void Opt()
        {
            int r1, r2, r3, r4, r5;
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
                    // 随机选择另外三个解执行变异操作生成临时解，且i≠r1≠r2≠r3，随机个体索引
                    r1 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r2 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r3 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (i == r1) r1 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r2 == r1 || r2 == i) r2 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r3 == r1 || r3 == r2 || r3 == i) r3 = G.rndInt_uni(0, N_Pop - 1, rd);

                    r4 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r5 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r4 == r1 || r4 == r2 || r4 == r3 || r4 == i) r4 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4 || r5 == i) r5 = G.rndInt_uni(0, N_Pop - 1, rd);

                    // 变异操作_目标个体
                    Mutation(r1, r2, r3, r4, r5, Mutation_Sel);

                    // 交叉操作_测试个体
                    Cross(i);

                    // 选择操作
                    Select(i);
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
                Console.WriteLine("DE_rand_1 " + FEs + " " + Best_Fun.ToString("f2") + " " + Best_Delta.ToString("f2"));
            }
        }

        /// <summary>
        /// 变异操作,rand开采能力更强，best局部寻优能力更强,实现DE/rand/1，DE/best/1，DE/rand/2，DE/best/2
        /// </summary>
        private void Mutation(int r1, int r2, int r3, int r4, int r5, int Mutation_Sel)
        {
            switch (Mutation_Sel)
            {
                case 1: // DE/rand/1
                    {
                        // 生成变异个体
                        for (int i = 0; i < X_Dim; i++)
                        {
                            t[i] = X[r1, i] + F * (X[r2, i] - X[r3, i]);
                            t[i] = G.BoundConstraint_DE(X[r1, i], t[i], X_Lb[i], X_Ub[i]);
                        }
                        break;
                    }
                case 2: // DE/best/1
                    {
                        // 从F中选出best个体_适应度最小
                        double value = Fun[0];
                        int b = 0;
                        for (int i = 1; i < N_Pop; i++)
                        {
                            if (value > Fun[i])
                            {
                                value = Fun[i];
                                b = i;
                            }
                        }
                        // 生成变异个体
                        for (int i = 0; i < X_Dim; i++)
                        {
                            t[i] = X[b, i] + F * (X[r2, i] - X[r3, i]);
                            t[i] = G.BoundConstraint_DE(X[b, i], t[i], X_Lb[i], X_Ub[i]);
                        }
                        break;
                    }
                case 3: // DE/rand/2
                    {
                        // 生成变异个体
                        for (int i = 0; i < X_Dim; i++)
                        {
                            t[i] = X[r1, i] + F * (X[r2, i] - X[r3, i]) + F * (X[r4, i] - X[r5, i]);
                            t[i] = G.BoundConstraint_DE(X[r1, i], t[i], X_Lb[i], X_Ub[i]);
                        }
                        break;
                    }
                case 4: // DE/best/2
                    {
                        // 从F中选出best个体_适应度最小
                        double value = Fun[0];
                        int b = 0;
                        for (int i = 1; i < N_Pop; i++)
                        {
                            if (value > Fun[i])
                            {
                                value = Fun[i];
                                b = i;
                            }
                        }
                        // 生成变异个体
                        for (int i = 0; i < X_Dim; i++)
                        {
                            t[i] = X[b, i] + F * (X[r1, i] - X[r2, i]) + F * (X[r3, i] - X[r4, i]);
                            t[i] = G.BoundConstraint_DE(X[b, i], t[i], X_Lb[i], X_Ub[i]);
                        }
                        break;
                    }
            }
        }

        /// <summary>
        /// 二项式交叉操作
        /// </summary>
        private void Cross(int Index)
        {
            for (int i = 0; i < X_Dim; i++)
            {
                if (G.rnd_uni(rd) < CR || i == X_Dim)
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
                for (int j = 0; j < X_Dim; j++)
                {
                    X[Index, j] = v[j];
                }
                Fun[Index] = Fun_tmp;
                Delta[Index] = Delta_tmp;
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
}