using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using InterfaceDeclare;

namespace OptimizerSet
{
    /// <summary>
    /// Differential Evolution with Composite Trial Vector Generation Strategies and Control Parameters, TEVC, 2011
    /// </summary>
    public class CoDE
    {
        #region 私有变量
        private Random rd; // 随机种子
        private int Max_FEs = 30000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private int N_Pop = 30; // 种群规模
        private double[] F = new double[3] { 1.0, 1.0, 0.8 }; // 缩放比例因子
        private double[] CR = new double[3] { 0.1, 0.9, 0.2 }; // 交叉概率因子
        private int X_Dim = 0; // 变量维度
        private double[,] X; // 种群
        private double[] X_Init; // 优化变量初始值
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] Fun; // 各个个体的适应度
        private double[] Delta; // 各个个体的约束违反度
        private double[] v; // 测试个体
        private double[] t1; // 目标个体1
        private double[] v1; // 测试个体1
        private double[] t2; // 目标个体2
        private double[] v2; // 测试个体2
        private double[] t3; // 目标个体3
        private double[] v3; // 测试个体3
        private double[] X_Opt; // 最优解
        private int Best_Index = 0; // 最优个体编号
        private double Fun_tmp = 0.0; // 测试个体的适应度
        private double Delta_tmp = 0.0; // 测试个体的约束违反度
        private double Fun_tmp1 = 0.0; // 测试个体1的适应度
        private double Delta_tmp1 = 0.0; // 测试个体1的约束违反度
        private double Fun_tmp2 = 0.0; // 测试个体2的适应度
        private double Delta_tmp2 = 0.0; // 测试个体2的约束违反度
        private double Fun_tmp3 = 0.0; // 测试个体3的适应度
        private double Delta_tmp3 = 0.0; // 测试个体3的约束违反度
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
            this.v1 = new double[X_Dim];
            this.t1 = new double[X_Dim];
            this.v2 = new double[X_Dim];
            this.t2 = new double[X_Dim];
            this.v3 = new double[X_Dim];
            this.t3 = new double[X_Dim];
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
        /// CoDE
        /// </summary>
        private void Opt()
        {
            int rndNum;
            int r11, r12, r13, r21, r22, r23, r24, r25, r31, r32, r33;
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
                    r11 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r12 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r13 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (i == r11) r11 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r12 == r11 || r12 == i) r12 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r13 == r11 || r13 == r12 || r13 == i) r13 = G.rndInt_uni(0, N_Pop - 1, rd);

                    r21 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r22 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r23 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r24 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r25 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (i == r21) r21 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r22 == r21 || r22 == i) r22 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r23 == r21 || r23 == r22 || r23 == i) r23 = G.rndInt_uni(0, N_Pop - 1, rd);
                    while (r24 == r21 || r24 == r22 || r24 == r23 || r24 == i) r24 = G.rndInt_uni(0, N_Pop - 1, rd);
                    if (N_Pop > 5)
                    {
                        while (r25 == r21 || r25 == r22 || r25 == r23 || r25 == r24 || r25 == i) r25 = G.rndInt_uni(0, N_Pop - 1, rd);
                    }

                    r31 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r32 = G.rndInt_uni(0, N_Pop - 1, rd);
                    r33 = G.rndInt_uni(0, N_Pop - 1, rd);

                    // rand/1/bin生成t1，v1
                    rndNum = G.rndInt_uni(0, 2, rd);
                    // 变异操作_目标个体
                    Mutation1(r11, r12, r13, F[rndNum]);
                    Cross(i, CR[rndNum], t1, v1);

                    // rand/2/bin生成t2，v2
                    rndNum = G.rndInt_uni(0, 2, rd);
                    // 变异操作_目标个体
                    Mutation2(r21, r22, r23, r24, r25, F[rndNum]);
                    Cross(i, CR[rndNum], t2, v2);

                    // current-to-rand/1生成t3，v3
                    rndNum = G.rndInt_uni(0, 2, rd);
                    // 变异操作_目标个体
                    Mutation3(i, r31, r32, r33, F[rndNum]);
                    for (int k = 0; k < t3.Length; k++)
                    {
                        v3[k] = t3[k];
                    }

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
                FEs = FEs + N_Pop * 3;
                //if (TolFun < Min_TolFun)
                //{
                //    break;
                //}
                Console.WriteLine("CoDE " + FEs + " " + Best_Fun.ToString("f2") + " " + Best_Delta.ToString("f2"));
            }
        }

        /// <summary>
        /// 变异操作rand/1
        /// </summary>
        private void Mutation1(int r1, int r2, int r3, double f)
        {
            // 生成变异个体
            for (int i = 0; i < X_Dim; i++)
            {
                t1[i] = X[r1, i] + f * (X[r2, i] - X[r3, i]);
                t1[i] = G.BoundConstraint_CoDE(t1[i], X_Lb[i], X_Ub[i]);
            }
         }

        /// <summary>
        /// 变异操作rand/2
        /// </summary>
        private void Mutation2(int r1, int r2, int r3, int r4, int r5, double f)
        {
            // 生成变异个体
            for (int i = 0; i < X_Dim; i++)
            {
                t2[i] = X[r1, i] + f * (X[r2, i] - X[r3, i]) + f * (X[r4, i] - X[r5, i]);
                t2[i] = G.BoundConstraint_CoDE(t2[i], X_Lb[i], X_Ub[i]);
            }
        }

        /// <summary>
        /// 变异操作current-to-rand/1
        /// </summary>
        private void Mutation3(int r, int r1, int r2, int r3, double f)
        {
            // 生成变异个体
            for (int i = 0; i < X_Dim; i++)
            {
                t3[i] = X[r, i] + G.rnd_uni(rd) * (X[r1, i] - X[r, i]) + f * (X[r2, i] - X[r3, i]);
                t3[i] = G.BoundConstraint_CoDE(t3[i], X_Lb[i], X_Ub[i]);
            }
        }

        /// <summary>
        /// 二项式交叉操作
        /// </summary>
        private void Cross(int Index, double cr, double[] t, double[] v)
        {
            for (int i = 0; i < X_Dim; i++)
            {
                if (G.rnd_uni(rd) < cr || i == X_Dim)
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
        private void Select(int index)
        {
            // 计算测试个体1的适应度与约束违反度
            Fun_tmp1 = ObjectFunction(v1, out Delta_tmp1);

            // 计算测试个体2的适应度与约束违反度
            Fun_tmp2 = ObjectFunction(v2, out Delta_tmp2);

            // 计算测试个体3的适应度与约束违反度
            Fun_tmp3 = ObjectFunction(v3, out Delta_tmp3);

            // 选择适应度与约束违反度小的个体
            Fun_tmp = Fun_tmp1;
            Delta_tmp = Delta_tmp1;
            for (int j = 0; j < X_Dim; j++)
            {
                v[j] = v1[j];
            }

            if ((Delta_tmp2 == Delta_tmp && Fun_tmp2 < Fun_tmp) || Delta_tmp2 < Delta_tmp)
            {
                for (int j = 0; j < X_Dim; j++)
                {
                    v[j] = v2[j];
                }
                Fun_tmp = Fun_tmp2;
                Delta_tmp = Delta_tmp2;
            }

            if ((Delta_tmp3 == Delta_tmp && Fun_tmp3 < Fun_tmp) || Delta_tmp3 < Delta_tmp)
            {
                for (int j = 0; j < X_Dim; j++)
                {
                    v[j] = v3[j];
                }
                Fun_tmp = Fun_tmp3;
                Delta_tmp = Delta_tmp3;
            }

            if ((Delta_tmp == Delta[index] && Fun_tmp < Fun[index]) || Delta_tmp < Delta[index])
            {
                for (int j = 0; j < X_Dim; j++)
                {
                    X[index, j] = v[j];
                }
                Fun[index] = Fun_tmp;
                Delta[index] = Delta_tmp;
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
                        X[i, k] = G.rnd_uni(rd) * (X_Ub[k] - X_Lb[k]) + X_Lb[k];
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