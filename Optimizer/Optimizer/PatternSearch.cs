using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using InterfaceDeclare;

namespace OptimizerSet
{
    public class PatternSearch
    {
        #region 私有变量
        private int Max_FEs = 3000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private double Init_Step_Ratio = 2; // 初始步长比
        private double Min_Step_Ratio = 1.0e4; // 最小步长比
        private double Alpha = 1; // 加速因子
        private double Beta = 0.5; // 收缩因子
        private int X_Dim = 0; // 变量个数
        private bool X_Comput_Flag = true; // X计算标志
        private bool Y_Comput_Flag = true; // Y计算标志
        private bool Step_Flag = true; // 步长更新标志
        private bool Out_Flag = false; // 跳出标志
        private double[] X; // 变量
        private double[] Y; // 变量
        private double[] Y_tmp; // 探测变量
        private double[] Last_X; // 上步变量
        private double[] Step; // 步长
        private double[] Min_Step; // 最小步长
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] X_Opt; // 最优解
        private double Fun_X = 0.0; // 目标值
        private double Last_Fun_X = 0.0; // 目标值
        private double Fun_Y = 0.0; // 目标值
        private double Fun_tmp = 0.0; // 探测目标值
        private double Fun_tmp_Last = 0.0; // 上步探测值
        private double Best_Fun = 0.0; // 最优目标值
        private double Last_Best_Fun = 0.0; // 上一步最优目标值
        private double Delta_X = 0.0; // 约束违反度
        private double Last_Delta_X = 0.0; // 约束违反度
        private double Delta_Y = 0.0; // 约束违反度
        private double Delta_tmp = 0.0; // 探测约束违反度
        private double Delta_tmp_Last = 0.0; // 上步探测约束违反度
        private double Best_Delta = 0.0; // 约束违反度
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
                this.Init_Step_Ratio = double.Parse(Config[2][1]); // 初始步长比
                this.Min_Step_Ratio = double.Parse(Config[3][1]); // 最小步长比
                this.Alpha = double.Parse(Config[4][1]); // 加速因子
                this.Beta = double.Parse(Config[5][1]); // 收缩因子
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.X_Dim = X_Init.Length;
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.X = new double[X_Dim]; // 变量
            this.Y = new double[X_Dim]; // 探测变量
            this.Y_tmp = new double[X_Dim]; // 探测变量
            this.Last_X = new double[X_Dim]; // 模式变量
            this.X_Opt = new double[X_Dim];
            this.Step = new double[X_Dim]; // 步长
            this.Min_Step = new double[X_Dim]; // 最小步长
            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
                this.X[i] = X_Init[i];
                this.Step[i] = (this.X_Ub[i] - this.X_Lb[i]) / this.Init_Step_Ratio;
                this.Min_Step[i] = (this.X_Ub[i] - this.X_Lb[i]) / this.Min_Step_Ratio;
            }
            Opt(); // 优化
            Fun = Best_Fun;
            Delta = Best_Delta;
            return X_Opt;
        }

        /// <summary>
        /// 模式搜索
        /// </summary>
        private void Opt()
        {
            int FEs = 0;
            int iter = 0;
            // Pattern Search
            for (int i = 0; i < X_Dim; i++)
            {
                Y[i] = X[i];
                Last_X[i] = X[i];
                X_Opt[i] = Last_X[i];
            }
            Last_Fun_X = ObjectFunction(X, out Last_Delta_X);
            Best_Fun = Last_Fun_X;
            Best_Delta = Last_Delta_X;
            Out_Flag = false;
            while (FEs < Max_FEs && !Out_Flag)
            {
                // 轴向探测
                if (X_Comput_Flag)
                {
                    //---------------------(F, Delta) = Function(X, obj_index)
                    Fun_X = ObjectFunction(X, out Delta_X);
                    FEs = FEs + 1;
                }
                if (Y_Comput_Flag)
                {
                    //---------------------(F, Delta) = Function(Y, obj_index)
                    Fun_Y = ObjectFunction(Y, out Delta_Y);
                    Fun_tmp_Last = Fun_Y;
                    Delta_tmp_Last = Delta_Y;
                    FEs = FEs + 1;
                }
                for (int i = 0; i < X_Dim; i++)
                {
                    Y_tmp[i] = Y[i];
                }
                for (int i = 0; i < X_Dim; i++)
                {
                    Y_tmp[i] = Y[i] + Step[i]; // 正向探测
                    Y_tmp[i] = G.BoundConstraint_PatternSearch(Y[i], Y_tmp[i], X_Lb[i], X_Ub[i]);
                    //---------------------(F_tmp, Delta_tmp) = Function(Y_tmp, obj_index)
                    Fun_tmp = ObjectFunction(Y_tmp, out Delta_tmp);
                    FEs = FEs + 1;
                    if ((Delta_tmp_Last == 0.0 && Delta_tmp == 0.0 && Fun_tmp_Last <= Fun_tmp) || ((Delta_tmp_Last <= Delta_tmp) && (Delta_tmp_Last + Delta_tmp) != 0.0))
                    {
                        Y_tmp[i] = Y[i] - Step[i]; // 反向探测
                        Y_tmp[i] = G.BoundConstraint_PatternSearch(Y[i], Y_tmp[i], X_Lb[i], X_Ub[i]);
                        //---------------------(F_tmp, Delta_tmp) = Function(Y_tmp, obj_index)
                        Fun_tmp = ObjectFunction(Y_tmp, out Delta_tmp);
                        FEs = FEs + 1;
                        if ((Delta_tmp_Last == 0.0 && Delta_tmp == 0.0 && Fun_tmp_Last <= Fun_tmp) || ((Delta_tmp_Last <= Delta_tmp) && (Delta_tmp_Last + Delta_tmp) != 0.0))
                        {
                            Fun_tmp = Fun_tmp_Last;
                            Delta_tmp = Delta_tmp_Last;
                            Y_tmp[i] = Y[i];
                            X_Comput_Flag = false;
                        }
                    }
                    Fun_tmp_Last = Fun_tmp;
                    Delta_tmp_Last = Delta_tmp;
                }
                // 模式移动
                if ((Delta_tmp_Last == Delta_X && Fun_tmp_Last < Fun_X) || Delta_tmp_Last < Delta_X)
                {
                    X_Comput_Flag = false;
                    Y_Comput_Flag = true;
                    for (int i = 0; i < X_Dim; i++)
                    {
                        Last_X[i] = X[i];
                        X[i] = Y_tmp[i];
                        Y[i] = X[i] + Alpha * (X[i] - Last_X[i]);
                        Y[i] = G.BoundConstraint_PatternSearch(X[i], Y[i], X_Lb[i], X_Ub[i]);
                    }
                    Last_Fun_X = Fun_X;
                    Last_Delta_X = Delta_X;
                    Fun_X = Fun_tmp_Last;
                    Delta_X = Delta_tmp_Last;
                }
                else if (Delta_tmp_Last == Delta_X && Fun_tmp_Last == Fun_X) // 更新步长
                {
                    Step_Flag = true;
                    Out_Flag = false;
                    for (int i = 0; i < X_Dim; i++)
                    {
                        if (X[i] != Y_tmp[i])
                        {
                            Step_Flag = false;
                        }
                    }
                    if (Step_Flag)
                    {
                        X_Comput_Flag = false;
                        Y_Comput_Flag = false;
                        for (int i = 0; i < X_Dim; i++)
                        {
                            Step[i] = Step[i] * Beta;
                            X[i] = Last_X[i];
                            Y[i] = Last_X[i];
                            if (Step[i] < Min_Step[i])
                            {
                                Out_Flag = true;
                            }
                        }
                        Fun_X = Last_Fun_X;
                        Delta_X = Last_Delta_X;
                        Fun_Y = Last_Fun_X;
                        Delta_Y = Last_Delta_X;
                        Fun_tmp_Last = Last_Fun_X;
                        Delta_tmp_Last = Last_Delta_X;
                    }
                    else
                    {
                        X_Comput_Flag = false;
                        Y_Comput_Flag = false;
                        for (int i = 0; i < X_Dim; i++)
                        {
                            X[i] = Y_tmp[i];
                            Y[i] = Y_tmp[i];
                            Last_X[i] = X[i];
                            if (Step[i] < Min_Step[i])
                            {
                                Out_Flag = true;
                            }
                        }
                        Last_Fun_X = Fun_tmp_Last;
                        Last_Delta_X = Delta_tmp_Last;
                        Fun_X = Fun_tmp_Last;
                        Delta_X = Delta_tmp_Last;
                        Fun_Y = Fun_tmp_Last;
                        Delta_Y = Delta_tmp_Last;
                    }
                }
                else
                {
                    Out_Flag = false;
                    X_Comput_Flag = false;
                    Y_Comput_Flag = false;
                    for (int i = 0; i < X_Dim; i++)
                    {
                        Step[i] = Step[i] * Beta;
                        X[i] = Last_X[i];
                        Y[i] = Last_X[i];
                        if (Step[i] < Min_Step[i])
                        {
                            Out_Flag = true;
                        }
                    }
                    Fun_X = Last_Fun_X;
                    Delta_X = Last_Delta_X;
                    Fun_Y = Last_Fun_X;
                    Delta_Y = Last_Delta_X;
                    Fun_tmp_Last = Last_Fun_X;
                    Delta_tmp_Last = Last_Delta_X;
                }
                if ((Delta_X == Best_Delta && Fun_X < Best_Fun) || Delta_X < Best_Delta)
                {
                    for (int i = 0; i < X_Dim; i++)
                    {
                        X_Opt[i] = X[i];
                    }
                    Best_Fun = Fun_X;
                    Best_Delta = Delta_X;
                }
                if (iter > 0)
                {
                    TolFun = Math.Abs((Last_Best_Fun - Best_Fun) / Last_Best_Fun); // 最优适应度相对变化
                }
                iter = iter + 1;
                Last_Best_Fun = Last_Fun_X; // 上一步目标值
                //if (TolFun < Min_TolFun)
                //{
                //    Out_Flag = true;
                //}
                Console.WriteLine("PatternSearch " + FEs + " " + Best_Fun.ToString("f2") + " " + Best_Delta.ToString("f2"));
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