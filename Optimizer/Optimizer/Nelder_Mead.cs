using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using InterfaceDeclare;

namespace OptimizerSet
{
    public class Nelder_Mead
    {
        #region 私有变量
        private int Max_FEs = 3000; // 最大评价次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private double Step_Ratio = 10; // 步长比
        private double Epsilon = 1.0e-4; // 结束条件
        private double Alpha = 1; // 反射因子
        private double Beta = 0.5; // 收缩因子
        private double Gamma = 2; // 扩大因子
        private int X_Dim = 0; // 变量个数
        private bool Out_Flag = false; // 跳出标志
        private double[,] X; // 变量
        private double[] Delta; // 约束违反度
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] X_Center; // 重心点
        private double[] X_R; // 反射点
        private double[] X_E; // 扩大点
        private double[] X_C; // 收缩点
        private double[] X_Opt; // 最优解
        private double[] Step; // 步长
        private double[] Fun; // 目标值
        private double Fun_Min = 0.0; // 最小点目标值
        private double Fun_Max = 0.0; // 最大点目标值
        private double Fun_SecMax = 0.0; // 次最大点目标值
        private double Fun_Center = 0.0; // 重心点目标值
        private double Fun_R = 0.0; // 反射点目标值
        private double Fun_E = 0.0; // 扩大点目标值
        private double Fun_C = 0.0; // 收缩点目标值
        private double Delta_Min = 0.0; // 最小点约束违反度
        private double Delta_Max = 0.0; // 最大点约束违反度
        private double Delta_SecMax = 0.0; // 次最大点约束违反度
        private double Delta_Center = 0.0; // 重心点约束违反度
        private double Delta_R = 0.0; // 反射点约束违反度
        private double Delta_E = 0.0; // 扩大点约束违反度
        private double Delta_C = 0.0; // 收缩点约束违反度
        private int Index_Min = 0;
        private int Index_Max = 0;
        private int Index_SecMax = 0;
        private double Best_Fun = 0.0; // 最优目标值
        private double Last_Best_Fun = 0.0; // 上一步最优目标值
        private double Best_Delta = 0.0; // 最优约束违反度
        private int Best_Index = 0;
        private double TolFun = 1; // 两代最优解之差
        private General G;
        private ObjectFunctions O = new ObjectFunctions();
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
                this.Step_Ratio = double.Parse(Config[2][1]); // 初始步长比
                this.Epsilon = double.Parse(Config[3][1]); // 结束条件
                this.Alpha = double.Parse(Config[4][1]); // 反射因子
                this.Beta = double.Parse(Config[5][1]); // 收缩因子
                this.Gamma = double.Parse(Config[6][1]); // 扩大因子
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.X_Dim = X_Init.Length;
            this.X = new double[X_Dim + 1, X_Dim]; // 变量
            this.Delta = new double[X_Dim + 1]; // 约束违反度
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.X_Center = new double[X_Dim]; // 重心点
            this.X_R = new double[X_Dim]; // 反射点
            this.X_E = new double[X_Dim]; // 扩大点
            this.X_C = new double[X_Dim]; // 收缩点
            this.X_Opt = new double[X_Dim]; // 最优解
            this.Step = new double[X_Dim]; // 步长
            this.Fun = new double[X_Dim + 1]; // 目标值
            
            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
                this.Step[i] = (this.X_Ub[i] - this.X_Lb[i]) / this.Step_Ratio;
                this.X[0, i] = X_Init[i];
                for (int j = 0; j < X_Dim; j++)
                {
                	if (i == j)
                	{
                        this.X[j + 1, i] = X_Init[i] + Step[i];
                        this.X[j + 1, i] = G.BoundConstraint(this.X[j + 1, i], this.X_Lb[i], this.X_Ub[i]);
                	}
                    else
                    {
                        this.X[j + 1, i] = X_Init[i];
                    }
                }  
            }
            Opt(); // 优化
            Fun = Best_Fun;
            Delta = Best_Delta;
            return X_Opt;
        }

        /// <summary>
        /// Nelder_Mead单纯形法
        /// </summary>
        private void Opt()
        {
            var M = Matrix<double>.Build;
            Matrix<double> X_Matrix; // 解矩阵
            int X_Rank = X_Dim; // 解矩阵的秩
            double[] X_tmp = new double[X_Dim];
            double SumEpsilon = 0.0;
            int FEs = 0;
            int iter = 0;
            Initial();
            Out_Flag = false;
            X_Matrix = M.DenseOfArray(X).Transpose();  // 解矩阵
            X_Rank = X_Matrix.Rank(); // 解矩阵的秩
            while (FEs < Max_FEs && !Out_Flag)
            {
                // 计算最小点
                Fun_Min = Fun[0];
                Delta_Min = Delta[0];
                Index_Min = 0;
                for (int i = 0; i < X_Dim + 1; i++)
                {
                    if ((Delta[i] == Delta_Min && Fun[i] < Fun_Min) || Delta[i] < Delta_Min)
                    {
                        Fun_Min = Fun[i];
                        Delta_Min = Delta[i];
                        Index_Min = i;
                    }
                }

                // 计算最大点
                Fun_Max = Fun[0];
                Delta_Max = Delta[0];
                Index_Max = 0;
                for (int i = 0; i < X_Dim + 1; i++)
                {
                    if ((Delta[i] == Delta_Max && Fun[i] > Fun_Max) || Delta[i] < Delta_Max)
                    {
                        Fun_Max = Fun[i];
                        Delta_Max = Delta[i];
                        Index_Max = i;
                    }
                }

                // 计算次最大点
                for (int i = 0; i < X_Dim + 1; i++)
                {
                    if (i != Index_Max)
                    {
                        Index_SecMax = i;
                        break;
                    }
                }
                Fun_SecMax = Fun[Index_SecMax];
                Delta_SecMax = Delta[Index_SecMax];
                for (int i = 0; i < X_Dim + 1; i++)
                {
                    if (i != Index_Max)
                    {
                        if ((Delta[i] == Delta_SecMax && Fun[i] > Fun_SecMax) || Delta[i] < Delta_SecMax)
                        {
                            Fun_SecMax = Fun[i];
                            Delta_SecMax = Delta[i];
                            Index_SecMax = i;
                        }
                    }
                }

                // 计算重心点
                for (int i = 0; i < X_Dim; i++)
                {
                    X_Center[i] = 0.0;
                    for (int j = 0; j < X_Dim + 1; j++)
                    {
                        if (j != Index_Max)
                        {
                            X_Center[i] = X_Center[i] + X[j, i];
                        }
                    }
                    X_Center[i] = X_Center[i] / X_Dim;
                    X_Center[i] = G.BoundConstraint(X_Center[i], this.X_Lb[i], this.X_Ub[i]);
                }
                Fun_Center = ObjectFunction(X_Center, out Delta_Center);
                FEs = FEs + 1;

                // 判断跳出
                SumEpsilon = 0.0;
                for (int i = 0; i < X_Dim + 1; i++)
                {
                    SumEpsilon = SumEpsilon + Math.Pow(Math.Pow(Fun[i] - Fun_Center, 2.0) / (X_Dim + 1), 0.5);
                }
                if (SumEpsilon < Epsilon)
                {
                    Out_Flag = true;
                }

                if (!Out_Flag)
                {
                    // 计算反射点
                    for (int i = 0; i < X_Dim; i++)
                    {
                        X_R[i] = X_Center[i] + Alpha * (X_Center[i] - X[Index_Max, i]);
                        X_R[i] = G.BoundConstraint(X_R[i], this.X_Lb[i], this.X_Ub[i]);
                    }
                    Fun_R = ObjectFunction(X_R, out Delta_R);
                    FEs = FEs + 1;

                    if ((Delta_R == Delta_Min && Fun_R < Fun_Min) || Delta_R < Delta_Min)
                    {
                        // 计算扩大点
                        for (int i = 0; i < X_Dim; i++)
                        {
                            X_E[i] = X_Center[i] + Gamma * (X_R[i] - X_Center[i]);
                            X_E[i] = G.BoundConstraint(X_E[i], this.X_Lb[i], this.X_Ub[i]);
                        }
                        Fun_E = ObjectFunction(X_E, out Delta_E);
                        FEs = FEs + 1;
                        if ((Delta_E == Delta_R && Fun_E < Fun_R) || Delta_E < Delta_R)
                        {
                            // 扩大点替换最大点
                            Fun[Index_Max] = Fun_E;
                            Delta[Index_Max] = Delta_E;
                            Fun_Max = Fun_E;
                            Delta_Max = Delta_E;
                            for (int i = 0; i < X_Dim; i++)
                            {
                                X[Index_Max, i] = X_E[i];
                            }
                        }
                        else
                        {
                            // 反射点替换最大点
                            Fun[Index_Max] = Fun_R;
                            Delta[Index_Max] = Delta_R;
                            Fun_Max = Fun_R;
                            Delta_Max = Delta_R;
                            for (int i = 0; i < X_Dim; i++)
                            {
                                X[Index_Max, i] = X_R[i];
                            }
                        }
                    }
                    else if ((Delta_R == Delta_SecMax && Fun_R < Fun_SecMax) || Delta_R < Delta_SecMax)
                    {
                        // 反射点替换最大点
                        Fun[Index_Max] = Fun_R;
                        Delta[Index_Max] = Delta_R;
                        Fun_Max = Fun_R;
                        Delta_Max = Delta_R;
                        for (int i = 0; i < X_Dim; i++)
                        {
                            X[Index_Max, i] = X_R[i];
                        }
                    }
                    else if ((Delta_R == Delta_Max && Fun_R < Fun_Max) || Delta_R < Delta_Max)
                    {
                        // 反射点替换最大点
                        Fun[Index_Max] = Fun_R;
                        Delta[Index_Max] = Delta_R;
                        Fun_Max = Fun_R;
                        Delta_Max = Delta_R;
                        for (int i = 0; i < X_Dim; i++)
                        {
                            X[Index_Max, i] = X_R[i];
                        }
                        // 计算收缩点
                        for (int i = 0; i < X_Dim; i++)
                        {
                            X_C[i] = X_Center[i] + Beta * (X[Index_Max, i] - X_Center[i]);
                            X_C[i] = G.BoundConstraint(X_C[i], this.X_Lb[i], this.X_Ub[i]);
                        }
                        Fun_C = ObjectFunction(X_C, out Delta_C);
                        FEs = FEs + 1;
                        if ((Delta_C == Delta_Max && Fun_C < Fun_Max) || Delta_C < Delta_Max)
                        {
                            // 收缩点替换最大点
                            Fun[Index_Max] = Fun_C;
                            Delta[Index_Max] = Delta_C;
                            Fun_Max = Fun_C;
                            Delta_Max = Delta_C;
                            for (int i = 0; i < X_Dim; i++)
                            {
                                X[Index_Max, i] = X_C[i];
                            }
                        }
                        else
                        {
                            // 缩边
                            for (int i = 0; i < X_Dim + 1; i++)
                            {
                                if (i != Index_Min)
                                {
                                    for (int j = 0; j < X_Dim; j++)
                                    {
                                        X[i, j] = X[Index_Min, j] + 0.5 * (X[i, j] - X[Index_Min, j]);
                                        X[i, j] = G.BoundConstraint(X[i, j], this.X_Lb[j], this.X_Ub[j]);
                                        X_tmp[j] = X[i, j];
                                    }
                                    Fun[i] = ObjectFunction(X_tmp, out Delta[i]);
                                    FEs = FEs + 1;
                                }
                            }
                        }
                    }
                    else
                    {
                        // 计算收缩点
                        for (int i = 0; i < X_Dim; i++)
                        {
                            X_C[i] = X_Center[i] + Beta * (X[Index_Max, i] - X_Center[i]);
                            X_C[i] = G.BoundConstraint(X_C[i], this.X_Lb[i], this.X_Ub[i]);
                        }
                        Fun_C = ObjectFunction(X_C, out Delta_C);
                        FEs = FEs + 1;
                        if ((Delta_C == Delta_Max && Fun_C < Fun_Max) || Delta_C < Delta_Max)
                        {
                            // 收缩点替换最大点
                            Fun[Index_Max] = Fun_C;
                            Delta[Index_Max] = Delta_C;
                            Fun_Max = Fun_C;
                            Delta_Max = Delta_C;
                            for (int i = 0; i < X_Dim; i++)
                            {
                                X[Index_Max, i] = X_C[i];
                            }
                        }
                        else
                        {
                            // 缩边
                            for (int i = 0; i < X_Dim + 1; i++)
                            {
                                if (i != Index_Min)
                                {
                                    for (int j = 0; j < X_Dim; j++)
                                    {
                                        X[i, j] = X[Index_Min, j] + 0.5 * (X[i, j] - X[Index_Min, j]);
                                        X[i, j] = G.BoundConstraint(X[i, j], this.X_Lb[j], this.X_Ub[j]);
                                        X_tmp[j] = X[i, j];
                                    }
                                    Fun[i] = ObjectFunction(X_tmp, out Delta[i]);
                                    FEs = FEs + 1;
                                }
                            }
                        }
                    }
                }

                // 得到最优适应度 
                for (int i = 0; i < X_Dim + 1; i++)
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

                if (iter > 0)
                {
                    TolFun = Math.Abs((Last_Best_Fun - Best_Fun) / Last_Best_Fun); // 最优适应度相对变化
                }
                iter = iter + 1;
                Last_Best_Fun = Best_Fun; // 上一步目标值
                //if (TolFun < Min_TolFun)
                //{
                //    break;
                //}
                X_Matrix = M.DenseOfArray(X).Transpose();  // 解矩阵
                X_Rank = X_Matrix.Rank(); // 解矩阵的秩
                Console.WriteLine("Nelder_Mead " + FEs + " " + Best_Fun.ToString("f2") + " " + Best_Delta.ToString("f2") + " " + X_Rank.ToString());
            }
        }

        /// <summary>初始化</summary>
        private void Initial()
        {
            double[] X_tmp = new double[X_Dim];
            for (int i = 0; i < X_Dim + 1; i++)
            {
                for (int j = 0; j < X_Dim; j++)
                {
                    X_tmp[j] = X[i, j];
                }

                Fun[i] = ObjectFunction(X_tmp, out Delta[i]); ;
            }
            Best_Fun = Fun[0];
            Best_Delta = Delta[0];
            for (int i = 0; i < X_Dim + 1; i++)
            {
                if ((Delta[i] == Best_Delta && Fun[i] < Best_Fun) || Delta[i] < Best_Delta)
                {
                    Best_Fun = Fun[i];
                    Best_Delta = Delta[i];
                    Best_Index = i;
                }
            }
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