using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using InterfaceDeclare;

namespace OptimizerSet
{
    public class GoldenSection
    {
        #region 私有变量
        private int MaX_Iteration = 300; // 最大迭代次数
        private double Min_TolFun = 1.0e-4; // 结束条件
        private int X_Dim = 0; // 变量维度
        private double[] X_Lb; // 优化变量下限
        private double[] X_Ub; // 优化变量上限
        private double[] X_Opt; // 最优解
        private double Best_Fun = 0.0;  // 最优目标函数值
        private double Best_Delta = 0.0;  // 最优约束违反度
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
                this.MaX_Iteration = int.Parse(Config[0][1]); // 最大迭代次数
                this.Min_TolFun = double.Parse(Config[1][1]); // 结束条件
            }
            catch (System.Exception ex)
            {
                throw new Exception("ConfigFile are not existed!");
            }
            this.X_Dim = X_Init.Count();
            this.X_Lb = new double[X_Dim];
            this.X_Ub = new double[X_Dim];
            this.X_Opt = new double[X_Dim];
            for (int i = 0; i < X_Dim; i++)
            {
                this.X_Lb[i] = Math.Max(X_Init[i] - X_MaxStep[i], X_Lb[i]);
                this.X_Ub[i] = Math.Min(X_Init[i] + X_MaxStep[i], X_Ub[i]);
            }
            Opt(); // 优化
            Fun = Best_Fun;
            Delta = Best_Delta;
            return X_Opt;
        }

        /// <summary>
        /// 黄金分割
        /// </summary>
        private void Opt()
        {
            int Count = 0; // 当前迭代次数
            double GoldenPoint = (Math.Sqrt(5.0) - 1.0) / 2.0;
            double Fun_Low_Try = 0.0;
            double Fun_High_Try = 0.0;
            double Delta_Low_Try = 0.0;
            double Delta_High_Try = 0.0;
            double[] X_Range = new double[X_Dim];
            double[] X_Low_Try = new double[X_Dim]; // 右试探点
            double[] X_High_Try = new double[X_Dim]; // 左试探点

            X_Range[0] = X_Ub[0] - X_Lb[0];
            X_Low_Try[0] = X_Lb[0] + (1.0 - GoldenPoint) * X_Range[0];
            X_High_Try[0] = X_Lb[0] + GoldenPoint * X_Range[0];

            while (Count < MaX_Iteration && Math.Abs(X_Range[0]) > Min_TolFun)
            {
                Fun_Low_Try = ObjectFunction(X_Low_Try, out Delta_Low_Try);
                Fun_High_Try = ObjectFunction(X_High_Try, out Delta_High_Try);
                if ((Delta_Low_Try == 0.0 && Delta_High_Try == 0.0 && Fun_Low_Try < Fun_High_Try) 
                    || (Delta_Low_Try > 0.0 && Delta_High_Try > 0.0 && Delta_Low_Try < Delta_High_Try)
                    || (Delta_Low_Try > 0.0 && Delta_High_Try > 0.0 && Delta_Low_Try == Delta_High_Try && Fun_Low_Try < Fun_High_Try))
                {
                    X_Ub[0] = X_High_Try[0]; // 更新右端点，左端点不动
                    X_Range[0] = X_Ub[0] - X_Lb[0]; // 更新x的搜索范围
                    X_Low_Try[0] = X_Lb[0] + (1.0 - GoldenPoint) * X_Range[0]; // 更新试探点
                    X_High_Try[0] = X_Lb[0] + GoldenPoint * X_Range[0]; // 更新试探点
                }
                else if (((Delta_Low_Try + Delta_High_Try == 0.0) && Fun_Low_Try >= Fun_High_Try) || ((Delta_Low_Try + Delta_High_Try != 0.0) && Delta_Low_Try >= Delta_High_Try))
                {
                    X_Lb[0] = X_Low_Try[0]; // 更新右端点，左端点不动
                    X_Range[0] = X_Ub[0] - X_Lb[0]; // 更新x的搜索范围
                    X_Low_Try[0] = X_Lb[0] + (1.0 - GoldenPoint) * X_Range[0]; // 更新试探点
                    X_High_Try[0] = X_Lb[0] + GoldenPoint * X_Range[0]; // 更新试探点
                }
                Count++;
                Console.WriteLine("GoldenSection Iter: " + Count + " " + System.DateTime.Now.ToString("yy-MM-dd hh:mm:ss"));
                Console.WriteLine("Function Low: " + X_Low_Try[0].ToString("f2") + " Function Low: " + Fun_Low_Try + " Delta Low: " + Delta_Low_Try);
                Console.WriteLine("Function High: " + X_High_Try[0].ToString("f2") + " Function High: " + Fun_High_Try + " Delta High: " + Delta_High_Try);
            }
            // 得到最优解
            X_Opt[0] = (X_Lb[0] + X_Ub[0]) / 2;
            Best_Fun = ObjectFunction(X_Opt, out Best_Delta);
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
        }
    }
}
