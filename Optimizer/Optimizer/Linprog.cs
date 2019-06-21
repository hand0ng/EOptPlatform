using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using lpsolve55;

namespace OptimizerSet
{
    public class Linprog
    {
        /// <summary>
        /// 开始优化
        /// </summary>
        /// <param name="f">目标函数系数向量</param>
        /// <param name="LB">优化变量下限</param>
        /// <param name="UB">优化变量上限</param>
        /// <returns><param name="Fun">最优目标函数</param></returns>
        /// <returns><param name="Flag">求解标志</param></returns>
        /// <param name="A">不等式约束左边系数，可选参数</param>
        /// <param name="b">不等式约束右边系数，可选参数</param>
        /// <summary>A*X<=b</summary>
        /// <param name="Aeq">等式约束左边系数，可选参数</param>
        /// <param name="Aeq">等式约束右边系数，可选参数</param>
        /// <summary>Aeq*X=beq</summary>
        /// <param name="isInteger">是否为整数，可选参数</param>
        /// <param name="isBinary">是否为01整数，可选参数</param>
        /// <returns>使适应度函数最小的变量值</returns>
        /// </summary>
        public double[] StartOpt(double[] f, double[] LB, double[] UB, out double Fun, out int Flag, double[,] A = null, double[] b = null, double[,] Aeq = null, double[] beq = null, int[] isInteger = null, int[] isBinary = null)
        {
            if (f.Length != LB.Length || f.Length != UB.Length || (A != null && f.Length != A.GetLength(1)) || (Aeq != null && f.Length != Aeq.GetLength(1)))
            {
                 throw new Exception("Variable number are not set correctly!");
            }
            if (A.GetLength(0) != b.Length)
            {
                throw new Exception("Inequality constraints are not set correctly!");
            }
            if (Aeq.GetLength(0) != beq.Length)
            {
                throw new Exception("Equality constraints are not set correctly!");
            }
            if (isInteger != null)
            {
                if (f.Length != isInteger.Length)
                {
                    throw new Exception("Integer variable flags are not set correctly!");
                }
                for (int i = 0; i < isInteger.Length; i++)
                {
                    if (isInteger[i] != 0 && isInteger[i] != 1)
                    {
                        throw new Exception("Integer variable flags are not set correctly!");
                    }
                }
            }
            if (isBinary != null)
            {
                if (f.Length != isBinary.Length)
                {
                    throw new Exception("Binary variable flags are not set correctly!");
                }
                for (int i = 0; i < isBinary.Length; i++)
                {
                    if (isBinary[i] != 0 && isBinary[i] != 1)
                    {
                        throw new Exception("Binary variable flags are not set correctly!");
                    }
                }       
            }
            IntPtr lp;
            int n = f.Length; // 变量数
            int m = A.GetLength(0) + Aeq.GetLength(0); // 约束数
            int m1 = A.GetLength(0); // 不等式约束数
            int m2 = Aeq.GetLength(0); // 等式约束数
            double[] Row = new double[n + 1];
            double[] X_Opt;
            Fun = 0.0;
            lp = lpsolve.make_lp(0, n); // n个变量
            // 不等式约束
            if (A != null)
            {
                for (int i = 0; i < m1; i++)
                {
                    for (int j = 0; j < n; j++)
                    {

                        Row[j + 1] = A[i, j];
                    }
                    Row[0] = 0;
                    lpsolve.add_constraint(lp, Row, lpsolve.lpsolve_constr_types.LE, b[i]);
                }
            }
            // 等式约束
            if (Aeq != null)
            {
                for (int i = 0; i < m2; i++)
                {
                    for (int j = 0; j < n; j++)
                    {

                        Row[j + 1] = Aeq[i, j];
                    }
                    Row[0] = 0;
                    lpsolve.add_constraint(lp, Row, lpsolve.lpsolve_constr_types.EQ, beq[i]);
                }
            }
            // 设置上下限
            for (int i = 0; i < n; i++)
            {
                lpsolve.set_lowbo(lp, i + 1, LB[i]);
            }
            for (int i = 0; i < n; i++)
            {

                lpsolve.set_upbo(lp, i + 1, UB[i]);
            }
            // 目标函数
            for (int i = 0; i < n; i++)
            {

                Row[i + 1] = f[i];
            }
            Row[0] = 0;
            lpsolve.set_obj_fn(lp, Row);
            // 整数变量
            if (isInteger != null)
            {
                for (int i = 0; i < isInteger.Length; i++)
                {
                    if (isInteger[i] == 1)
                    {
                        lpsolve.set_int(lp, i, true);  // sets variable 1 to integer
                    }
                }
            }
            // 01整数变量
            if (isBinary != null)
            {
                for (int i = 0; i < isBinary.Length; i++)
                {
                    if (isBinary[i] == 1)
                    {
                        lpsolve.set_binary(lp, i, true); // sets variable 1 to integer
                    }
                }
            }
            // 求解
            lpsolve.solve(lp);
            
            // 结果输出
            X_Opt = new double[lpsolve.get_Ncolumns(lp)];
            lpsolve.get_variables(lp, X_Opt);
            Fun = lpsolve.get_objective(lp);
            Flag = lpsolve.get_status(lp); // 问题是否有解

            lpsolve.delete_lp(lp); // 删除定义的pl问题
            return X_Opt;
        }
    }
}