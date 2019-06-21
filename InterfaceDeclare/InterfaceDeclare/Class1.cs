using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace InterfaceDeclare//声明接口
{
    public interface IModel
    {
        /// <summary>
        /// 设置模型参数
        /// </summary>
        /// <param name="strVal">参数路径</param>
        /// <returns>成功返回true，失败返回false</returns>
        bool SetPara(string strVal);

        /// <summary>
        /// 判断模型参数是否设置完成
        /// </summary>
        /// <returns>完成则返回true，未完成返回false</returns>
        bool HasParaSet();

        /// <summary>
        /// 设置模型输入
        /// </summary>
        /// <param name="dblVal">模型输入</param>
        void SetInput(double[] dblVal);

        /// <summary>
        /// 模型计算
        /// </summary>
        void ModelCalculate();

        /// <summary>
        /// 获得模型输出
        /// </summary>
        /// <returns></returns>
        double[] GetOutput();

        /// <summary>
        /// 获得产品ID
        /// </summary>
        /// <returns></returns>
        int[] GetProductID();

        /// <summary>
        /// 关闭模型
        /// </summary>
        /// <returns></returns>
        void ModelClose();
    }

    /// <summary>
    /// 适应度函数
    /// </summary>
    public interface IFitness
    {
        /// <summary>
        /// 设置价格因子
        /// </summary>
        /// <param name="dblVal">价格因子</param>
        /// <param name="nProductID">价格因子对应的产品ID</param>
        void SetCoefficient(double[] dblVal, int[] nProductID = null);

        /// <summary>
        /// 设置模型逼近的输出
        /// </summary>
        /// <param name="dblVal">模型逼近的输出</param>
        void SetOutput(double[] dblVal);

        /// <summary>
        /// 设置约束参数
        /// </summary>
        /// <param name="dblVal">约束参数</param>
        void SetConstraint(double[] dblConstraint);

        /// <summary>
        /// 适应度函数
        /// </summary>
        /// <param name="dblInputArg">函数输入</param>
        /// <param name="model">模型接口</param>
        /// <param name="dblViolation">违反度，输出参数，未违反约束时为0</param>
        /// <returns>适应度函数值</returns>
        double Fitness(double[] dblInputArg, IModel[] model, out double dblViolation);
    }

    public interface IMultiFitness
    {
        /// <summary>
        /// 设置价格因子
        /// </summary>
        /// <param name="dblVal">价格因子</param>
        /// <param name="nProductID">价格因子对应的产品ID</param>
        void SetCoefficient(double[] dblVal, int[] nProductID = null);

        /// <summary>
        /// 设置模型逼近的输出
        /// </summary>
        /// <param name="dblVal">模型逼近的输出</param>
        void SetOutput(double[] dblVal);

        /// <summary>
        /// 设置约束参数
        /// </summary>
        /// <param name="dblVal">约束参数</param>
        void SetConstraint(double[] dblConstraint);

        /// <summary>
        /// 获得目标个数
        /// </summary>
        int GetObjectiveNumber();

        /// <summary>
        /// 适应度函数
        /// </summary>
        /// <param name="dblInputArg">函数输入</param>
        /// <param name="model">模型接口</param>
        /// <param name="dblViolation">违反度，输出参数，未违反约束时为0</param>
        /// <returns>适应度函数值</returns>
        double[] Fitness(double[] dblInputArg, IModel[] model, out double dblViolation);
    }
}
