using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InterfaceDeclare;
using System.IO;

namespace FitnessSet
{
    public class RevisedCOT:IFitness
    {
        private double dblOutputArg;

        public void SetCoefficient(double[] dblVal, int[] nProductID = null)
        {
            return;
        }

        public void SetOutput(double[] dblVal)
        {
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            dblOutputArg = dblVal[0];
        }

        public void SetConstraint(double[] dblConstraint)
        {
            return;
        }

        public double Fitness(double[] dblInputArg, IModel[] model, out double dblViolation)
        {
            try
            {
                model[0].SetInput(dblInputArg);
                model[0].ModelCalculate();
                double[] dblOutput = model[0].GetOutput();
                dblViolation = 0.0;
                return -1.0 * Math.Abs(dblOutput[0] - dblOutputArg);
            }
            catch (System.Exception ex)
            {
                throw new Exception(ex.Message);
            }
        }
    }

    public class ObjectFunc:IFitness
    {
        private double[] dblCoefficient;//价格因子
        private int[] nIndex;//价格因子对应的产品ID

        public void SetCoefficient(double[] dblVal, int[] nProductID = null)
        {
            if (dblVal == null)
            {
                throw new Exception("Input Coefficient is null!");
            }
            if (nProductID != null)
            {
                if (nProductID.Length != dblVal.Length)
                {
                    throw new Exception("Input ProductID's dimension is not match with Coefficient!");
                }
            }
            if (nProductID == null)
            {
                throw new Exception("Input ProductID is null!");
            }
            dblCoefficient = new double[dblVal.Length];
            nIndex = new int[nProductID.Length];
            for (int i = 0; i < dblVal.Length; i++ )
            {
                dblCoefficient[i] = dblVal[i];
                nIndex[i] = nProductID[i];
            }
        }

        public void SetOutput(double[] dblVal)
        {
            return;
        }

        public void SetConstraint(double[] dblConstraint)
        {
            return;
        }

        public double Fitness(double[] dblInputArg, IModel[] model, out double dblViolation)
        {
            //模型0为coilsim收率模型，1为燃料气模型，2为SS蒸汽模型
            //输入参数
            //0----COT
            //1----DHR
            //2----打靶目标序号
            //3----Feed
            //4----COP
            //5----炉管数
            //6----排烟温度
            //7----氧含量
            try
            {
                double[] dblInput0 = new double[6];
                for (int i = 0; i < 6; i++ )
                {
                    dblInput0[i] = dblInputArg[i];
                }
                model[0].SetInput(dblInput0);
                model[0].ModelCalculate();
                double[] dblOutput0 = model[0].GetOutput();
                int[] nProductID0 = model[0].GetProductID();
                Dictionary<int, double> dic = new Dictionary<int, double>();
                for (int i = 0; i < nProductID0.Length; i++ )
                {
                    dic.Add(nProductID0[i], dblOutput0[i] * dblInputArg[3]);//收率乘以进料量
                }
                dic.Add(-1, dblInputArg[3]);//进料量
                double[] dblInput1 = new double[5];
                dblInput1[0] = dblInputArg[0];
                dblInput1[1] = dblInputArg[1];
                dblInput1[2] = dblInputArg[3];
                dblInput1[3] = dblInputArg[6];
                dblInput1[4] = dblInputArg[7];
                model[1].SetInput(dblInput1);
                model[1].ModelCalculate();
                double[] dblOutput1 = model[1].GetOutput();
                dic.Add(-2, dblOutput1[0]);//燃料气量
                dic.Add(-3, dblInputArg[3] * dblInputArg[1]);//DS量
                double[] dblInput2 = new double[3];
                dblInput2[0] = dblInputArg[0];
                dblInput2[1] = dblInputArg[1];
                dblInput2[2] = dblInputArg[2];
                model[2].SetInput(dblInput2);
                model[2].ModelCalculate();
                double[] dblOutput2 = model[2].GetOutput();
                dic.Add(1000, dblOutput2[0]);//SS蒸汽量
                double dblResult = 0.0;
                if (dblCoefficient == null)
                {
                    throw new Exception("Coefficient is not set!");
                }
                if (nIndex == null)
                {
                    throw new Exception("ProductID is not set!");
                }
                for (int i = 0; i < nIndex.Length; i++)
                {
                    if (dic.ContainsKey(nIndex[i]))
                    {
                        dblResult += dic[nIndex[i]] * dblCoefficient[i];
                    }
                }
                dblViolation = 0.0;
                return dblResult;
            }
            catch (System.Exception ex)
            {
                throw new Exception(ex.Message);
            }
        }
    }

    public class AspenCorrection : IFitness
    {
        private double[] dblOutputArg;

        public void SetCoefficient(double[] dblVal, int[] nProductID = null)
        {
            return;
        }

        public void SetOutput(double[] dblVal)
        {
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            dblOutputArg = new double[dblVal.Length];
            for (int i = 0; i < dblVal.Length; i++)
            {
                dblOutputArg[i] = dblVal[i];
            }            
        }

        public void SetConstraint(double[] dblConstraint)
        {
            return;
        }

        public double Fitness(double[] dblInputArg, IModel[] model, out double dblViolation)
        {
            try
            {
                double dblFitness = 0.0;
                model[0].SetInput(dblInputArg);
                model[0].ModelCalculate();
                double[] dblOutput = model[0].GetOutput();
                dblViolation = 0.0;
                for (int i = 0; i < dblOutput.Length; i++)
                {
                    dblFitness = dblFitness + Math.Pow(dblOutput[i] - dblOutputArg[i], 2);
                }
                return dblFitness;
            }
            catch (System.Exception ex)
            {
                dblViolation = double.MaxValue;
                return double.MaxValue;
            }
        }
    }

    public class AspenMASSFLOW : IFitness
    {
        private double dblCon;

        public void SetCoefficient(double[] dblVal, int[] nProductID = null)
        {
            return;
        }

        public void SetOutput(double[] dblVal)
        {
            return;
        }

        public void SetConstraint(double[] dblConstraint)
        {
            dblCon = dblConstraint[0];
            return;
        }

        public double Fitness(double[] dblInputArg, IModel[] model, out double dblViolation)
        {
            try
            {
                double dblFitness = 0.0;
                dblViolation = 0.0;
                model[0].SetInput(dblInputArg);
                model[0].ModelCalculate();
                double[] dblOutput = model[0].GetOutput();

                dblFitness = dblFitness + dblOutput[0] * -1.0;
                dblViolation = dblCon - dblOutput[1] + dblOutput[2] + dblOutput[3];
                dblViolation = dblViolation < 0.0 ? dblViolation : 0.0;
                return dblFitness;
            }
            catch (System.Exception ex)
            {
                dblViolation = double.MaxValue;
                return double.MaxValue;
            }
        }
    }
}