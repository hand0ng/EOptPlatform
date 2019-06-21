using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InterfaceDeclare;
using System.IO;
using System.Diagnostics;
using System.Threading;
using System.ServiceProcess;
using AspenCustomModelerLibrary;

namespace ModelSet
{
    public class BPNN : IModel
    {
        private bool blErr = false;//设置参数错误标志
        private int nInputNum;//输入节点个数
        private int nMediumNum;//中间节点个数
        private int nOutputNum;//输出节点个数
        private int nMidFuncType;//中间层激活函数
        private int nOutFuncType;//输出层激活函数
        private int nInputNormalType;//输入归一化方式
        private int nOutputNormalType;//输出反归一化方式
        private double[] dblInput;//输入
        private double[] dblOutput;//输出
        private double[] dblInputW;//输入层到中间层权值
        private double[] dblOutputW;//中间层到输出层权值
        private double[] dblMediumB;//中间层偏置
        private double[] dblOutputB;//输出层偏置
        public double[] dblInputH;//输入上限
        public double[] dblInputL;//输入下限
        private double[] dblOutputH;//输出上限
        private double[] dblOutputL;//输出下限
        public string strParaPath;//配置文件路径


        #region 激活函数

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double tansig(double x)
        {
            double y = 0.0;

            y = 2.0 / (1.0 + Math.Exp(-2.0 * x)) - 1.0;

            return y;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double purelin(double x)
        {
            double y = 0.0;
            y = x;

            return y;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double logsig(double x)
        {
            double y = 0.0;
            y = 1.0 / (1.0 + Math.Exp(-x));

            return y;
        }
        #endregion


        public double[] GetOutput()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblOutput == null)
            {
                throw new Exception("Model Output is null!");
            }
            double[] dblResult = new double[nOutputNum];
            for (int i = 0; i < nOutputNum; i++ )
            {
                dblResult[i] = dblOutput[i];
            }
            return dblResult;
        }

        public void SetInput(double[] dblVal)
        {
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            if (dblVal.Length != nInputNum)
            {
                throw new Exception("Input Value's dimension is not match with the model!");
            }
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            for (int i = 0; i < nInputNum; i++ )
            {
                dblInput[i] = dblVal[i];
            }
        }

        public bool SetPara(string strVal)
        {
            if (!File.Exists(strVal))
            {
                blErr = true;
                throw new Exception("Input File Directory is not valid!");
            }
            List<string[]> lsPara = CSVOperation.ReadCSV(strVal);
            if (lsPara.Count <= 0)
            {
                blErr = true;
                throw new Exception("Input File is empty!");
            }
            nInputNum = Convert.ToInt32(lsPara[1][0]);
            nMediumNum = Convert.ToInt32(lsPara[1][1]);
            nOutputNum = Convert.ToInt32(lsPara[1][2]);
            nMidFuncType = Convert.ToInt32(lsPara[1][3]);
            nOutFuncType = Convert.ToInt32(lsPara[1][4]);
            nInputNormalType = Convert.ToInt32(lsPara[1][5]);
            nOutputNormalType = Convert.ToInt32(lsPara[1][6]);

            dblInput = new double[nInputNum];
            dblOutput = new double[nOutputNum];
            dblInputW = new double[nInputNum * nMediumNum];
            dblOutputW = new double[nMediumNum * nOutputNum];
            dblMediumB = new double[nMediumNum];
            dblOutputB = new double[nOutputNum];
            dblInputH = new double[nInputNum];
            dblInputL = new double[nInputNum];
            dblOutputH = new double[nOutputNum];
            dblOutputL = new double[nOutputNum];

            for (int i = 0; i < nMediumNum; i++)
            {
                for (int j = 0; j < nInputNum; j++)
                {
                    dblInputW[i * nInputNum + j] = Convert.ToDouble(lsPara[i + 3][j + 1]);
                }
                dblMediumB[i] = Convert.ToDouble(lsPara[i + 3][nInputNum + 1]);
            }
            for (int i = 0; i < nOutputNum; i++)
            {
                for (int j = 0; j < nMediumNum; j++)
                {
                    dblOutputW[i * nMediumNum + j] = Convert.ToDouble(lsPara[j + 3][nInputNum + i + 2]);
                }
                dblOutputB[i] = Convert.ToDouble(lsPara[nMediumNum + 3][nInputNum + i + 2]);

                dblOutputH[i] = Convert.ToDouble(lsPara[nMediumNum + nInputNum + 5 + i][1]);
                dblOutputL[i] = Convert.ToDouble(lsPara[nMediumNum + nInputNum + 5 + i][2]);
            }

            for (int i = 0; i < nInputNum; i++)
            {
                dblInputH[i] = Convert.ToDouble(lsPara[nMediumNum + 5 + i][1]);
                dblInputL[i] = Convert.ToDouble(lsPara[nMediumNum + 5 + i][2]);
            }
            blErr = false;
            strParaPath = strVal;
            return true;
        }

        public bool HasParaSet()
        {
            return !blErr;
        }

        public void ModelCalculate()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            //输入归一化
            double[] dblNormalInput = new double[nInputNum];
            for (int i = 0; i < nInputNum; i++)
            {
                double dblInputVal = dblInput[i];
                double dblInputHVal = dblInputH[i];
                double dblInputLVal = dblInputL[i];
                if (dblInputVal > dblInputHVal)
                {
                    dblInputVal = dblInputHVal;
                }
                if (dblInputVal < dblInputLVal)
                {
                    dblInputVal = dblInputLVal;
                }
                switch (nInputNormalType)
                {
                    case 0:
                        {
                            dblNormalInput[i] = (dblInputVal - dblInputLVal) / (dblInputHVal - dblInputLVal);
                        }
                        break;
                    case 1:
                        {
                            dblNormalInput[i] = 2.0 * (dblInputVal - dblInputLVal) / (dblInputHVal - dblInputLVal) - 1.0;
                        }
                        break;
                    default:
                        {
                            dblNormalInput[i] = (dblInputVal - dblInputLVal) / (dblInputHVal - dblInputLVal);
                        }
                        break;
                }
            }

            //计算中间节点
            double[] dblMedium = new double[nMediumNum];
            for (int i = 0; i < nMediumNum; i++)
            {
                double mediumInput = 0;
                for (int j = 0; j < nInputNum; j++)
                {
                    mediumInput += dblInputW[j + nInputNum * i] * dblNormalInput[j];
                }
                mediumInput = mediumInput + dblMediumB[i];
                switch (nMidFuncType)
                {
                    case 0:
                        dblMedium[i] = tansig(mediumInput);
                        break;
                    case 1:
                        dblMedium[i] = purelin(mediumInput);
                        break;
                    case 2:
                        dblMedium[i] = logsig(mediumInput);
                        break;
                    default:
                        dblMedium[i] = tansig(mediumInput);
                        break;
                }
            }

            //计算输出层节点
            double dblOutputVal = 0.0;
            for (int i = 0; i < nOutputNum; i++)
            {
                double dblOutputHVal = dblOutputH[i];
                double dblOutputLVal = dblOutputL[i];
                double opLayerInput = 0;
                for (int j = 0; j < nMediumNum; j++)
                {
                    opLayerInput += dblOutputW[j + nMediumNum * i] * dblMedium[j];
                }
                opLayerInput = opLayerInput + dblOutputB[i];

                switch (nOutFuncType)
                {
                    case 0:
                        dblOutputVal = tansig(opLayerInput);
                        break;
                    case 1:
                        dblOutputVal = purelin(opLayerInput);
                        break;
                    case 2:
                        dblOutputVal = logsig(opLayerInput);
                        break;
                    default:
                        dblOutputVal = tansig(opLayerInput);
                        break;
                }

                //输出反归一化
                double dblFinalOutput = 0.0;
                switch (nOutputNormalType)
                {
                    case 0:
                        {
                            dblFinalOutput = dblOutputVal * (dblOutputHVal - dblOutputLVal) + dblOutputLVal;
                        }
                        break;
                    case 1:
                        {
                            dblFinalOutput = (dblOutputVal + 1.0) * (dblOutputHVal - dblOutputLVal) / 2.0 + dblOutputLVal;
                        }
                        break;
                    default:
                        {
                            dblFinalOutput = dblOutputVal * (dblOutputHVal - dblOutputLVal) + dblOutputLVal;
                        }
                        break;
                }
                if (dblFinalOutput > dblOutputHVal)
                {
                    dblFinalOutput = dblOutputHVal;
                }
                if (dblFinalOutput < dblOutputLVal)
                {
                    dblFinalOutput = dblOutputLVal;
                }
                dblOutput[i] = (double)dblFinalOutput;
            }
        }

        public int[] GetProductID()
        {
            return null;
        }

        public void ModelClose()
        {
            return;
        }
    }

    public class Coilsim : IModel
    {
        private double[] dblInput;//输入
        private double[] dblOutput;//输出
        private int[] nProductID;//产品ID
        private string strExePath;//exe文件路径
        private string strCasePath;//case文件路径
        private string strCaseName;//case名
        private bool blErr = false;//设置参数错误标志

        public double[] GetOutput()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblOutput == null)
            {
                throw new Exception("Model Output is null!");
            }
            double[] dblResult = new double[dblOutput.Length];
            for (int i = 0; i < dblOutput.Length; i++)
            {
                dblResult[i] = dblOutput[i];
            }
            return dblResult;
        }

        public void SetInput(double[] dblVal)
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            if (dblVal.Length != 6)
            {
                throw new Exception("Input Value's dimension is not match with model!");
            }
            dblInput = new double[dblVal.Length];
            for (int i = 0; i < dblVal.Length; i++ )
            {
                dblInput[i] = dblVal[i];
            }
        }

        public bool SetPara(string strVal)
        {
            if (!Directory.Exists(strVal))
            {
                blErr = true;
                throw new Exception("Input Directory is not valid!");
            }
            strCasePath = strVal;
            string[] strSlot = strVal.Split('\\');
            string[] strTubeSlot = new string[strSlot.Length - 4];
            string[] strExeSlot = new string[strSlot.Length - 2];
            for (int i = 0; i < strSlot.Length - 4; i++)
            {
                strTubeSlot[i] = strSlot[i];
            }
            for (int i = 0; i < strSlot.Length - 2; i++)
            {
                strExeSlot[i] = strSlot[i];
            }
            strExePath = String.Join("\\", strExeSlot);
            strCaseName = strSlot[strSlot.Length - 1];
            if (!Directory.Exists(strExePath))
            {
                blErr = true;
                throw new Exception("Input Directory has no coilsim.exe Path!");
            }
            blErr = false;
            return true;
        }

        public bool HasParaSet()
        {
            return !blErr;
        }

        public void ModelCalculate()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            //获得炉管数
            int nTubeNum = (int)dblInput[5];

            //设置操作条件
            List<string[]> lsExp = CSVOperation.ReadCSV(strCasePath + "\\exp.txt");
            lsExp[1][0] = ((int)dblInput[2]).ToString();//打靶目标序号
            lsExp[2][0] = dblInput[0].ToString();//打靶目标值
            lsExp[6][0] = dblInput[4].ToString();//COP
            lsExp[8][0] = (dblInput[3] * 1000.0 / (double)nTubeNum).ToString();//进料量
            lsExp[9][0] = dblInput[1].ToString();//汽烃比
            CSVOperation.WriteCSV(strCasePath + "\\exp.txt", false, lsExp);

            //设置Simulation文件
            //List<string[]> lsSimulation = CSVOperation.ReadCSV(strExePath + "\\Projects\\Simulation.txt");
            //lsSimulation[2][0] = strCaseName;
            List<string[]> lsSimulation = new List<string[]>();
            lsSimulation.Add(new string[1] { "1" });
            lsSimulation.Add(new string[1] { "0" });
            lsSimulation.Add(new string[1] { strCaseName });
            CSVOperation.WriteCSV(strExePath + "\\Projects\\Simulation.txt", false, lsSimulation);

            //调用coilsim
            Process process = new Process();
            process.StartInfo.FileName = strExePath + "\\Coilsim.exe";
            process.StartInfo.WindowStyle = ProcessWindowStyle.Normal;
            process.StartInfo.WorkingDirectory = strExePath;
            process.Start();
            Thread taskKill = new Thread(new ParameterizedThreadStart(KillCoilsim));
            taskKill.Start(process);
            process.WaitForExit();
             

            //获得产品收率
            if (!File.Exists(strCasePath + "\\yields.csv"))
            {
                nProductID = null;
                dblOutput = null;
                throw new Exception("Model Failed to Converge!");
            }
            List<string[]> lsYields = CSVOperation.ReadCSV(strCasePath + "\\yields.csv");
            int nProductNum = lsYields.Count;
            nProductID = new int[nProductNum];
            dblOutput = new double[nProductNum];
            for (int i = 0; i < nProductNum; i++)
            {
                nProductID[i] = Convert.ToInt32(lsYields[i][0]);
                dblOutput[i] = Convert.ToDouble(lsYields[i][2]);
            }
        }

        public int[] GetProductID()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (nProductID == null)
            {
                throw new Exception("Model ProductID is null!");
            }
            int[] nResult = new int[nProductID.Length];
            for (int i = 0; i < nProductID.Length; i++)
            {
                nResult[i] = nProductID[i];
            }
            return nResult;
        }

        public void ModelClose()
        {
            return;
        }

        /// <summary>
        /// 两分钟还没有算完则认为没有收敛，kill进程
        /// </summary>
        /// <param name="oProcess"></param>
        private void KillCoilsim(object oProcess)
        {
            Process p = (Process)oProcess;
            Thread.Sleep(2 * 60 * 1000);
            if (!p.HasExited)
            {
                p.Kill();
            }
        }
    }

    public class FuelGas : IModel
    {

        private double[] dblInput;//输入
        private double[] dblOutput;//输出
        private int[] nProductID;//产品ID
        //private string strExePath;//exe文件路径
        //private string strCasePath;//case文件路径
        //private string strCaseName;//case名
        private bool blErr = false;//设置参数错误标志

        public double[] GetOutput()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblOutput == null)
            {
                throw new Exception("Model Output is null!");
            }
            double[] dblResult = new double[dblOutput.Length];
            for (int i = 0; i < dblOutput.Length; i++)
            {
                dblResult[i] = dblOutput[i];
            }
            return dblResult;
        }

        public void SetInput(double[] dblVal)
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            if (dblVal.Length != 5)
            {
                throw new Exception("Input Value's dimension is not match with model!");
            }
            dblInput = new double[dblVal.Length];
            for (int i = 0; i < dblVal.Length; i++)
            {
                dblInput[i] = dblVal[i];
            }
        }

        public bool SetPara(string strVal)
        {
            if (!File.Exists(strVal))
            {
                blErr = true;
                throw new Exception("Input Directory is not valid!");
            }
            blErr = false;
            return true;
        }

        public bool HasParaSet()
        {
            return !blErr;
        }

        public void ModelCalculate()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            nProductID = new int[1] { -2 };
            dblOutput = new double[1];
            dblOutput[0] = dblInput.Sum(); 
        }

        public int[] GetProductID()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (nProductID == null)
            {
                throw new Exception("Model ProductID is null!");
            }
            int[] nResult = new int[nProductID.Length];
            for (int i = 0; i < nProductID.Length; i++)
            {
                nResult[i] = nProductID[i];
            }
            return nResult;
        }

        public void ModelClose()
        {
            return;
        }
    }

    public class SS : IModel
    {
        private double[] dblInput;//输入
        private double[] dblOutput;//输出
        private int[] nProductID;//产品ID
        //private string strExePath;//exe文件路径
        //private string strCasePath;//case文件路径
        //private string strCaseName;//case名
        private bool blErr = false;//设置参数错误标志

        public double[] GetOutput()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblOutput == null)
            {
                throw new Exception("Model Output is null!");
            }
            double[] dblResult = new double[dblOutput.Length];
            for (int i = 0; i < dblOutput.Length; i++)
            {
                dblResult[i] = dblOutput[i];
            }
            return dblResult;
        }


        public int[] GetProductID()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (nProductID == null)
            {
                throw new Exception("Model ProductID is null!");
            }
            int[] nResult = new int[nProductID.Length];
            for (int i = 0; i < nProductID.Length; i++)
            {
                nResult[i] = nProductID[i];
            }
            return nResult;
        }

        public bool HasParaSet()
        {
            return !blErr;
        }

        public void ModelCalculate()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            nProductID = new int[1] { 1000 };
            dblOutput = new double[1];
            dblOutput[0] = dblInput.Sum();
        }

        public void SetInput(double[] dblVal)
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            if (dblVal.Length != 3)
            {
                throw new Exception("Input Value's dimension is not match with model!");
            }
            dblInput = new double[dblVal.Length];
            for (int i = 0; i < dblVal.Length; i++)
            {
                dblInput[i] = dblVal[i];
            }
        }

        public bool SetPara(string strVal)
        {
            if (!File.Exists(strVal))
            {
                blErr = true;
                throw new Exception("Input Directory is not valid!");
            }
            blErr = false;
            return true;
        }
        public void ModelClose()
        {

        }

    }

    public class Aspen : IModel
    {
        private Happ.HappLS aspen;
        private AspenChromatography aspenc = new AspenChromatography();
        
        //private ServiceController scSQL; 
        private int InputNum = 0; // 输入个数
        private int OutputNum = 0; // 输出个数
        private double[] dblInput; // 输入
        private double[] dblOutput; // 输出
        private string Path_AspenFile = ""; // 输入变量路径
        private string[] Path_Input; // 输入变量路径
        private string[] Path_Output; // 输出变量路径
        private bool blErr = false; // 设置参数错误标志
        private string Store_strVal = "";
        private System.Threading.Thread timer = null;

        public bool SetPara(string strVal)
        {
            aspenc.ActiveDocument("asdsad");
            aspenc.Activate();
            
            
            if (!File.Exists(strVal))
            {
                blErr = true;
                throw new Exception("Input Directory is not valid!");
            }
            int tmpInputNum = 0;
            int tmpOutputNum = 0;
            int inputnum = 0;
            int outputnum = 0;
            List<string[]> Config = new List<string[]>();
            Config = CSVOperation.ReadCSV(strVal);
            Path_AspenFile = Config[0][0];
            for (int i = 0; i < Config.Count - 1; i++)
            {
                if (Config[1 + i][0] == "Input")
                {
                    tmpInputNum = tmpInputNum + 1;
                }
                else if (Config[1 + i][0] == "Output")
                {
                    tmpOutputNum = tmpOutputNum + 1;
                }
            }
            InputNum = tmpInputNum;
            OutputNum = tmpOutputNum;
            Path_Input = new string[InputNum];
            Path_Output = new string[OutputNum];
            dblInput = new double[InputNum];
            dblOutput = new double[OutputNum];
            for (int i = 0; i < Config.Count - 1; i++)
            {
                if (Config[1 + i][0] == "Input")
                {
                    Path_Input[inputnum] = Config[1 + i][1];
                    inputnum = inputnum + 1;
                }
                else if (Config[1 + i][0] == "Output")
                {
                    Path_Output[outputnum] = Config[1 + i][1];
                    outputnum = outputnum + 1;
                }
            }
            //scSQL = new ServiceController();
            //scSQL.MachineName = "localhost";
            //scSQL.ServiceName = "MSSQLSERVER";
            //if (scSQL.Status == ServiceControllerStatus.Stopped)
            //{
            //    scSQL.Start();
            //}
            aspen = new Happ.HappLS();
            aspen.InitFromFile2(Path_AspenFile);
            aspen.Visible = true;
            aspen.SuppressDialogs = 1;
            blErr = false;
            timer = null;
            Store_strVal = strVal;
            return true;
        }

        public bool HasParaSet()
        {
            return !blErr;
        }

        public void SetInput(double[] dblVal)
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            if (dblVal == null)
            {
                throw new Exception("Input Value is null!");
            }
            if (dblVal.Length != InputNum)
            {
                throw new Exception("Input Value's dimension is not match with model!");
            }
            for (int i = 0; i < dblVal.Length; i++)
            {
                dblInput[i] = dblVal[i];
            }
        }

        public void ModelCalculate()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model Parameters are not set!");
            }
            bool restart_flg = false;
            double[] tmpInput = new double[dblInput.Length];
            for (int i = 0; i < InputNum; i++)
            {
                tmpInput[i] = dblInput[i];
                aspen.Tree.FindNode(Path_Input[i]).Value = tmpInput[i];
            }
            System.Threading.Thread.Sleep(100);
            try
            {
                timer = new System.Threading.Thread(new System.Threading.ThreadStart(timecount));
                timer.Start();
                aspen.Run2();
                for (int i = 0; i < OutputNum; i++)
                {
                    dblOutput[i] = aspen.Tree.FindNode(Path_Output[i]).Value;
                }
            }
            catch (System.Exception ex)
            {
                Console.WriteLine("Crash, " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss"));
                if (timer != null && timer.ThreadState == System.Threading.ThreadState.WaitSleepJoin)
                {
                    timer.Abort();
                }
                timer = null;
                //scSQL.Stop();
                //scSQL.WaitForStatus(ServiceControllerStatus.Stopped);
                //scSQL.Close();
                aspen = null;
                SetPara(Store_strVal);
                Console.WriteLine("Restart, " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss"));
                restart_flg = true;
                System.Threading.Thread.Sleep(5000);
                throw new Exception(ex.Message);
            }
            if (restart_flg)
            {
                for (int i = 0; i < InputNum; i++)
                {
                    aspen.Tree.FindNode(Path_Input[i]).Value = tmpInput[i];
                }
                aspen.Run2();
                for (int i = 0; i < OutputNum; i++)
                {
                    dblOutput[i] = aspen.Tree.FindNode(Path_Output[i]).Value;
                }
                restart_flg = false;
            }
            if (timer != null && timer.ThreadState == System.Threading.ThreadState.WaitSleepJoin)
            {
                timer.Abort();
            }
            timer = null;
        }

        private void timecount()
        {
            System.Threading.Thread.Sleep(60 * 1000 *15);
            System.Diagnostics.Process[] proc = System.Diagnostics.Process.GetProcesses();
            foreach (System.Diagnostics.Process item in proc)
            {
                if (item.ProcessName == "AspenPlus")
                {
                    item.Kill();
                }
            }
        }

        public double[] GetOutput()
        {
            return dblOutput;
        }

        public int[] GetProductID()
        {
            return null;
        }

        public void ModelClose()
        {
            if (!HasParaSet())
            {
                throw new Exception("Model are not open!");
            }
            aspen.Quit();
            System.Threading.Thread.Sleep(5000);
            aspen.Quit();
            System.Diagnostics.Process[] proc = System.Diagnostics.Process.GetProcesses();
            foreach (System.Diagnostics.Process item in proc)
            {
                if (item.ProcessName == "AspenPlus")
                {
                    item.Kill();
                }
            }
            //scSQL.Stop();
            //scSQL.WaitForStatus(ServiceControllerStatus.Stopped);
            //scSQL.Close();  
        }
    }

    public class MyClass: IModel
    {

        public double[] GetOutput()
        {
            throw new NotImplementedException();
        }

        public int[] GetProductID()
        {
            throw new NotImplementedException();
        }

        public bool HasParaSet()
        {
            throw new NotImplementedException();
        }

        public void ModelCalculate()
        {
            throw new NotImplementedException();
        }

        public void SetInput(double[] dblVal)
        {
            throw new NotImplementedException();
        }

        public bool SetPara(string strVal)
        {
            throw new NotImplementedException();
        }

        public void ModelClose()
        {
            throw new NotImplementedException();
        }

    }
}