using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OptimizerSet;
using MathNet.Numerics.Statistics;
using System.Diagnostics;
using System.Reflection;
using Microsoft.VisualBasic.Devices;

namespace ConsoleApplication1
{
    public class test3
    {
        public void test()
        {
            string cfile = "C:\\Users\\Administrator\\Desktop\\Optimizer\\";
            General G = new General();
            NSGAIII NSGAIII = new NSGAIII();
            Stopwatch stopwatch = new Stopwatch();
            Random rd = new Random();
            int[] X_Dim = new int[] { 30 };
            int TotalRuns = 30;
            //double X_Range = 32.768;

            for (int ii = 0; ii < X_Dim.Length; ii++)
            {
                double[,] X_Opt; // 最优解
                double[] X_I = new double[X_Dim[ii]];
                double[] X_MaxStep = new double[X_Dim[ii]];
                double[] X_Lb = new double[X_Dim[ii]];
                double[] X_Ub = new double[X_Dim[ii]];

                for (int k = 0; k < X_Dim[ii]; k++)
                {
                    X_MaxStep[k] = 100.0;
                    X_Lb[k] = 0.0;
                    X_Ub[k] = 1.0;
                }

                double[,] NSGAIII_Fun;
                double[] NSGAIII_Time = new double[TotalRuns];
                string NSGAIII_Config = cfile + "NSGAIII.csv";
                Computer MyComputer = new Computer();
                
                for (int i = 0; i < TotalRuns; i++)
                {
                    for (int k = 0; k < X_Dim[ii]; k++)
                    {
                        X_I[k] = G.rnd_uni(rd) * (X_Ub[k] - X_Lb[k]) + X_Lb[k];
                    }

                    stopwatch.Start();
                    X_Opt = NSGAIII.StartOpt(NSGAIII_Config, X_I, X_MaxStep, X_Lb, X_Ub, out NSGAIII_Fun);
                    stopwatch.Stop();
                    NSGAIII_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("NSGAIII: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();
                    List<string[]> comeout = new List<string[]>();
                    for (int k = 0; k < NSGAIII_Fun.GetLength(0); k++)
                    {
                        string[] Fun = new string[2];
                        for (int j = 0; j < NSGAIII_Fun.GetLength(1); j++)
                        {
                            Fun[j] = NSGAIII_Fun[k, j].ToString();
                        }
                        comeout.Add(Fun);
                    }
                    string comeoutPath = cfile + "test3\\test3_" + i.ToString() + ".csv";
                    G.WriteCSV(comeoutPath, false, comeout);
                }
            }
        }

        private string[] data_process(string Algorithm, double[] Fun, double[] Time)
        {
            string[] Data = new string[8];
            Data[0] = Algorithm;
            Data[1] = Statistics.Mean(Fun).ToString();
            Data[2] = Statistics.StandardDeviation(Fun).ToString();
            Data[3] = Statistics.Minimum(Fun).ToString();
            Data[4] = Statistics.Maximum(Fun).ToString();
            Data[5] = Statistics.Mean(Time).ToString();
            Data[6] = Statistics.Minimum(Time).ToString();
            Data[7] = Statistics.Maximum(Time).ToString();
            return Data;
        }
    }
}
