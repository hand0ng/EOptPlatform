using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OptimizerSet;
using MathNet.Numerics.Statistics;
using System.Diagnostics;
using System.Reflection;

namespace ConsoleApplication1
{
    public class test1
    {
        public void test()
        {
            string cfile = "C:\\Users\\Administrator\\Desktop\\Optimizer\\";
            General G = new General();
            JADE JADE = new JADE();
            Nelder_Mead Nelder_Mead = new Nelder_Mead();
            CMA_ES CMA_ES = new CMA_ES();
            CoDE CoDE = new CoDE();
            PSO PSO = new PSO();
            PatternSearch PatternSearch = new PatternSearch();
            DE_rand_1 DE_rand_1 = new DE_rand_1();
            DE_best_1 DE_best_1 = new DE_best_1();
            DE_rand_2 DE_rand_2 = new DE_rand_2();
            DE_best_2 DE_best_2 = new DE_best_2();
            Stopwatch stopwatch = new Stopwatch();
            Random rd = new Random();
            List<string[]> comeout = new List<string[]>();
            int[] X_Dim = new int[] { 5, 10, 15, 20, 25, 30 };
            int TotalRuns = 30;
            double X_Range = 32.768;

            for (int ii = 0; ii < X_Dim.Length; ii++)
            {
                double[] X_Opt = new double[X_Dim[ii]]; // 最优解
                double[] X_I = new double[X_Dim[ii]];
                double[] X_MaxStep = new double[X_Dim[ii]];
                double[] X_Lb = new double[X_Dim[ii]];
                double[] X_Ub = new double[X_Dim[ii]];

                for (int k = 0; k < X_Dim[ii]; k++)
                {
                    X_MaxStep[k] = 100.0;
                    X_Lb[k] = -X_Range;
                    X_Ub[k] = X_Range;
                }

                double[] JADE_Fun = new double[TotalRuns];
                double[] JADE_Time = new double[TotalRuns];
                string JADE_Config = cfile + "JADE.csv";

                double[] Nelder_Mead_Fun = new double[TotalRuns];
                double[] Nelder_Mead_Time = new double[TotalRuns];
                string Nelder_Mead_Config = cfile + "Nelder_Mead.csv";

                double[] CMA_ES_Fun = new double[TotalRuns];
                double[] CMA_ES_Time = new double[TotalRuns];
                string CMA_ES_Config = cfile + "CMA_ES.csv";

                double[] CoDE_Fun = new double[TotalRuns];
                double[] CoDE_Time = new double[TotalRuns];
                string CoDE_Config = cfile + "CoDE.csv";

                double[] PSO_Fun = new double[TotalRuns];
                double[] PSO_Time = new double[TotalRuns];
                string PSO_Config = cfile + "PSO.csv";

                double[] PatternSearch_Fun = new double[TotalRuns];
                double[] PatternSearch_Time = new double[TotalRuns];
                string PatternSearch_Config = cfile + "PatternSearch.csv";

                double[] DE_rand_1_Fun = new double[TotalRuns];
                double[] DE_rand_1_Time = new double[TotalRuns];
                string DE_rand_1_Config = cfile + "DE_rand_1.csv";

                double[] DE_best_1_Fun = new double[TotalRuns];
                double[] DE_best_1_Time = new double[TotalRuns];
                string DE_best_1_Config = cfile + "DE_best_1.csv";

                double[] DE_rand_2_Fun = new double[TotalRuns];
                double[] DE_rand_2_Time = new double[TotalRuns];
                string DE_rand_2_Config = cfile + "DE_rand_2.csv";

                double[] DE_best_2_Fun = new double[TotalRuns];
                double[] DE_best_2_Time = new double[TotalRuns];
                string DE_best_2_Config = cfile + "DE_best_2.csv";

                for (int i = 0; i < TotalRuns; i++)
                {
                    for (int k = 0; k < X_Dim[ii]; k++)
                    {
                        X_I[k] = rnd_uni(rd) * X_Range;
                    }

                    stopwatch.Start();
                    X_Opt = JADE.StartOpt(JADE_Config, X_I, X_MaxStep, X_Lb, X_Ub, out JADE_Fun[i]);
                    stopwatch.Stop();
                    JADE_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("JADE: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = Nelder_Mead.StartOpt(Nelder_Mead_Config, X_I, X_MaxStep, X_Lb, X_Ub, out Nelder_Mead_Fun[i]);
                    stopwatch.Stop();
                    Nelder_Mead_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("Nelder_Mead: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = CMA_ES.StartOpt(CMA_ES_Config, X_I, X_MaxStep, X_Lb, X_Ub, out CMA_ES_Fun[i]);
                    stopwatch.Stop();
                    CMA_ES_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("CMA_ES: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = CoDE.StartOpt(CoDE_Config, X_I, X_MaxStep, X_Lb, X_Ub, out CoDE_Fun[i]);
                    stopwatch.Stop();
                    CoDE_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("CoDE: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = PSO.StartOpt(PSO_Config, X_I, X_MaxStep, X_Lb, X_Ub, out PSO_Fun[i]);
                    stopwatch.Stop();
                    PSO_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("PSO: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = PatternSearch.StartOpt(PatternSearch_Config, X_I, X_MaxStep, X_Lb, X_Ub, out PatternSearch_Fun[i]);
                    stopwatch.Stop();
                    PatternSearch_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("PatternSearch: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = DE_rand_1.StartOpt(DE_rand_1_Config, X_I, X_MaxStep, X_Lb, X_Ub, out DE_rand_1_Fun[i]);
                    stopwatch.Stop();
                    DE_rand_1_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("DE_rand_1: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = DE_best_1.StartOpt(DE_best_1_Config, X_I, X_MaxStep, X_Lb, X_Ub, out DE_best_1_Fun[i]);
                    stopwatch.Stop();
                    DE_best_1_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("DE_best_1: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = DE_rand_2.StartOpt(DE_rand_2_Config, X_I, X_MaxStep, X_Lb, X_Ub, out DE_rand_2_Fun[i]);
                    stopwatch.Stop();
                    DE_rand_2_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("DE_rand_2: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();

                    stopwatch.Start();
                    X_Opt = DE_best_2.StartOpt(DE_best_2_Config, X_I, X_MaxStep, X_Lb, X_Ub, out DE_best_2_Fun[i]);
                    stopwatch.Stop();
                    DE_best_2_Time[i] = stopwatch.Elapsed.TotalSeconds;
                    Console.WriteLine("DE_best_2: run: " + (i + 1) + ", time(s): " + stopwatch.Elapsed.TotalSeconds);
                    stopwatch.Reset();
                }
                string[] Dim = new string[2] { "X_DIM", X_Dim[ii].ToString() };
                string[] Table_Head1 = new string[8] { "Algorithm", "Mean", "StandardDeviation", "Min", "Max", "MeanTime(s)", "MinTime(s)", "MaxTime(s)" };
                comeout.Add(Dim);
                comeout.Add(Table_Head1);
                comeout.Add(data_process("JADE", JADE_Fun, JADE_Time));
                comeout.Add(data_process("Nelder_Mead", Nelder_Mead_Fun, Nelder_Mead_Time));
                comeout.Add(data_process("CMA_ES", CMA_ES_Fun, CMA_ES_Time));
                comeout.Add(data_process("CoDE", CoDE_Fun, CoDE_Time));
                comeout.Add(data_process("PSO", PSO_Fun, PSO_Time));
                comeout.Add(data_process("PatternSearch", PatternSearch_Fun, PatternSearch_Time));
                comeout.Add(data_process("DE_rand_1", DE_rand_1_Fun, DE_rand_1_Time));
                comeout.Add(data_process("DE_best_1", DE_best_1_Fun, DE_best_1_Time));
                comeout.Add(data_process("DE_rand_2", DE_rand_2_Fun, DE_rand_2_Time));
                comeout.Add(data_process("DE_best_2", DE_best_2_Fun, DE_best_2_Time));

                string comeoutPath = cfile + "test1_d" + X_Dim[ii].ToString() + ".csv";
                G.WriteCSV(comeoutPath, false, comeout);
            }
        }

        /// <summary>生成-1~1之间的随机实数</summary>
        private double rnd_uni(Random rd)
        {
            double randnumber = rd.Next(-10000, 10000);
            randnumber = randnumber / 10000;
            return randnumber;
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
