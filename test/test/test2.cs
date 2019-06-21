using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OptimizerSet;
using MathNet.Numerics.Statistics;
using System.Diagnostics;
using System.Reflection;
using System.IO;

namespace ConsoleApplication1
{
    public class test2
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
            List<string[]> fun = new List<string[]>();
            List<string[]> time = new List<string[]>();
            int[] X_Dim = new int[20];
            int TotalRuns = 30;
            for (int i = 0; i < X_Dim.Length; i++)
            {
                X_Dim[i] = 2 * (i + 3);
            }

            string[] Table_Head2 = new string[11] { "X_Dim", "JADE", "Nelder_Mead", "CMA_ES", "CoDE", "PSO", "PatternSearch", "DE_rand_1", "DE_best_1", "DE_rand_2", "DE_best_2" };
            fun.Add(Table_Head2);
            time.Add(Table_Head2);

            for (int ii = 10; ii < X_Dim.Length; ii++)
            {
                Console.WriteLine("X_Dim: " + X_Dim[ii]);
                double[] X_Opt = new double[X_Dim[ii]]; // 最优解
                double[] X_I = new double[X_Dim[ii]];
                double[] X_MaxStep = new double[X_Dim[ii]];
                double[] X_Lb = new double[X_Dim[ii]];
                double[] X_Ub = new double[X_Dim[ii]];

                for (int k = 0; k < X_Dim[ii] - 1; k++)
                {
                    X_MaxStep[k] = 100.0;
                    X_Lb[k] = 0;
                }
                X_MaxStep[X_Dim[ii] - 1] = 100.0;
                X_Lb[X_Dim[ii] - 1] = Math.PI;
                for (int k = 0; k < X_Dim[ii] / 2 - 1; k++)
                {
                    X_Ub[k] = 1;
                }
                X_Ub[X_Dim[ii] / 2 - 1] = 0;
                for (int k = X_Dim[ii] / 2; k < X_Dim[ii]; k++)
                {
                    X_Ub[k] = Math.PI;
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
                        X_I[k] = G.rnd_uni(rd) * (X_Ub[k] - X_Lb[k]) + X_Lb[k];
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

                if (true)
                {
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

                    string comeoutPath = cfile + "comeout_" + X_Dim[ii].ToString() + ".csv";
                    G.WriteCSV(comeoutPath, false, comeout);

                    fun.Add(data_process1(X_Dim[ii], JADE_Fun, Nelder_Mead_Fun, CMA_ES_Fun, CoDE_Fun, PSO_Fun, PatternSearch_Fun, DE_rand_1_Fun, DE_best_1_Fun, DE_rand_2_Fun, DE_best_2_Fun));
                    time.Add(data_process1(X_Dim[ii], JADE_Time, Nelder_Mead_Time, CMA_ES_Time, CoDE_Time, PSO_Time, PatternSearch_Time, DE_rand_1_Time, DE_best_1_Time, DE_rand_2_Time, DE_best_2_Time));
                }
            }
            string comeoutPath1 = cfile + "fun.csv";
            G.WriteCSV(comeoutPath1, false, fun);
            string comeoutPath2 = cfile + "time.csv";
            G.WriteCSV(comeoutPath2, false, time);
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

        private string[] data_process1(int dim, double[] Fun1, double[] Fun2, double[] Fun3, double[] Fun4, double[] Fun5, double[] Fun6, double[] Fun7, double[] Fun8, double[] Fun9, double[] Fun10)
        {
            string[] Data = new string[11];
            Data[0] = dim.ToString();
            Data[1] = Statistics.Mean(Fun1).ToString();
            Data[2] = Statistics.Mean(Fun2).ToString();
            Data[3] = Statistics.Mean(Fun3).ToString();
            Data[4] = Statistics.Mean(Fun4).ToString();
            Data[5] = Statistics.Mean(Fun5).ToString();
            Data[6] = Statistics.Mean(Fun6).ToString();
            Data[7] = Statistics.Mean(Fun7).ToString();
            Data[8] = Statistics.Mean(Fun8).ToString();
            Data[9] = Statistics.Mean(Fun9).ToString();
            Data[10] = Statistics.Mean(Fun10).ToString();
            return Data;
        }
    }
}
