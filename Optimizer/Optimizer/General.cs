using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace OptimizerSet
{
    public class General
    {
        /// <summary>边界约束</summary>
        public double BoundConstraint_CoDE(double a, double amin, double amax)
        {
            a = (a < amin) ? Math.Min(amax, 2 * amin - a) : a;
            a = (a > amax) ? Math.Max(amin, 2 * amax - a) : a;
            return a;
        }

        /// <summary>边界约束</summary>
        public double BoundConstraint_CMA_ES(double a, double amin, double amax)
        {
            if (a < amin)
            {
                a = 2 * amin - a;
                if (a > amax)
                {
                    a = amax;
                }
            }
            if (a > amax)
            {
                a = 2 * amax - a;
                if (a < amin)
                {
                    a = amin;
                }
            }
            return a;
        }
        
        /// <summary>边界约束</summary>
        public double BoundConstraint_JADE(double old_a, double new_a, double amin, double amax)
        {
            new_a = (new_a < amin) ? (amin + old_a) / 2 : new_a;
            new_a = (new_a > amax) ? (amax + old_a) / 2 : new_a;
            return new_a;
        }

        /// <summary>边界约束</summary>
        public double BoundConstraint_PatternSearch(double old_a, double new_a, double amin, double amax)
        {
            //new_a = (new_a < amin) ? old_a - (old_a - amin) * 0.8 : new_a;
            //new_a = (new_a > amax) ? old_a + (amax - old_a) * 0.8 : new_a;
            new_a = (new_a < amin) ? Math.Min(amax, 2 * amin - new_a) : new_a;
            new_a = (new_a > amax) ? Math.Max(amin, 2 * amax - new_a) : new_a;
            return new_a;
        }

        /// <summary>边界约束</summary>
        public double BoundConstraint_DE(double old_a, double new_a, double amin, double amax)
        {
            new_a = (new_a < amin) ? Math.Min(amax, 2 * amin - new_a) : new_a;
            new_a = (new_a > amax) ? Math.Max(amin, 2 * amax - new_a) : new_a;
            return new_a;
        }

        /// <summary>边界约束</summary>
        public double BoundConstraint(double a, double amin, double amax)
        {
            a = (a < amin) ? amin : a;
            a = (a > amax) ? amax : a;
            return a;
        }

        /// <summary>数组最大值编号</summary>
        public int maxIndex(double[] a)
        {
            double a_max = a.Max();
            int Index = 0;
            for (int i = 0; i < a.Length; i++)
            {
                if (a_max == a[i])
                {
                    Index = i;
                }
            }
            return Index;
        }

        /// <summary>数组最大值编号数组</summary>
        public ArrayList maxIndexArray(int[] a)
        {
            int a_max = a.Max();
            ArrayList Index = new ArrayList();
            for (int i = 0; i < a.Length; i++)
            {
                if (a_max == a[i])
                {
                    Index.Add(i);
                }
            }
            return Index;
        }

        /// <summary>数组最小值编号</summary>
        public int minIndex(double[] a)
        {
            double a_min = a.Min();
            int Index = 0;
            for (int i = 0; i < a.Length; i++)
            {
                if (a_min == a[i])
                {
                    Index = i;
                }
            }
            return Index;
        }

        /// <summary>数组最小值编号数组</summary>
        public ArrayList minIndexArray(int[] a)
        {
            int a_min = a.Min();
            ArrayList Index = new ArrayList();
            for (int i = 0; i < a.Length; i++)
            {
                if (a_min == a[i])
                {
                    Index.Add(i);
                }
            }
            return Index;
        }

        /// <summary>
        /// 生成0~1之间的随机实数
        /// </summary>
        public double rnd_uni(Random rd)
        {
            //double randnumber = rd.Next(0, 10000);
            //randnumber = randnumber / 10000;
            //return randnumber;
            //double randnumber = rndInt_uni(0, 10000, rd);
            //randnumber = randnumber / 10000;
            //return randnumber;
            double randnumber = rd.NextDouble();
            return randnumber;
        }

        /// <summary>生成x1~x2之间的随机整数</summary>
        public int rndInt_uni(int x1, int x2, Random rd)
        {
            int randnumber;
            randnumber = rd.Next() % (x2 - x1 + 1) + x1;
            return randnumber;
        }

        /// <summary>生成正态分布随机实数</summary>
        public double rnd_normal(double averageValue, double standardDeviation, Random rd)
        {
            return averageValue + standardDeviation * (float)
                (Math.Sqrt(-2.0 * Math.Log(1.0 - rd.NextDouble())) * Math.Sin((Math.PI + Math.PI) * rd.NextDouble()));
        }

        /// <summary>随机打乱数组</summary>
        public double[] arrrandom(double[] arr, Random rd)
        {
            int k = 0;
            double dbltmp = 0.0;
            for (int i = 0; i < arr.Length; i++)
            {
                k = rndInt_uni(0, arr.Length - 1, rd);
                if (k != i)
                {
                    dbltmp = arr[i];
                    arr[i] = arr[k];
                    arr[k] = dbltmp;
                }
            }
            return arr;
        }

        /// <summary>
        /// 冒泡排序
        /// </summary>
        private void BubbleSort(double[] Fun, double[] Delta, int[] Index)
        {
            bool flag;
            double Fun_tmp = 0.0;
            double Delta_tmp = 0.0;
            int Index_tmp = 0;
            int N = Fun.GetLength(0);
            for (int i = N - 1; i > 0; i--)
            {
                flag = true;
                for (int j = 0; j < i; j++)
                {
                    if ((Delta[j] == Delta[j + 1] && Fun[j] > Fun[j + 1]) || Delta[j] > Delta[j + 1])
                    {
                        Fun_tmp = Fun[j];
                        Fun[j] = Fun[j + 1];
                        Fun[j + 1] = Fun_tmp;
                        Delta_tmp = Delta[j];
                        Delta[j] = Delta[j + 1];
                        Delta[j + 1] = Delta_tmp;
                        Index_tmp = Index[j];
                        Index[j] = Index[j + 1];
                        Index[j + 1] = Index_tmp;
                        flag = false;
                    }
                }
                if (flag) break;
            }
        }

        /// <summary>读取CSV文件</summary>
        public List<string[]> ReadCSV(string filePathName)
        {
            List<string[]> ls = new List<string[]>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string[] strSlot = filePathName.Split('.');
            string strFormat = strSlot[strSlot.Length - 1];
            string strLine = "";
            while (strLine != null)
            {
                strLine = fileReader.ReadLine();
                if (strLine != null && strLine.Length > 0)
                {
                    if ("csv" == strFormat)
                    {
                        ls.Add(strLine.Split(','));
                    }
                    if ("txt" == strFormat)
                    {
                        ls.Add(strLine.Split(' '));
                    }
                }
            }
            fileReader.Close();
            GC.Collect();
            return ls;
        }

        /// <summary>写CSV文件</summary>
        public void WriteCSV(string filePathName, bool append, List<string[]> ls)
        {
            StreamWriter fileWriter = new StreamWriter(filePathName, append, Encoding.Default); // true为追加，false为重写
            string[] strSlot = filePathName.Split('.');
            string strFormat = strSlot[strSlot.Length - 1];
            foreach (string[] strArr in ls)
            {
                if ("csv" == strFormat)
                {
                    fileWriter.WriteLine(String.Join(",", strArr));
                }
                if ("txt" == strFormat)
                {
                    fileWriter.WriteLine(String.Join(" ", strArr));
                }
            }
            fileWriter.Flush();
            fileWriter.Close();
        }
    }
}