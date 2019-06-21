using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;

namespace ModelSet
{
    class CSVOperation
    {
        //write a file, existed file will be overwritten if append = false
        public static void WriteCSV(string filePathName, bool append, List<string[]> ls)
        {
            StreamWriter fileWriter = new StreamWriter(filePathName, append, Encoding.Default);
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

        //read a file.
        public static List<string[]> ReadCSV(string filePathName)
        {
            List<string[]> ls = new List<string[]>();
            FileStream fs = new FileStream(filePathName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            StreamReader fileReader = new StreamReader(fs,Encoding.Default);
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
                        string[] str1 = strLine.Split(new char[2] { '\t', ' ' });
                        List<string> lsStr = new List<string>();
                        for (int i = 0; i < str1.Length; i++)
                        {
                            if (str1[i] != "")
                            {
                                lsStr.Add(str1[i]);
                            }
                        }
                        if (lsStr.Count > 0)
                        {
                            ls.Add(lsStr.ToArray());
                        }
                    }
                }
            }
            fileReader.Close();
            return ls;
        }
    }
}
