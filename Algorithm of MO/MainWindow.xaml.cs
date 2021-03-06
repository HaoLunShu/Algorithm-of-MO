﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;
using JMetalCSharp.Core;
using JMetalCSharp.QualityIndicator;
using JMetalCSharp.Utils;

namespace Algorithm_of_MO
{
    /// <summary>
    /// MainWindow.xaml 的互動邏輯
    /// </summary>
    public partial class MainWindow : Window
    {
        #region Constructor

        public MainWindow()
        {
            InitializeComponent();
        }

        #endregion

        #region Main Function

        public void MyAlgorithm()
        {
            Problem problem = null; // The problem to solve
            Algorithm algorithm = null; // The algorithm to use
            Operator crossover = null; // Crossover operator
            Operator mutation = null; // Mutation operator
            Operator selection= null; // Selection operator

            QualityIndicator indicators; // Object to get quality indicators

            FileRead fileRead = new FileRead();
            fileRead.fileread();
            problem = fileRead.GetProblem();
            algorithm = fileRead.GetAlgorithm();
            //crossover = fileRead.GetCrossover();
            mutation = fileRead.GetMutation();
            selection = fileRead.GetSelection();

            // Quality Indicators Operator
            //indicators = new QualityIndicator(problem, "DTLZ1.3D.pf");
            indicators = new QualityIndicator(problem, fileRead.Qi);
            //indicators = new QualityIndicator(problem, "LZ09_F1.pf");

            // Add the operators to the algorithm
            algorithm.AddOperator("crossover", crossover);
            algorithm.AddOperator("mutation", mutation);
            algorithm.AddOperator("selection", selection);

            // Add the indicator object to the algorithm
            algorithm.SetInputParameter("indicators", indicators);

            for (int i = 1; i <= int.Parse(fileRead.Rt); i++)
            {

                // Logger object and file to store log messages
                var logger = Logger.Log;

                var appenders = logger.Logger.Repository.GetAppenders();
                var fileAppender = appenders[0] as log4net.Appender.FileAppender;
                fileAppender.File = fileRead.DirPath + "/" + fileRead.Al + i +".log";
                fileAppender.ActivateOptions();

                string filevar = fileRead.DirPath + "/VAR" + i;
                string filefun = fileRead.DirPath + "/FUN" + i;

                FileStream file = new FileStream(fileRead.DirPath + "/" + fileRead.Al + "new" + i + ".txt", FileMode.OpenOrCreate, FileAccess.ReadWrite);
                StreamWriter newlogger = new StreamWriter(file);

                // Execute the Algorithm
                long initTime = Environment.TickCount;
                SolutionSet population = algorithm.Execute();
                long estimatedTime = Environment.TickCount - initTime;

                // Result messages 
                logger.Info("Total execution time: " + estimatedTime + "ms");
                logger.Info("Variables values have been writen to file " + filevar);
                newlogger.WriteLine("Total execution time: " + estimatedTime + "ms" + "\n");
                newlogger.WriteLine("Variables values have been writen to file " + filevar + "\n");
                population.PrintVariablesToFile(filevar);
                logger.Info("Objectives values have been writen to file " + filefun);
                population.PrintObjectivesToFile(filefun);
                Console.WriteLine("Time: " + estimatedTime);
                newlogger.WriteLine("Time: " + estimatedTime + "\n");
                Console.ReadLine();
                if (indicators != null)
                {
                    logger.Info("Quality indicators");
                    logger.Info("Hypervolume: " + indicators.GetHypervolume(population));
                    logger.Info("GD         : " + indicators.GetGD(population));
                    logger.Info("IGD        : " + indicators.GetIGD(population));
                    logger.Info("Spread     : " + indicators.GetSpread(population));
                    logger.Info("Epsilon    : " + indicators.GetEpsilon(population));

                    newlogger.WriteLine("Quality indicators");
                    newlogger.WriteLine("Hypervolume: " + indicators.GetHypervolume(population).ToString("F18") + "\n");
                    newlogger.WriteLine("GD         : " + indicators.GetGD(population).ToString("F18") + "\n");
                    newlogger.WriteLine("IGD        : " + indicators.GetIGD(population).ToString("F18") + "\n");
                    newlogger.WriteLine("Spread     : " + indicators.GetSpread(population).ToString("F18") + "\n");
                    newlogger.WriteLine("Epsilon    : " + indicators.GetEpsilon(population).ToString("F18") + "\n");

                    int evaluations = (int)algorithm.GetOutputParameter("evaluations");
                    logger.Info("Speed      : " + evaluations + "     evaluations");
                    newlogger.WriteLine("Speed      : " + evaluations + "     evaluations" + "\n");
                }
                newlogger.Close();
                file.Close();
                algorithm.AddOperator("mutation", mutation);
            }
        }

        #endregion

        #region Methods

        public void calculateStatistics(string filepath, string algorithm, int numbers)
        {
            string p = filepath + "/" + algorithm + "new";
            string[] s = { "Time", "Hypervolume", "GD", "IGD", "Spread" , "Epsilon", "Speed" };
            double[] allParameters = new double[s.Length];
            int[] T = new int[numbers];
            double[,] QI = new double[5, numbers];
            List<double> Q = new List<double>();
            int[] spe = new int[numbers];
            int flag = 0;
            for(int i = 1; i <= numbers; i++)
            {
                string text = System.IO.File.ReadAllText(p + i + ".txt");
                string text1;
                System.IO.StreamReader sr = new System.IO.StreamReader(p + i + ".txt");
                int t = 0;
                while ((text1 = sr.ReadLine()) != null)
                {
                    string[] words = text1.Split(' ');
                    if (words[0] == "Time:")
                    {
                        t = int.Parse(words[1]);
                        break;
                    }
                }
                T[flag] = t;
                allParameters[0] = allParameters[0] + t;

                bool[] b = new bool[s.Length];
                string[] d = new string[s.Length];
                double[] para = new double[s.Length];
                for (int j = 1; j < s.Length - 1; j++)
                {
                    b[j] = text.Contains(s[j]);
                    if (b[j])
                    {
                        int index = text.IndexOf(s[j]);
                        if (index >= 0)
                        {
                            d[j] = text.Substring(index + 13, 20);
                            para[j] = double.Parse(d[j]);
                            QI[j-1, flag] = para[j];
                            allParameters[j] = allParameters[j] + para[j];
                        }
                    }
                }

                bool bs = text.Contains(s[6]);
                if (bs)
                {
                    int indexspeed = text.IndexOf(s[6]);
                    if (indexspeed >= 0)
                    {
                        string speed = text.Substring(indexspeed + 13, 5);
                        int sp = Int32.Parse(speed);
                        spe[flag] = sp;
                        allParameters[6] = allParameters[6] + sp;
                    }
                }

                flag++;

            }

            double[] med = new double[s.Length];
            med = median(T, QI, spe);
            double[] iqr = new double[s.Length];
            iqr = IQR(T, QI, spe);
            double[] sd = new double[s.Length];
            sd = SD(T, QI, spe);

            for (int k = 0; k < allParameters.Length; k++)
            {
                Console.WriteLine("Average " + s[k] + ": " + allParameters[k] / numbers);
                File.AppendAllText(filepath + "/finalQI" + ".txt", "Average " + s[k] + ": " + allParameters[k] / numbers + "\n");
                Console.WriteLine("Standard Deviation " + s[k] + ": " + sd[k]);
                //System.Data.Objects.EntityFunctions.StandardDeviationP(QI)
                File.AppendAllText(filepath + "/finalQI" + ".txt", "Standard Deviation " + s[k] + ": " + sd[k] + "\n");
                Console.WriteLine("Median " + s[k] + ": " + med[k]);
                File.AppendAllText(filepath + "/finalQI" + ".txt", "Median " + s[k] + ": " + med[k] + "\n");
                Console.WriteLine("IQR " + s[k] + ": " + iqr[k]);
                File.AppendAllText(filepath + "/finalQI" + ".txt", "IQR " + s[k] + ": " + iqr[k] + "\n");
            }

            File.AppendAllText(filepath + "/finalQI" + ".txt", "\n");

        }

        private double[] SD(int[] Time, double[,] Q, int[] SP)
        {
            double[] hy = new double[Time.Length];
            double[] gd = new double[Time.Length];
            double[] igd = new double[Time.Length];
            double[] spr = new double[Time.Length];
            double[] eps = new double[Time.Length];
            double[] avg = new double[7];
            double[] sum = new double[7];
            double[] result = new double[7];

            for (int l = 0; l < Time.Length; l++)
            {
                hy[l] = Q[0, l];
                gd[l] = Q[1, l];
                igd[l] = Q[2, l];
                spr[l] = Q[3, l];
                eps[l] = Q[4, l];

            }

            avg[0] = Time.Average();
            avg[1] = hy.Average();
            avg[2] = gd.Average();
            avg[3] = igd.Average();
            avg[4] = spr.Average();
            avg[5] = eps.Average();
            avg[6] = SP.Average();

            double[] sumOfSquaresOfDifferences = new double[7];
            
            sumOfSquaresOfDifferences[0] = Time.Sum(val => (val - avg[0]) * (val - avg[0]));
            sumOfSquaresOfDifferences[1] = hy.Sum(val => (val - avg[1]) * (val - avg[1]));
            sumOfSquaresOfDifferences[2] = gd.Sum(val => (val - avg[2]) * (val - avg[2]));
            sumOfSquaresOfDifferences[3] = igd.Sum(val => (val - avg[3]) * (val - avg[3]));
            sumOfSquaresOfDifferences[4] = spr.Sum(val => (val - avg[4]) * (val - avg[4]));
            sumOfSquaresOfDifferences[5] = eps.Sum(val => (val - avg[5]) * (val - avg[5]));
            sumOfSquaresOfDifferences[6] = SP.Sum(val => (val - avg[6]) * (val - avg[6]));

            for(int y = 0; y < 7; y++)
            {
                result[y] = Math.Sqrt(sumOfSquaresOfDifferences[y] / Time.Length);
            }
            
            return result;

        }

        private double[] median(int[] Time, double[,] Q, int[] SP)
        {
            double[] hy = new double[Time.Length];
            double[] gd = new double[Time.Length];
            double[] igd = new double[Time.Length];
            double[] spr = new double[Time.Length];
            double[] eps = new double[Time.Length];
            double[] result = new double[7];

            for(int l = 0; l < Time.Length; l++)
            {
                hy[l] = Q[0, l];
                gd[l] = Q[1, l];
                igd[l] = Q[2, l];
                spr[l] = Q[3, l];
                eps[l] = Q[4, l];

            }

            Array.Sort(Time);
            Array.Sort(hy);
            Array.Sort(gd);
            Array.Sort(igd);
            Array.Sort(spr);
            Array.Sort(eps);
            Array.Sort(SP);
            
            if(Time.Length % 2 == 0)
            {
                int y = Time.Length / 2 - 1;
                result[0] = (Time[y] + Time[y + 1]) / 2;
                result[1] = (hy[y] + hy[y + 1]) / 2;
                result[2] = (gd[y] + gd[y + 1]) / 2;
                result[3] = (igd[y] + igd[y + 1]) / 2;
                result[4] = (spr[y] + spr[y + 1]) / 2;
                result[5] = (eps[y] + eps[y + 1]) / 2;
                result[6] = (SP[y] + SP[y + 1]) / 2;
            }
            else
            {
                int x = Time.Length / 2;
                result[0] = Time[x];
                result[1] = hy[x];
                result[2] = gd[x];
                result[3] = igd[x];
                result[4] = spr[x];
                result[5] = eps[x];
                result[6] = SP[x];
            }

            return result;
            
        }

        private double[] IQR(int[] Time, double[,] Q, int[] SP)
        {
            double[] hy = new double[Time.Length];
            double[] gd = new double[Time.Length];
            double[] igd = new double[Time.Length];
            double[] spr = new double[Time.Length];
            double[] eps = new double[Time.Length];
            double[] Q1 = new double[7];
            double[] Q3 = new double[7];
            double[] result = new double[7];

            for (int l = 0; l < Time.Length; l++)
            {
                hy[l] = Q[0, l];
                gd[l] = Q[1, l];
                igd[l] = Q[2, l];
                spr[l] = Q[3, l];
                eps[l] = Q[4, l];

            }

            Array.Sort(Time);
            Array.Sort(hy);
            Array.Sort(gd);
            Array.Sort(igd);
            Array.Sort(spr);
            Array.Sort(eps);
            Array.Sort(SP);

            if (Time.Length % 4 == 0)
            {
                int y = Time.Length / 4 - 1;
                Q1[0] = (Time[y] + Time[y + 1]) / 2;
                Q1[1] = (hy[y] + hy[y + 1]) / 2;
                Q1[2] = (gd[y] + gd[y + 1]) / 2;
                Q1[3] = (igd[y] + igd[y + 1]) / 2;
                Q1[4] = (spr[y] + spr[y + 1]) / 2;
                Q1[5] = (eps[y] + eps[y + 1]) / 2;
                Q1[6] = (SP[y] + SP[y + 1]) / 2;
            }
            else
            {
                int x = Time.Length / 4;
                Q1[0] = Time[x];
                Q1[1] = hy[x];
                Q1[2] = gd[x];
                Q1[3] = igd[x];
                Q1[4] = spr[x];
                Q1[5] = eps[x];
                Q1[6] = SP[x];
            }

            if ((Time.Length * 3) % 4 == 0)
            {
                int y = Time.Length / 4 * 3 - 1;
                Q3[0] = (Time[y] + Time[y + 1]) / 2;
                Q3[1] = (hy[y] + hy[y + 1]) / 2;
                Q3[2] = (gd[y] + gd[y + 1]) / 2;
                Q3[3] = (igd[y] + igd[y + 1]) / 2;
                Q3[4] = (spr[y] + spr[y + 1]) / 2;
                Q3[5] = (eps[y] + eps[y + 1]) / 2;
                Q3[6] = (SP[y] + SP[y + 1]) / 2;
            }
            else
            {
                int x = Time.Length / 4 * 3;
                Q3[0] = Time[x];
                Q3[1] = hy[x];
                Q3[2] = gd[x];
                Q3[3] = igd[x];
                Q3[4] = spr[x];
                Q3[5] = eps[x];
                Q3[6] = SP[x];
            }

            for(int z = 0; z < result.Length; z++)
            {
                result[z] = Q3[z] - Q1[z];
            }

            return result;

        }

        #endregion

        #region Windows function

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            MyAlgorithm();
        }

        private void Button_Click_1(object sender, RoutedEventArgs e)
        {
            FileRead file = new FileRead();
            file.fileread();
            if (file.Co2 != "")
                calculateStatistics(file.DirPath, file.Al, int.Parse(file.Rt));
            else
                calculateStatistics(file.DirPath, file.Al, int.Parse(file.Rt));
        }

        private void Button_Click_2(object sender, RoutedEventArgs e)
        {
            FileRead read = new FileRead();
            read.fileread();
            Problem problem = read.GetProblem();
            Algorithm algorithm = read.GetAlgorithm();
            QualityIndicator indicators = new QualityIndicator(problem, read.Qi);

            NewAlgorithmTest test = new NewAlgorithmTest();

            for (int i = 1; i <= int.Parse(read.Rt); i++)
            {

                // Logger object and file to store log messages
                var logger = Logger.Log;

                var appenders = logger.Logger.Repository.GetAppenders();
                var fileAppender = appenders[0] as log4net.Appender.FileAppender;
                fileAppender.File = read.DirPath + "/" + read.Al + i + ".log";
                fileAppender.ActivateOptions();

                string filevar = read.DirPath + "/VAR" + i;
                string filefun = read.DirPath + "/FUN" + i;

                FileStream file = new FileStream(read.DirPath + "/MOEADnew" + i + ".txt", FileMode.OpenOrCreate, FileAccess.ReadWrite);
                StreamWriter newlogger = new StreamWriter(file);

                // Execute the Algorithm
                long initTime = Environment.TickCount;
                SolutionSet population = test.Mix();
                long estimatedTime = Environment.TickCount - initTime;

                // Result messages 
                logger.Info("Total execution time: " + estimatedTime + "ms");
                logger.Info("Variables values have been writen to file " + filevar);
                newlogger.WriteLine("Total execution time: " + estimatedTime + "ms" + "\n");
                newlogger.WriteLine("Variables values have been writen to file " + filevar + "\n");
                population.PrintVariablesToFile(filevar);
                logger.Info("Objectives values have been writen to file " + filefun);
                population.PrintObjectivesToFile(filefun);
                Console.WriteLine("Time: " + estimatedTime);
                newlogger.WriteLine("Time: " + estimatedTime + "\n");
                Console.ReadLine();
                if (indicators != null)
                {
                    logger.Info("Quality indicators");
                    logger.Info("Hypervolume: " + indicators.GetHypervolume(population));
                    logger.Info("GD         : " + indicators.GetGD(population));
                    logger.Info("IGD        : " + indicators.GetIGD(population));
                    logger.Info("Spread     : " + indicators.GetSpread(population));
                    logger.Info("Epsilon    : " + indicators.GetEpsilon(population));

                    newlogger.WriteLine("Quality indicators");
                    newlogger.WriteLine("Hypervolume: " + indicators.GetHypervolume(population).ToString("F18") + "\n");
                    newlogger.WriteLine("GD         : " + indicators.GetGD(population).ToString("F18") + "\n");
                    newlogger.WriteLine("IGD        : " + indicators.GetIGD(population).ToString("F18") + "\n");
                    newlogger.WriteLine("Spread     : " + indicators.GetSpread(population).ToString("F18") + "\n");
                    newlogger.WriteLine("Epsilon    : " + indicators.GetEpsilon(population).ToString("F18") + "\n");

                    int evaluations = test.GetEvaluations();
                    logger.Info("Speed      : " + evaluations + "     evaluations");
                    newlogger.WriteLine("Speed      : " + evaluations + "     evaluations" + "\n");
                }
                newlogger.Close();
                file.Close();
            }
        }

        #endregion
    }
}
