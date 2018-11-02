using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using JMetalCSharp;
using JMetalCSharp.Core;
using JMetalCSharp.Operators.Crossover;
using JMetalCSharp.Operators.Mutation;
using JMetalCSharp.Operators.Selection;
using JMetalCSharp.Problems;
using JMetalCSharp.Problems.DTLZ;
using JMetalCSharp.Problems.ZDT;
using JMetalCSharp.Problems.LZ09;
using JMetalCSharp.QualityIndicator;
using JMetalCSharp.Utils;

namespace Algorithm_of_MO
{
    /// <summary>
    /// MainWindow.xaml 的互動邏輯
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        public void MyAlgorithm()
        {
            Problem problem; // The problem to solve
            Algorithm algorithm; // The algorithm to use
            Operator crossover; // Crossover operator
            Operator mutation; // Mutation operator
            Operator selection; // Selection operator

            Dictionary<string, object> parameters; // Operator parameters

            QualityIndicator indicators; // Object to get quality indicators

            indicators = null;
            // Default problem
            //problem = new Kursawe("Real", 3);
            //problem = new Kursawe("BinaryReal", 3);
            //problem = new Water("Real");
            //problem = new ZDT3("ArrayReal", 30);
            problem = new LZ09_F1("Real");
            //problem = new ConstrEx("Real");
            //problem = new DTLZ1("Real", 10, 3);
            //problem = new OKA2("Real") ;

            //algorithm = new JMetalCSharp.Metaheuristics.NSGAII.NSGAII(problem);
            //algorithm = new ssNSGAII(problem);
            algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);

            // Algorithm parameters
            algorithm.SetInputParameter("populationSize", 101);
            algorithm.SetInputParameter("maxEvaluations", 30000);
            algorithm.SetInputParameter("iterationsNumber", 250);

            algorithm.SetInputParameter("dataDirectory", "Data/Parameters/Weight");

            //algorithm.SetInputParameter("finalSize", 300); // used by MOEAD_DRA

            algorithm.SetInputParameter("T", 20);
            algorithm.SetInputParameter("delta", 0.9);
            algorithm.SetInputParameter("nr", 2);

            // Mutation and Crossover for Real codification 
            parameters = new Dictionary<string, object>();
            parameters.Add("CR", 1.0);
            parameters.Add("F", 0.5);
            parameters.Add("K", 1.0);
            crossover = CrossoverFactory.GetCrossoverOperator("DifferentialEvolutionCrossover", parameters);
            /*parameters.Add("probability", 1.0);
            parameters.Add("distributionIndex", 20.0);
            crossover = CrossoverFactory.GetCrossoverOperator("SBXCrossover", parameters);*/

            parameters = new Dictionary<string, object>();
            parameters.Add("probability", 1.0 / problem.NumberOfVariables);
            parameters.Add("distributionIndex", 20.0);
            mutation = MutationFactory.GetMutationOperator("PolynomialMutation", parameters);

            // Selection Operator 
            parameters = null;
            selection = SelectionFactory.GetSelectionOperator("BinaryTournament2", parameters);

            // Quality Indicators Operator
            //indicators = new QualityIndicator(problem, "DTLZ1.3D.pf");
            //indicators = new QualityIndicator(problem, "ZDT3.pf");
            indicators = new QualityIndicator(problem, "LZ09_F1.pf");

            // Add the operators to the algorithm
            algorithm.AddOperator("crossover", crossover);
            algorithm.AddOperator("mutation", mutation);
            algorithm.AddOperator("selection", selection);

            // Add the indicator object to the algorithm
            algorithm.SetInputParameter("indicators", indicators);

            for (int i = 1; i <= 100; i++)
            {

                // Logger object and file to store log messages
                var logger = Logger.Log;

                var appenders = logger.Logger.Repository.GetAppenders();
                var fileAppender = appenders[0] as log4net.Appender.FileAppender;
                fileAppender.File = "Result/MOEAD_DE/LZ09_F1/MOEAD_DE" + i +".log";
                fileAppender.ActivateOptions();

                string filevar = "Result/MOEAD_DE/LZ09_F1/VAR" + i;
                string filefun = "Result/MOEAD_DE/LZ09_F1/FUN" + i;

                FileStream file = new FileStream("Result/MOEAD_DE/LZ09_F1/MOEAD_DEnew" + i + ".txt", FileMode.OpenOrCreate, FileAccess.ReadWrite);
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
                    newlogger.WriteLine("Hypervolume: " + indicators.GetHypervolume(population) + "\n");
                    newlogger.WriteLine("GD         : " + indicators.GetGD(population) + "\n");
                    newlogger.WriteLine("IGD        : " + indicators.GetIGD(population) + "\n");
                    newlogger.WriteLine("Spread     : " + indicators.GetSpread(population) + "\n");
                    newlogger.WriteLine("Epsilon    : " + indicators.GetEpsilon(population) + "\n");

                    int evaluations = (int)algorithm.GetOutputParameter("evaluations");
                    logger.Info("Speed      : " + evaluations + "     evaluations");
                    newlogger.WriteLine("Speed      : " + evaluations + "     evaluations" + "\n");
                }
                newlogger.Close();
                file.Close();
            }
        }

        public void calculateStatistics(string filepath, string algorithm, int numbers)
        {
            string p = filepath + "/" + algorithm + "new";
            string[] s = { "Time", "Hypervolume", "GD", "IGD", "Spread" , "Epsilon", "Speed" };
            double[] allParameters = new double[s.Length];
            int[] T = new int[numbers];
            double[,] QI = new double[5, numbers];
            int[] spe = new int[numbers];
            int flag = 0;
            for(int i = 1; i <= numbers; i++)
            {
                System.IO.StreamReader sr = new System.IO.StreamReader(p + i + ".txt");
                string text;
                while ((text = sr.ReadLine()) != null)
                {
                    bool bt = text.Contains(s[0]);
                    if (bt)
                    {
                        int indextime = text.LastIndexOf(s[0]);
                        if (indextime >= 0)
                        {
                            string time = text.Substring(indextime + 6, 4);
                            int t = Int32.Parse(time);
                            T[flag] = t;
                            allParameters[0] = allParameters[0] + t;
                        }
                    }

                    bool[] b = new bool[s.Length];
                    string[] d = new string[s.Length];
                    double[] para = new double[s.Length];
                    for (int j = 1; j < s.Length - 1; j++)
                    {
                        b[j] = text.Contains(s[j]);
                        if (b[j])
                        {
                            int index = text.LastIndexOf(s[j]);
                            if (index >= 0)
                            {
                                d[j] = text.Substring(index + 13);
                                para[j] = double.Parse(d[j]);
                                QI[j-1, flag] = para[j];
                                allParameters[j] = allParameters[j] + para[j];
                            }
                        }
                    }

                    bool bs = text.Contains(s[6]);
                    if (bs)
                    {
                        int indexspeed = text.LastIndexOf(s[6]);
                        if (indexspeed >= 0)
                        {
                            string speed = text.Substring(indexspeed + 13, 5);
                            int sp = Int32.Parse(speed);
                            spe[flag] = sp;
                            allParameters[6] = allParameters[6] + sp;
                        }
                    }
                }

                flag++;

            }

            double[] r = new double[s.Length];
            r = median(T, QI, spe);
            double[] re = new double[s.Length];
            re = IQR(T, QI, spe);

            for (int k = 0; k < allParameters.Length; k++)
            {
                Console.WriteLine("Average " + s[k] + ": " + allParameters[k] / numbers);
                File.AppendAllText("Result/MOEAD_DE/LZ09_F1/finalQI" + ".txt", "Average " + s[k] + ": " + allParameters[k] / numbers + "\n");
                Console.WriteLine("Median " + s[k] + ": " + r[k]);
                File.AppendAllText("Result/MOEAD_DE/LZ09_F1/finalQI" + ".txt", "Median " + s[k] + ": " + r[k] + "\n");
                Console.WriteLine("IQR " + s[k] + ": " + re[k]);
                File.AppendAllText("Result/MOEAD_DE/LZ09_F1/finalQI" + ".txt", "IQR " + s[k] + ": " + re[k] + "\n");
            }

        }

        private void SD(double average, int n)
        {

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
                eps[l] = Q[4, 1];

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

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            MyAlgorithm();
        }

        private void Button_Click_1(object sender, RoutedEventArgs e)
        {
            calculateStatistics("Result/MOEAD_DE/LZ09_F1", "MOEAD_DE", 100);
        }
    }
}
