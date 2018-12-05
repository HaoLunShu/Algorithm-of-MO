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
using System.Data.Objects;
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
using JMetalCSharp.Problems.Kursawe;
using JMetalCSharp.Problems.Schaffer;
using JMetalCSharp.Problems.Fonseca;

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

        Dictionary<string, string> setting = new Dictionary<string, string>();

        public void MyAlgorithm()
        {
            Problem problem = null; // The problem to solve
            Algorithm algorithm = null; // The algorithm to use
            Operator crossover = null; // Crossover operator
            Operator mutation = null; // Mutation operator
            Operator selection= null; // Selection operator

            Dictionary<string, object> parameters; // Operator parameters

            QualityIndicator indicators; // Object to get quality indicators

            string text;
            string pb = "";
            string st = "";
            string nov = "";
            string noo = "";
            string al = "";
            string ps = "";
            string me = "";
            string itn = "";
            string dad = "";
            string t = "";
            string delta = "";
            string nr = "";
            string co = "";
            string poc = "";
            string dioc = "";
            string cr = "";
            string f = "";
            string k = "";
            string zelta = "";
            string dev = "";
            string mu = "";
            //string pom = "";
            string diom = "";
            string s = "";
            string qi = "";
            string rt = "";
            System.IO.StreamReader sr = new System.IO.StreamReader(@"Data\Parameters\setting.txt");
            while ((text = sr.ReadLine()) != null)
            {
                string[] words = text.Split(' ');
                if (true == setting.ContainsKey(words[0]))
                    setting[words[0]] = words[1];
                else
                    setting.Add(words[0], words[1]);
            }
            sr.Close();

            if (true == setting.ContainsKey("problem"))
                pb = setting["problem"];
            if (true == setting.ContainsKey("solutionType"))
                st = setting["solutionType"];
            if (true == setting.ContainsKey("numberOfVariables"))
                nov = setting["numberOfVariables"];
            if (true == setting.ContainsKey("numberOfObjectives"))
                noo = setting["numberOfObjectives"];
            if (true == setting.ContainsKey("algorithm"))
                al = setting["algorithm"];
            if (true == setting.ContainsKey("populationSize"))
                ps = setting["populationSize"];
            if (true == setting.ContainsKey("maxEvaluations"))
                me = setting["maxEvaluations"];
            if (true == setting.ContainsKey("iterationsNumber"))
                itn = setting["iterationsNumber"];
            if (true == setting.ContainsKey("dataDirectory"))
                dad = setting["dataDirectory"];
            if (true == setting.ContainsKey("T"))
                t = setting["T"];
            if (true == setting.ContainsKey("delta"))
                delta = setting["delta"];
            if (true == setting.ContainsKey("nr"))
                nr = setting["nr"];
            if (true == setting.ContainsKey("Crossover"))
                co = setting["Crossover"];
            if (true == setting.ContainsKey("probabilityOfCrossover"))
                poc = setting["probabilityOfCrossover"];
            if (true == setting.ContainsKey("distributionIndexOfCrossover"))
                dioc = setting["distributionIndexOfCrossover"];
            if (true == setting.ContainsKey("CR"))
                cr = setting["CR"];
            if (true == setting.ContainsKey("F"))
                f = setting["F"];
            if (true == setting.ContainsKey("K"))
                k = setting["K"];
            if (true == setting.ContainsKey("zelta"))
                zelta = setting["zelta"];
            if (true == setting.ContainsKey("DEVariant"))
                dev = setting["DEVariant"];
            if (true == setting.ContainsKey("Mutation"))
                mu = setting["Mutation"];
            //if (true == setting.ContainsKey("probabilityOfMutation"))
                //pom = setting["probabilityOfMutation"];
            if (true == setting.ContainsKey("distributionIndexOfMutation"))
                diom = setting["distributionIndexOfMutation"];
            if (true == setting.ContainsKey("Selection"))
                s = setting["Selection"];
            if (true == setting.ContainsKey("QualityIndicator"))
                qi = setting["QualityIndicator"];
            if (true == setting.ContainsKey("repeatTimes"))
                rt = setting["repeatTimes"];
            switch (pb)
            {
                case "ZDT1":
                    problem = new ZDT1(st, int.Parse(nov));
                    break;
                case "ZDT2":
                    problem = new ZDT2(st, int.Parse(nov));
                    break;
                case "ZDT3":
                    problem = new ZDT3(st, int.Parse(nov));
                    break;
                case "ZDT4":
                    problem = new ZDT4(st, int.Parse(nov));
                    break;
                case "ZDT5":
                    problem = new ZDT5(st, int.Parse(nov));
                    break;
                case "ZDT6":
                    problem = new ZDT6(st, int.Parse(nov));
                    break;
                case "DTLZ1":
                    problem = new DTLZ1(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ2":
                    problem = new DTLZ2(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ3":
                    problem = new DTLZ3(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ4":
                    problem = new DTLZ4(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ5":
                    problem = new DTLZ5(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ6":
                    problem = new DTLZ6(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ7":
                    problem = new DTLZ7(st, int.Parse(nov), int.Parse(noo));
                    break;
                case "Fonseca":
                    problem = new Fonseca(st);
                    break;
                case "Kursawe":
                    problem = new Kursawe(st, int.Parse(nov));
                    break;
                case "LZ09_F1":
                    problem = new LZ09_F1(st);
                    break;
                case "LZ09_F2":
                    problem = new LZ09_F2(st);
                    break;
                case "LZ09_F3":
                    problem = new LZ09_F3(st);
                    break;
                case "LZ09_F4":
                    problem = new LZ09_F4(st);
                    break;
                case "LZ09_F5":
                    problem = new LZ09_F5(st);
                    break;
                case "LZ09_F6":
                    problem = new LZ09_F6(st);
                    break;
                case "LZ09_F7":
                    problem = new LZ09_F7(st);
                    break;
                case "LZ09_F8":
                    problem = new LZ09_F8(st);
                    break;
                case "LZ09_F9":
                    problem = new LZ09_F9(st);
                    break;
                case "Schaffer":
                    problem = new Schaffer(st);
                    break;
                default:
                    break;
            }

            string dirPath = "Result/" + al + "_" + co + "/" + pb + "_" + st;
            if (Directory.Exists(dirPath))
            {
                Console.WriteLine("The directory {0} already exists.", dirPath);
            }
            else
            {
                Directory.CreateDirectory(dirPath);
                Console.WriteLine("The directory {0} was created.", dirPath);
            }

            string filepath = dirPath + "/Parameter.txt";
            string[] line1 = { "numberOfVariables " + nov, "numberOfObjectives " + noo, "populationSize " + ps, "maxEvaluations " + me, "iterationsNumber " + itn };
            string[] line2 = { "T " + t, "delta " + delta, "nr " + nr };
            string[] line3 = { "probabilityOfCrossover " + poc, "distributionIndexOfCrossover " + dioc };
            string[] line4 = { "CR " + cr, "F " + f, "K " + k };
            string[] line5 = { "zelta " + zelta };
            string[] line6 = { "probabilityOfMutation " + 1.0 / problem.NumberOfVariables, "distributionIndexOfMutation " + diom, "" };
            File.AppendAllLines(filepath, line1);

            switch (al)
            {
                case "NSGAII":
                    algorithm = new JMetalCSharp.Metaheuristics.NSGAII.NSGAII(problem);
                    break;
                case "MOEAD":
                    algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);
                    algorithm.SetInputParameter("T", int.Parse(t));
                    algorithm.SetInputParameter("delta", double.Parse(delta));
                    algorithm.SetInputParameter("nr", int.Parse(nr));
                    File.AppendAllLines(filepath, line2);
                    break;
            }

            indicators = null;
            // Default problem
            //problem = new Kursawe("Real", 3);
            //problem = new Kursawe("BinaryReal", 3);
            //problem = new Water("Real");
            //problem = new ZDT3("ArrayReal", 30);
            //problem = new LZ09_F1("Real");
            //problem = new ConstrEx("Real");
            //problem = new DTLZ1("Real", 10, 3);
            //problem = new OKA2("Real") ;

            //algorithm = new JMetalCSharp.Metaheuristics.NSGAII.NSGAII(problem);
            //algorithm = new ssNSGAII(problem);
            //algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);

            // Algorithm parameters
            algorithm.SetInputParameter("populationSize", int.Parse(ps));
            algorithm.SetInputParameter("maxEvaluations", int.Parse(me));
            algorithm.SetInputParameter("iterationsNumber", int.Parse(itn));

            algorithm.SetInputParameter("dataDirectory", "Data/Parameters/Weight");

            //algorithm.SetInputParameter("finalSize", 300); // used by MOEAD_DRA



            // Mutation and Crossover for Real codification 
            parameters = new Dictionary<string, object>();
            switch(co)
            {
                case "SBXCrossover":
                    parameters.Add("probability", double.Parse(poc));
                    parameters.Add("distributionIndex", double.Parse(dioc));
                    crossover = CrossoverFactory.GetCrossoverOperator("SBXCrossover", parameters);
                    File.AppendAllLines(filepath, line3);
                    break;
                case "DifferentialEvolutionCrossover":
                    parameters.Add("CR", double.Parse(cr));
                    parameters.Add("F", double.Parse(f));
                    parameters.Add("K", double.Parse(k));
                    crossover = CrossoverFactory.GetCrossoverOperator("DifferentialEvolutionCrossover", parameters);
                    File.AppendAllLines(filepath, line4);
                    break;
                case "ACOR":
                    parameters.Add("zelta", double.Parse(zelta));
                    crossover = CrossoverFactory.GetCrossoverOperator("ACOR", parameters);
                    File.AppendAllLines(filepath, line5);
                    break;
            }

            parameters = new Dictionary<string, object>();
            parameters.Add("probability", 1.0 / problem.NumberOfVariables);
            parameters.Add("distributionIndex", double.Parse(diom));
            mutation = MutationFactory.GetMutationOperator("PolynomialMutation", parameters);
            File.AppendAllLines(filepath, line6);
            //File.AppendAllText(filepath, "");

            // Selection Operator 
            parameters = null;
            selection = SelectionFactory.GetSelectionOperator("BinaryTournament2", parameters);

            // Quality Indicators Operator
            //indicators = new QualityIndicator(problem, "DTLZ1.3D.pf");
            indicators = new QualityIndicator(problem, qi);
            //indicators = new QualityIndicator(problem, "LZ09_F1.pf");

            // Add the operators to the algorithm
            algorithm.AddOperator("crossover", crossover);
            algorithm.AddOperator("mutation", mutation);
            algorithm.AddOperator("selection", selection);

            // Add the indicator object to the algorithm
            algorithm.SetInputParameter("indicators", indicators);

            for (int i = 1; i <= int.Parse(rt); i++)
            {

                // Logger object and file to store log messages
                var logger = Logger.Log;

                var appenders = logger.Logger.Repository.GetAppenders();
                var fileAppender = appenders[0] as log4net.Appender.FileAppender;
                fileAppender.File = dirPath + "/" + al + i +".log";
                fileAppender.ActivateOptions();

                string filevar = dirPath + "/VAR" + i;
                string filefun = dirPath + "/FUN" + i;

                FileStream file = new FileStream(dirPath + "/" + al + "new" + i + ".txt", FileMode.OpenOrCreate, FileAccess.ReadWrite);
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
            }
        }

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
                    bool bt = text.Contains(s[0]);
                    if (bt)
                    {
                        int indextime = text.IndexOf(s[0]);
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
                eps[l] = Q[4, 1];

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
            string pb = "";
            string st = "";
            string co = "";
            string al = "";
            string rt = "";
            if (true == setting.ContainsKey("problem"))
                pb = setting["problem"];
            if (true == setting.ContainsKey("solutionType"))
                st = setting["solutionType"];
            if (true == setting.ContainsKey("algorithm"))
                al = setting["algorithm"];
            if (true == setting.ContainsKey("Crossover"))
                co = setting["Crossover"];
            if (true == setting.ContainsKey("repeatTimes"))
                rt = setting["repeatTimes"];
            calculateStatistics("Result/" + al + "_" + co + "/" + pb + "_" + st, al, int.Parse(rt));
        }
    }
}
