using JMetalCSharp.Core;
using JMetalCSharp.Operators.Crossover;
using JMetalCSharp.Operators.Mutation;
using JMetalCSharp.Operators.Selection;
using JMetalCSharp.Problems.DTLZ;
using JMetalCSharp.Problems.Fonseca;
using JMetalCSharp.Problems.Kursawe;
using JMetalCSharp.Problems.LZ09;
using JMetalCSharp.Problems.Schaffer;
using JMetalCSharp.Problems.ZDT;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithm_of_MO
{
    class FileRead
    {
        Problem problem = null; // The problem to solve
        Algorithm algorithm = null; // The algorithm to use
        Operator crossover = null; // Crossover operator
        Operator crossover2 = null; // Crossover operator
        Operator mutation = null; // Mutation operator
        Operator selection = null; // Selection operator

        string text;
        public string Pb { get; private set; } = "";
        public string St { get; private set; } = "";
        string nov = "";
        string noo = "";
        public string Al { get; set; } = "";
        public string Ps { get; private set; } = "";
        public string Itn { get; private set; } = "";
        string dad = "";
        public string T { get; private set; } = "";
        public string Delta { get; private set; } = "";
        public string Nr { get; private set; } = "";
        public string Co { get; set; } = "";
        public string Co2 { get; set; } = "";
        public string Ra { get; set; } = "";
        string poc = "";
        string dioc = "";
        string cr = "";
        string f = "";
        string k = "";
        string zeta = "";
        string q = "";
        string mu = "";
        string diom = "";
        string s = "";
        public string Qi { get; private set; } = "";
        public string Rt { get; private set; } = "";
        public string DirPath { get; private set; } = "";

        public void fileread()
        {
            Dictionary<string, string> setting = new Dictionary<string, string>();
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
                Pb = setting["problem"];
            if (true == setting.ContainsKey("solutionType"))
                St = setting["solutionType"];
            if (true == setting.ContainsKey("numberOfVariables"))
                nov = setting["numberOfVariables"];
            if (true == setting.ContainsKey("numberOfObjectives"))
                noo = setting["numberOfObjectives"];
            if (true == setting.ContainsKey("algorithm"))
                Al = setting["algorithm"];
            if (true == setting.ContainsKey("populationSize"))
                Ps = setting["populationSize"];
            if (true == setting.ContainsKey("iterationsNumber"))
                Itn = setting["iterationsNumber"];
            if (true == setting.ContainsKey("dataDirectory"))
                dad = setting["dataDirectory"];
            if (true == setting.ContainsKey("T"))
                T = setting["T"];
            if (true == setting.ContainsKey("delta"))
                Delta = setting["delta"];
            if (true == setting.ContainsKey("nr"))
                Nr = setting["nr"];
            if (true == setting.ContainsKey("Crossover"))
                Co = setting["Crossover"];
            if (true == setting.ContainsKey("Crossover2"))
                Co2 = setting["Crossover2"];
            if (true == setting.ContainsKey("Ratio"))
                Ra = setting["Ratio"];
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
            if (true == setting.ContainsKey("zeta"))
                zeta = setting["zeta"];
            if (true == setting.ContainsKey("q"))
                q = setting["q"];
            if (true == setting.ContainsKey("Mutation"))
                mu = setting["Mutation"];
            else mu = null;
            if (true == setting.ContainsKey("distributionIndexOfMutation"))
                diom = setting["distributionIndexOfMutation"];
            if (true == setting.ContainsKey("Selection"))
                s = setting["Selection"];
            else s = null;
            if (true == setting.ContainsKey("QualityIndicator"))
                Qi = setting["QualityIndicator"];
            if (true == setting.ContainsKey("repeatTimes"))
                Rt = setting["repeatTimes"];
        }

        public Problem GetProblem()
        {
            switch (Pb)
            {
                case "ZDT1":
                    problem = new ZDT1(St, int.Parse(nov));
                    break;
                case "ZDT2":
                    problem = new ZDT2(St, int.Parse(nov));
                    break;
                case "ZDT3":
                    problem = new ZDT3(St, int.Parse(nov));
                    break;
                case "ZDT4":
                    problem = new ZDT4(St, int.Parse(nov));
                    break;
                case "ZDT5":
                    problem = new ZDT5(St, int.Parse(nov));
                    break;
                case "ZDT6":
                    problem = new ZDT6(St, int.Parse(nov));
                    break;
                case "DTLZ1":
                    problem = new DTLZ1(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ2":
                    problem = new DTLZ2(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ3":
                    problem = new DTLZ3(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ4":
                    problem = new DTLZ4(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ5":
                    problem = new DTLZ5(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ6":
                    problem = new DTLZ6(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "DTLZ7":
                    problem = new DTLZ7(St, int.Parse(nov), int.Parse(noo));
                    break;
                case "Fonseca":
                    problem = new Fonseca(St);
                    break;
                case "Kursawe":
                    problem = new Kursawe(St, int.Parse(nov));
                    break;
                case "LZ09_F1":
                    problem = new LZ09_F1(St);
                    break;
                case "LZ09_F2":
                    problem = new LZ09_F2(St);
                    break;
                case "LZ09_F3":
                    problem = new LZ09_F3(St);
                    break;
                case "LZ09_F4":
                    problem = new LZ09_F4(St);
                    break;
                case "LZ09_F5":
                    problem = new LZ09_F5(St);
                    break;
                case "LZ09_F6":
                    problem = new LZ09_F6(St);
                    break;
                case "LZ09_F7":
                    problem = new LZ09_F7(St);
                    break;
                case "LZ09_F8":
                    problem = new LZ09_F8(St);
                    break;
                case "LZ09_F9":
                    problem = new LZ09_F9(St);
                    break;
                case "Schaffer":
                    problem = new Schaffer(St);
                    break;
                default:
                    break;
            }

            return problem;
        }

        public Algorithm GetAlgorithm()
        {
            if(Co2 != "")
                DirPath = "Result/" + Al + "_" + Co + "_" + Co2 + "/" + Pb + "_" + St;
            else
                DirPath = "Result/" + Al + "_" + Co + "/" + Pb + "_" + St;
            if (Directory.Exists(DirPath))
            {
                Console.WriteLine("The directory {0} already exists.", DirPath);
            }
            else
            {
                Directory.CreateDirectory(DirPath);
                Console.WriteLine("The directory {0} was created.", DirPath);
            }

            string filepath = DirPath + "/Parameter.txt";
            string[] line1 = { "numberOfVariables " + nov, "numberOfObjectives " + noo, "populationSize " + Ps, "iterationsNumber " + Itn };
            string[] line2 = { "T " + T, "delta " + Delta, "nr " + Nr };
            File.AppendAllLines(filepath, line1);

            switch (Al)
            {
                case "NSGAII":
                    algorithm = new JMetalCSharp.Metaheuristics.NSGAII.NSGAII(problem);
                    break;
                case "MOEAD":
                    algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);
                    algorithm.SetInputParameter("T", int.Parse(T));
                    algorithm.SetInputParameter("delta", double.Parse(Delta));
                    algorithm.SetInputParameter("nr", int.Parse(Nr));
                    File.AppendAllLines(filepath, line2);
                    break;
                default:
                    break;
            }
            algorithm.SetInputParameter("populationSize", int.Parse(Ps));
            algorithm.SetInputParameter("q", double.Parse(q));
            algorithm.SetInputParameter("iterationsNumber", int.Parse(Itn));

            algorithm.SetInputParameter("dataDirectory", "Data/Parameters/Weight");

            return algorithm;
        }

        public Operator GetCrossover()
        {
            if (Co2 != "")
                DirPath = "Result/" + Al + "_" + Co + "_" + Co2 + "/" + Pb + "_" + St;
            else
                DirPath = "Result/" + Al + "_" + Co + "/" + Pb + "_" + St;
            if (Directory.Exists(DirPath))
            {
                Console.WriteLine("The directory {0} already exists.", DirPath);
            }
            else
            {
                Directory.CreateDirectory(DirPath);
                Console.WriteLine("The directory {0} was created.", DirPath);
            }

            string filepath = DirPath + "/Parameter.txt";
            string[] line3 = { "probabilityOfCrossover " + poc, "distributionIndexOfCrossover " + dioc };
            string[] line4 = { "CR " + cr, "F " + f, "K " + k };
            string[] line6 = { "probabilityOfMutation " + 1.0 / problem.NumberOfVariables, "distributionIndexOfMutation " + diom, "" };
            string[] line5 = { "zeta " + zeta, "q " + q };

            Dictionary<string, object>  parameters = new Dictionary<string, object>();

            switch (Co)
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
                    parameters.Add("zeta", double.Parse(zeta));
                    crossover = CrossoverFactory.GetCrossoverOperator("ACOR", parameters);
                    File.AppendAllLines(filepath, line5);
                    break;
                case null:
                    parameters.Add("probability", double.Parse(poc));
                    parameters.Add("distributionIndex", double.Parse(dioc));
                    parameters.Add("zeta", double.Parse(zeta));
                    crossover = CrossoverFactory.GetCrossoverOperator("null", parameters);
                    File.AppendAllLines(filepath, line3);
                    File.AppendAllLines(filepath, line5);
                    break;
                default:
                    break;
            }

            return crossover;
        }

        public Operator GetCrossover2()
        {
            Dictionary<string, object> parameters = new Dictionary<string, object>();

            switch (Co2)
            {
                case "SBXCrossover":
                    parameters.Add("probability", double.Parse(poc));
                    parameters.Add("distributionIndex", double.Parse(dioc));
                    crossover2 = CrossoverFactory.GetCrossoverOperator("SBXCrossover", parameters);
                    break;
                case "DifferentialEvolutionCrossover":
                    parameters.Add("CR", double.Parse(cr));
                    parameters.Add("F", double.Parse(f));
                    parameters.Add("K", double.Parse(k));
                    crossover2 = CrossoverFactory.GetCrossoverOperator("DifferentialEvolutionCrossover", parameters);
                    break;
                case "ACOR":
                    parameters.Add("zeta", double.Parse(zeta));
                    crossover2 = CrossoverFactory.GetCrossoverOperator("ACOR", parameters);
                    break;
                case null:
                    parameters.Add("probability", double.Parse(poc));
                    parameters.Add("distributionIndex", double.Parse(dioc));
                    parameters.Add("zeta", double.Parse(zeta));
                    crossover2 = CrossoverFactory.GetCrossoverOperator("null", parameters);
                    break;
                default:
                    break;
            }

            return crossover2;
        }

        public Operator GetMutation()
        {
            if (Co2 != "")
                DirPath = "Result/" + Al + "_" + Co + "_" + Co2 + "/" + Pb + "_" + St;
            else
                DirPath = "Result/" + Al + "_" + Co + "/" + Pb + "_" + St;
            if (Directory.Exists(DirPath))
            {
                Console.WriteLine("The directory {0} already exists.", DirPath);
            }
            else
            {
                Directory.CreateDirectory(DirPath);
                Console.WriteLine("The directory {0} was created.", DirPath);
            }

            string filepath = DirPath + "/Parameter.txt";
            string[] line6 = { "probabilityOfMutation " + 1.0 / problem.NumberOfVariables, "distributionIndexOfMutation " + diom, "" };

            Dictionary<string, object> parameters = new Dictionary<string, object>();

            switch (mu)
            {
                case "PolynomialMutation":
                    parameters = new Dictionary<string, object>();
                    parameters.Add("probability", 1.0 / problem.NumberOfVariables);
                    //parameters.Add("probability", 0.2);
                    parameters.Add("distributionIndex", double.Parse(diom));
                    mutation = MutationFactory.GetMutationOperator("PolynomialMutation", parameters);
                    File.AppendAllLines(filepath, line6);
                    break;
                case "DynamicPolynomialMutation":
                    parameters = new Dictionary<string, object>();
                    parameters.Add("probability", 1.0);
                    parameters.Add("distributionIndex", double.Parse(diom));
                    mutation = MutationFactory.GetMutationOperator("DynamicPolynomialMutation", parameters);
                    File.AppendAllLines(filepath, line6);
                    break;
                case "BoxMuller":
                    parameters.Add("zeta", double.Parse(zeta));
                    parameters.Add("probability", 1.0 / problem.NumberOfVariables);
                    crossover = MutationFactory.GetMutationOperator("BoxMuller", parameters);
                    File.AppendAllLines(filepath, line6);
                    break;
                case null:
                    break;
                default:
                    break;
            }

            return mutation;
        }

        public Operator GetSelection()
        {
            Dictionary<string, object> parameters = new Dictionary<string, object>();

            switch (s)
            {
                case "BinaryTournament2":
                    parameters = null;
                    selection = SelectionFactory.GetSelectionOperator("BinaryTournament2", parameters);
                    break;
                default:
                    break;
            }

            return selection;
        }
    }
}
