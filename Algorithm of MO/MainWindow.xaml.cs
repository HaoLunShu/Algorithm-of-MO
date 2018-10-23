using System;
using System.Collections.Generic;
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
                problem = new ZDT1("ArrayReal", 30);
            //problem = new ConstrEx("Real");
            //problem = new DTLZ2("Real", 10, 3);
            //problem = new OKA2("Real") ;

            algorithm = new JMetalCSharp.Metaheuristics.NSGAII.NSGAII(problem);
            //algorithm = new ssNSGAII(problem);
            //algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);

            // Algorithm parameters
            algorithm.SetInputParameter("populationSize", 100);
            algorithm.SetInputParameter("maxEvaluations", 25000);
            algorithm.SetInputParameter("iterationsNumber", 250);

            algorithm.SetInputParameter("dataDirectory", "Data/Parameters/Weight");

            //algorithm.SetInputParameter("finalSize", 300); // used by MOEAD_DRA

            //algorithm.SetInputParameter("T", 20);
            //algorithm.SetInputParameter("delta", 0.9);
            //algorithm.SetInputParameter("nr", 2);

            // Mutation and Crossover for Real codification 
            parameters = new Dictionary<string, object>();
            /*parameters.Add("CR", 1.0);
            parameters.Add("F", 0.5);
            crossover = CrossoverFactory.GetCrossoverOperator("DifferentialEvolutionCrossover", parameters);*/
            parameters.Add("probability", 1.0);
            parameters.Add("distributionIndex", 20.0);
            crossover = CrossoverFactory.GetCrossoverOperator("SBXCrossover", parameters);

            parameters = new Dictionary<string, object>();
            parameters.Add("probability", 1.0 / problem.NumberOfVariables);
            parameters.Add("distributionIndex", 20.0);
            mutation = MutationFactory.GetMutationOperator("PolynomialMutation", parameters);

            // Selection Operator 
            parameters = null;
            selection = SelectionFactory.GetSelectionOperator("BinaryTournament2", parameters);

            // Quality Indicators Operator
            indicators = new QualityIndicator(problem, "ZDT1.pf");

            // Add the operators to the algorithm
            algorithm.AddOperator("crossover", crossover);
            algorithm.AddOperator("mutation", mutation);
            algorithm.AddOperator("selection", selection);

            // Add the indicator object to the algorithm
            algorithm.SetInputParameter("indicators", indicators);

            for (int i = 1; i <= 30; i++)
            {

                // Logger object and file to store log messages
                var logger = Logger.Log;

                var appenders = logger.Logger.Repository.GetAppenders();
                var fileAppender = appenders[0] as log4net.Appender.FileAppender;
                fileAppender.File = "Result/NSGAII/ZDT1/NSGAII" + i +".log";
                fileAppender.ActivateOptions();

                string filevar = "Result/NSGAII/ZDT1/VAR" + i;
                string filefun = "Result/NSGAII/ZDT1/FUN" + i;

                // Execute the Algorithm
                long initTime = Environment.TickCount;
                SolutionSet population = algorithm.Execute();
                long estimatedTime = Environment.TickCount - initTime;

                // Result messages 
                logger.Info("Total execution time: " + estimatedTime + "ms");
                logger.Info("Variables values have been writen to file " + filevar);
                population.PrintVariablesToFile(filevar);
                logger.Info("Objectives values have been writen to file " + filefun);
                population.PrintObjectivesToFile(filefun);
                Console.WriteLine("Time: " + estimatedTime);
                Console.ReadLine();
                if (indicators != null)
                {
                    logger.Info("Quality indicators");
                    logger.Info("Hypervolume: " + indicators.GetHypervolume(population));
                    logger.Info("GD         : " + indicators.GetGD(population));
                    logger.Info("IGD        : " + indicators.GetIGD(population));
                    logger.Info("Spread     : " + indicators.GetSpread(population));
                    logger.Info("Epsilon    : " + indicators.GetEpsilon(population));

                    int evaluations = (int)algorithm.GetOutputParameter("evaluations");
                    logger.Info("Speed      : " + evaluations + " evaluations");
                }
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            MyAlgorithm();
        }
    }
}
