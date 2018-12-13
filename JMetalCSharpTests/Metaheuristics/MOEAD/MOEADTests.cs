using Microsoft.VisualStudio.TestTools.UnitTesting;
using JMetalCSharp.Metaheuristics.MOEAD;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using JMetalCSharp.Problems.ZDT;
using JMetalCSharp.Operators.Crossover;
using JMetalCSharp.Core;

namespace JMetalCSharp.Metaheuristics.MOEAD.Tests
{
    [TestClass()]
    public class MOEADTests
    {
        [TestMethod()]
        public void ACOrSelectionTest()
        {
            Problem problem = new ZDT1("Real", 10);
            int populationSize = 10;
            int t = 2;
            SolutionSet population = new SolutionSet(populationSize);
            MOEAD moead = new MOEAD(problem);
            int p = 0;
            double[][] expectedLambda = new double[10][];
            double[][] dis = new double[10][];
            double[][] expectedneigh = new double[10][];
            double[] z = new double[2];

            for (int i = 0; i < populationSize; i++)
            {
                double a = 1.0 * i / (populationSize - 1);
                expectedLambda[i][0] = a;
                expectedLambda[i][1] = 1 - a;
            }

            for(int i = 0; i < populationSize; i++)
            {
                for(int j = 0; j < populationSize; j++)
                {
                    dis[i][j] = Math.Sqrt(Math.Pow((expectedLambda[i][0] - expectedLambda[j][0]), 2) + Math.Pow((expectedLambda[i][1] - expectedLambda[j][1]), 2));
                }
                expectedneigh[i] = MinSort(dis[i], t, populationSize);
            }

            for (int i = 0; i < populationSize; i++)
            {
                Solution newSolution = new Solution(problem);

                problem.Evaluate(newSolution);
                population.Add(newSolution);
            }

            for(int i = 0; i < 2; i++)
            {
                z[i] = 1.0e+30;
            }

            for(int i = 0; i < populationSize; i++)
            {
                Solution temp = population.Get(i);
                for(int j = 0; j < 2; j++)
                {
                    if(temp.Objective[j]  < z[j])
                    {
                        z[j] = temp.Objective[j];
                    }
                }
            }
            /*Algorithm algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);
            algorithm.SetInputParameter("T", 20);
            algorithm.SetInputParameter("delta", 0.5);
            algorithm.SetInputParameter("nr", 2);*/

            moead.set();

            for (int i = 0; i < populationSize; i++)
            {
                p = moead.ACOrSelection(i, 0);
            }

            Console.WriteLine(p);

            /*algorithm.SetInputParameter("populationSize", 100);
            algorithm.SetInputParameter("iterationsNumber", 1);
            algorithm.SetInputParameter("dataDirectory", "Data/Parameters/Weight");
            Dictionary<string, object> parameters = new Dictionary<string, object>();
            parameters.Add("zelta", 0.85);
            Operator crossover = CrossoverFactory.GetCrossoverOperator("ACOR", parameters);
            algorithm.AddOperator("crossover", crossover);*/

            Assert.Fail();
        }

        public double[] MinSort(double[] array, int m, int n)
        {
            for(int i = 0; i < m; i++)
            {
                for(int j = i + 1; j < n; j++)
                {
                    if(array[i] > array[j])
                    {
                        double temp = array[i];
                        array[i] = array[j];
                        array[j] = temp;
                    }
                }
            }

            return array;
        }
    }
}