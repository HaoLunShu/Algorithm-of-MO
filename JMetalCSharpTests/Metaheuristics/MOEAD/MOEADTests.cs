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
using JMetalCSharp.Encoding.Variable;

namespace JMetalCSharp.Metaheuristics.MOEAD.Tests
{
    [TestClass()]
    public class MOEADTests
    {
        [TestMethod()]
        public void ACOrSelectionTest()
        {
            Problem problem = new ZDT1("Real", 10);
            int populationSize = 11;
            SolutionSet population = new SolutionSet(populationSize);
            MOEAD moead = new MOEAD(problem);
            //Operator crossover = null;
            ACOR crossover = null;
            Dictionary<string, object> parameters = new Dictionary<string, object>();
            parameters.Add("zelta", 0.85);
            crossover = new ACOR(parameters);
            int p = 0;
            int[] expectedp = new int[11] { 1, 1, 2, 3, 4, 4, 5, 6, 8, 8, 9 };
            double[][] expectedLambda = new double[11][] { new double[] { 0, 1 }, new double[] { 0.1, 0.9 }, new double[] { 0.2, 0.8 }, new double[] { 0.3, 0.7 }, new double[] { 0.4, 0.6 }, new double[] { 0.5, 0.5 }, new double[] { 0.6, 0.4 }, new double[] { 0.7, 0.3 }, new double[] { 0.8, 0.2 }, new double[] { 0.9, 0.1 }, new double[] { 1, 0 } };
            int[][] expectedneigh = new int[11][] { new int[] { 0, 1 }, new int[] { 1, 0 }, new int[] { 2, 1 }, new int[] { 3, 2 }, new int[] { 4, 5 }, new int[] { 5, 4 }, new int[] { 6, 5 }, new int[] { 7, 6 }, new int[] { 8, 9 }, new int[] { 9, 8 }, new int[] { 10, 9 } };
            //double[] expectedz = new double[2] { 0, 1 };
            double[] expectedz = new double[2] { 1.024e-07, 0.974113029 };
            double k = 1;
            double[] expectedstdDev = new double[11] { 0.8, 0.8, 0.16, 0.032, 0.00128, 0.00128, 0.000256, 0.0000512, 0.000002048, 0.000002048, 4.096e-07 };
            double[][] expectedpro = new double[11][] { new double[] { 0.155241307, 0.844758693 }, new double[] { 0.844758693, 0.155241307 }, new double[] { 0.875915826, 0.124084174 }, new double[] { 0.966565547, 0.033434453 }, new double[] { 0.910514895, 0.089485105 }, new double[] { 0.128479172, 0.871520828 }, new double[] { 0.37023011, 0.62976989 }, new double[] { 0.4515275, 0.5484725 }, new double[] { 0.508750379, 0.491249621 }, new double[] { 0.491249621, 0.508750379 }, new double[] { 0.496138879, 0.503861121 } };
            double[][] expectedVariable = new double[11][] { new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10] };
            double[][] expectedvariable = new double[11][] { new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10], new double[10] };

            for (int i = 0; i < populationSize; i++)
            {
                //double k = (double) i / 10.0;
                k = Math.Pow(0.2 , i);
                Variable[] variables = new Variable[10];
                for(int j = 0; j < variables.Length; j++)
                {
                    variables[j] = new Real(0, 1, k);
                    expectedVariable[i][j] = k;
                }
                Solution newSolution = new Solution(problem, variables);

                problem.Evaluate(newSolution);
                population.Add(i , newSolution);
            }

            moead.AutoSet(population);

            int[][] n = (int[][]) moead.Get("neighborhood");
            double[][] l = (double[][]) moead.Get("lambda");
            double[] z = (double[]) moead.Get("z");

            for(int i = 0; i < 11; i++)
            {
                for(int j = 0; j < 2; j++)
                {
                    Assert.AreEqual(expectedneigh[i][j], n[i][j]);
                    Assert.AreEqual(expectedLambda[i][j], l[i][j], 0.0001);
                    Assert.AreEqual(expectedz[j], z[j], 0.0001);
                }
            }
            
            for (int i = 0; i < populationSize; i++)
            {
                double[] rand = new double[10];
                p = moead.ACOrSelection(i, 1);
                //double[] pro = (double[])moead.Get("probability");
                Assert.AreEqual(expectedp[i], p);
                for(int j = 0; j < expectedneigh[i].Length; j++)
                {
                    Solution child = (Solution) crossover.Execute(population.Get(n[i][j]));
                    rand = (double[]) crossover.get();
                    //Assert.AreEqual(0.1, moead.Get(i, j), 0.00001);
                    Assert.AreEqual(expectedstdDev[i], moead.Get(i, j), 0.000001);
                    //Assert.AreEqual(expectedpro[i][j], pro[j], 0.0001);
                    for(int m = 0; m < 10; m++)
                    {
                        expectedvariable[n[i][j]][m] = expectedVariable[n[i][j]][m] + 0.85 * expectedstdDev[n[i][j]] * rand[m];
                        if (expectedvariable[n[i][j]][m] > 1)
                            expectedvariable[n[i][j]][m] = 1;
                        if (expectedvariable[n[i][j]][m] < 0)
                            expectedvariable[n[i][j]][m] = 0;
                        Assert.AreEqual(expectedvariable[n[i][j]][m], Double.Parse(child.Variable[m].ToString()), 0.000001);
                    }
                }
            }
        }
    }
}