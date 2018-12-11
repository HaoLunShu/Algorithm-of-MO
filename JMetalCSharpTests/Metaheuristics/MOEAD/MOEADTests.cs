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
            Problem problem = new ZDT1("Real", 30);
            Algorithm algorithm = new JMetalCSharp.Metaheuristics.MOEAD.MOEAD(problem);
            algorithm.SetInputParameter("T", 20);
            algorithm.SetInputParameter("delta", 0.9);
            algorithm.SetInputParameter("nr", 2);
            algorithm.SetInputParameter("populationSize", 100);
            algorithm.SetInputParameter("iterationsNumber", 1);
            algorithm.SetInputParameter("dataDirectory", "Data/Parameters/Weight");
            Dictionary<string, object> parameters = new Dictionary<string, object>();
            parameters.Add("zelta", 0.85);
            Operator crossover = CrossoverFactory.GetCrossoverOperator("ACOR", parameters);
            algorithm.AddOperator("crossover", crossover);

            // Execute the Algorithm
            long initTime = Environment.TickCount;
            SolutionSet population = algorithm.Execute();
            long estimatedTime = Environment.TickCount - initTime;

            Assert.Fail();
        }
    }
}