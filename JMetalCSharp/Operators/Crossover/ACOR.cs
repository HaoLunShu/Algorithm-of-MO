using JMetalCSharp.Core;
using JMetalCSharp.Encoding.SolutionType;
using JMetalCSharp.Utils;
using JMetalCSharp.Utils.Wrapper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JMetalCSharp.Operators.Crossover
{
    public class ACOR : Crossover
    {
        /// <summary>
		/// DEFAULT_CR defines a default CR (crossover operation control) value
		/// </summary>
		private static readonly double DEFAULT_ZELTA = 0.5;

        private double zelta;

        public ACOR(Dictionary<string, object> parameters)
            : base(parameters)
        {
            zelta = DEFAULT_ZELTA;
            Utils.Utils.GetDoubleValueFromParameter(parameters, "zelta", ref zelta);
        }

        private Solution DoACOr(Solution parent)
        {
            Solution current = new Solution(parent);

            Solution child;

            child = new Solution(current);

            XReal xCurrent = new XReal(current);
            XReal xChild = new XReal(child);

            int numberOfVariables = xCurrent.GetNumberOfDecisionVariables();

            for (int j = 0; j < numberOfVariables; j++)
            {
                double value;
                //    value = xParent2.GetValue(j) + f * (xParent0.GetValue(j) - xParent1.GetValue(j));

                double u1 = JMetalRandom.NextDouble(0, 1);
                double u2 = JMetalRandom.NextDouble(0, 1);
                double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
                value = xCurrent.GetValue(j) + current.stdDev[j] * randStdNormal;

                if (value < xChild.GetLowerBound(j))
                {
                    value = xChild.GetLowerBound(j);
                }
                if (value > xChild.GetUpperBound(j))
                {
                    value = xChild.GetUpperBound(j);
                }

                xChild.SetValue(j, value);
            }
            return child;
        }
    }
}
