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
		private static readonly double DEFAULT_ZETA = 0.85;

        private double zeta;

        private double[] randStdNormal;

        public ACOR(Dictionary<string, object> parameters)
            : base(parameters)
        {
            zeta = DEFAULT_ZETA;
            Utils.Utils.GetDoubleValueFromParameter(parameters, "zeta", ref zeta);
        }

        /// <summary>
		/// Valid solution types to apply this operator
		/// </summary>
		private static readonly List<Type> VALID_TYPES = new List<Type>()
        {
            typeof(RealSolutionType),
            typeof(ArrayRealSolutionType)
        };

        private Solution DoACOr(Solution parent1, Solution parent2)
        {
            Solution current = new Solution(parent1);

            Solution child;

            child = new Solution(current);

            Solution refer = new Solution(parent2);

            XReal xCurrent = new XReal(current);
            XReal xRefer = new XReal(refer);
            XReal xChild = new XReal(child);

            int numberOfVariables = xCurrent.GetNumberOfDecisionVariables();

            randStdNormal = new double[numberOfVariables];

            for (int j = 0; j < numberOfVariables; j++)
            {
                double value;
                //    value = xParent2.GetValue(j) + f * (xParent0.GetValue(j) - xParent1.GetValue(j));

                double u1 = JMetalRandom.NextDouble(0, 1);
                double u2 = JMetalRandom.NextDouble(0, 1);
                randStdNormal[j] = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
                value = xCurrent.GetValue(j) + zeta * xCurrent.GetStdDev(j) * randStdNormal[j];

                if (value < xChild.GetLowerBound(j))
                {
                    value = xChild.GetLowerBound(j);
                    //value = JMetalRandom.NextDouble(xChild.GetLowerBound(j), xChild.GetUpperBound(j));
                }
                if (value > xChild.GetUpperBound(j))
                {
                    value = xChild.GetUpperBound(j);
                    //value = JMetalRandom.NextDouble(xChild.GetLowerBound(j), xChild.GetUpperBound(j));
                }

                xChild.SetValue(j, value);
            }
            return child;
        }

        #region Public Overrides

        /// <summary>
        /// Executes the operation
        /// </summary>
        /// <param name="obj">An object containing one parent</param>
        /// <returns>An object containing the offSprings</returns>
        public override object Execute(object obj)
        {
            Solution[] parents = (Solution[])obj;

            if (parents.Length != 2)
            {
                Logger.Log.Error("Exception in " + this.GetType().FullName + ".Execute()");
                Console.WriteLine("Exception in " + this.GetType().FullName + ".Execute()");
                throw new Exception("Exception in " + this.GetType().FullName + ".Execute()");
            }

            if (!(VALID_TYPES.Contains(parents[0].Type.GetType())
                    && VALID_TYPES.Contains(parents[1].Type.GetType())))
            {
                Logger.Log.Error("Exception in " + this.GetType().FullName + ".Execute()");
                Console.WriteLine("Exception in " + this.GetType().FullName + ".Execute()");
                throw new Exception("Exception in " + this.GetType().FullName + ".Execute()");
            }

            Solution offSpring;
            offSpring = DoACOr(parents[0], parents[1]);

            return offSpring;
        }

        #endregion

        public object get()
        {
            return this.randStdNormal;
        }

    }
}