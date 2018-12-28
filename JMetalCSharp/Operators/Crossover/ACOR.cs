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
		private static readonly double DEFAULT_ZELTA = 0.85;

        private double zelta;

        private double[] randStdNormal;

        public ACOR(Dictionary<string, object> parameters)
            : base(parameters)
        {
            zelta = DEFAULT_ZELTA;
            Utils.Utils.GetDoubleValueFromParameter(parameters, "zelta", ref zelta);
        }

        /// <summary>
		/// Valid solution types to apply this operator
		/// </summary>
		private static readonly List<Type> VALID_TYPES = new List<Type>()
        {
            typeof(RealSolutionType),
            typeof(ArrayRealSolutionType)
        };

        private Solution DoACOr(Solution parent)
        {
            Solution current = new Solution(parent);

            Solution child;

            child = new Solution(current);

            XReal xCurrent = new XReal(current);
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
                value = xCurrent.GetValue(j) + zelta * xCurrent.GetStdDev(j) * randStdNormal[j];

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
            Solution parents = (Solution)obj;

            if (!(VALID_TYPES.Contains(parents.Type.GetType())))
            {
                Logger.Log.Error("Exception in " + this.GetType().FullName + ".Execute()");
                Console.WriteLine("Exception in " + this.GetType().FullName + ".Execute()");
                throw new Exception("Exception in " + this.GetType().FullName + ".Execute()");
            }

            Solution offSpring;
            offSpring = DoACOr(parents);

            return offSpring;
        }

        #endregion

        public object get()
        {
            return this.randStdNormal;
        }

    }
}
