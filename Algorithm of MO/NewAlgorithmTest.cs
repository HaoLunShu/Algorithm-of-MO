using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using JMetalCSharp;
using JMetalCSharp.Core;
using JMetalCSharp.Encoding.SolutionType;
using JMetalCSharp.Utils;
using JMetalCSharp.Utils.Wrapper;
using JMetalCSharp.Operators.Crossover;
using JMetalCSharp.Operators.Mutation;
using JMetalCSharp.Operators.Selection;
using JMetalCSharp.Problems;
using JMetalCSharp.Problems.DTLZ;
using JMetalCSharp.Problems.ZDT;
using JMetalCSharp.Problems.LZ09;
using JMetalCSharp.QualityIndicator;
using JMetalCSharp.Metaheuristics.MOEAD;
using System.IO;

namespace Algorithm_of_MO
{
    class NewAlgorithmTest
    {
        #region Private attributes

        private FileRead fileRead;

        private int populationSize;

        Problem problem;

        /// <summary>
        /// Store the population
        /// </summary>
        private SolutionSet population;

        /// <summary>
        /// Z vector (ideal point)
        /// </summary>
        private double[] z;

        /// <summary>
		/// Z nad vector (ideal point)
		/// </summary>
		//private double[] znad;

        /// <summary>
        /// Lambda vectors
        /// </summary>
        private double[][] lambda;

        /// <summary>
        /// neighbour size
        /// </summary>
        private int t;

        /// <summary>
        /// Neighborhood
        /// </summary>
        private int[][] neighborhood;

        /// <summary>
        /// probability that parent solutions are selected from neighbourhood
        /// </summary>
        private double delta;

        /// <summary>
        /// probability that parent solutions are selected from neighbourhood
        /// </summary>
        private double gamma;

        /// <summary>
        /// maximal number of solutions replaced by each child solution
        /// </summary>
        private int nr;

        private Solution[] indArray;

        private int evaluations;

        private int iterationsNumber;
        private int iteration;

        private Operator[] crossover = new Operator[3];
        private Operator mutation;

        private string dataDirectory;

        #endregion

        public NewAlgorithmTest()
        {
            fileRead = new FileRead();
            fileRead.fileread();
            problem = fileRead.GetProblem();

            crossover = fileRead.GetCrossover();
            mutation = fileRead.GetMutation();
        }

        #region Main Functions

        public SolutionSet Mix()
        {
            QualityIndicator indicators = new QualityIndicator(problem, fileRead.Qi); // QualityIndicator object
            int requiredEvaluations = 0; // Use in the example of use of the
                                         // indicators object (see below)

            evaluations = 0;

            iteration = 0;
            populationSize = int.Parse(fileRead.Ps);
            iterationsNumber = int.Parse(fileRead.Itn);
            dataDirectory = "Data/Parameters/Weight";

            int allR = int.Parse(fileRead.DERa) + int.Parse(fileRead.SBXRa) + int.Parse(fileRead.ACORa);
            double ACOR = (double) int.Parse(fileRead.ACORa) / allR;
            double DER = (double) int.Parse(fileRead.DERa) / allR;



            Logger.Log.Info("POPSIZE: " + populationSize);
            Console.WriteLine("POPSIZE: " + populationSize);

            population = new SolutionSet(populationSize);
            indArray = new Solution[problem.NumberOfObjectives];

            t = int.Parse(fileRead.T);
            nr = int.Parse(fileRead.Nr);
            delta = double.Parse(fileRead.Delta);
            gamma = double.Parse(fileRead.Gamma);

            neighborhood = new int[populationSize][];
            for (int i = 0; i < populationSize; i++)
            {
                neighborhood[i] = new int[t];
            }

            z = new double[problem.NumberOfObjectives];
            //znad = new double[Problem.NumberOfObjectives];

            lambda = new double[populationSize][];
            for (int i = 0; i < populationSize; i++)
            {
                lambda[i] = new double[problem.NumberOfObjectives];
            }

            /**/string dir = "Result/" + fileRead.Al;
            if (fileRead.DERa != "0")
                dir = dir + "_DifferentialEvolutionCrossover";
            if (fileRead.SBXRa != "0")
                dir = dir + "_SBXCrossover";
            if (fileRead.ACORa != "0")
                dir = dir + "_ACOR";
            dir = dir + "/" + fileRead.Pb + "_" + fileRead.St + "/Record/DE(" + int.Parse(fileRead.DERa) + ")+SBX(" + int.Parse(fileRead.SBXRa) + ")+ACOR(" + int.Parse(fileRead.ACORa) + ")_NoMutation/nr=" + fileRead.Nr + "/ACOR+SBX";
            if (Directory.Exists(dir))
            {
                Console.WriteLine("The directory {0} already exists.", dir);
            }
            else
            {
                Directory.CreateDirectory(dir);
                Console.WriteLine("The directory {0} was created.", dir);
            }

            //Step 1. Initialization
            //Step 1.1 Compute euclidean distances between weight vectors and find T
            InitUniformWeight();

            InitNeighborhood();

            //Step 1.2 Initialize population
            InitPoputalion();

            //Step 1.3 Initizlize z
            InitIdealPoint();

            //Step 2 Update
            for(int a = 1; a <= iterationsNumber; a++)
            {
                int[] permutation = new int[populationSize];
                JMetalCSharp.Metaheuristics.MOEAD.Utils.RandomPermutation(permutation, populationSize);

                Solution[] parents = new Solution[2];
                int t = 0;

                if(a <= (ACOR) * iterationsNumber)
                {
                    //ACOR
                    for (int i = 0; i < populationSize; i++)
                    {

                        int n = permutation[i];
                        // or int n = i;


                        int type;
                        double rnd = JMetalRandom.NextDouble();

                        // STEP 2.1. ACOR selection based on probability
                        if (rnd < gamma) // if (rnd < realb)    
                        {
                            type = 1;   // minmum
                            //parents[0] = population.Get(ACOrSelection2(n, type, pro_T));
                        }
                        else
                        {
                            type = 2;   // whole neighborhood probability
                            //parents[0] = population.Get(ACOrSelection2(n, type, pro_A));
                        }
                        GetStdDev(neighborhood);
                        //GetStdDev1(neighborhood, type);

                        //List<int> p = new List<int>();
                        //MatingSelection(p, n, 1, type);

                        // STEP 2.2. Reproduction
                        Solution child;

                        parents[0] = population.Get(ACOrSelection(n, type));
                        parents[1] = population.Get(n);
                        //parents[0] = population.Get(p[0]);

                        // Apply ACOR crossover 
                        child = (Solution)crossover[2].Execute(parents);

                        child.NumberofReplace = t;

                        // Apply mutation
                        // mutation.Execute(child);

                        // Evaluation
                        problem.Evaluate(child);

                        evaluations++;

                        // STEP 2.3. Repair. Not necessary

                        // STEP 2.4. Update z_
                        UpdateReference(child);

                        // STEP 2.5. Update of solutions
                        t = UpdateProblemWithReplace(child, n, 1);
                    }
                }
                else if(a >= (1 - DER) * iterationsNumber)
                {
                    //DE
                    for (int i = 0; i < populationSize; i++)
                    {
                        int n = permutation[i]; // or int n = i;

                        int type;
                        double rnd = JMetalRandom.NextDouble();

                        // STEP 2.1. Mating selection based on probability
                        if (rnd < delta) // if (rnd < realb)    
                        {
                            type = 1;   // neighborhood
                        }
                        else
                        {
                            type = 2;   // whole population
                        }
                        List<int> p = new List<int>();
                        MatingSelection(p, n, 2, type);

                        // STEP 2.2. Reproduction
                        Solution child;
                        Solution[] parent = new Solution[3];

                        parent[0] = population.Get(p[0]);
                        parent[1] = population.Get(p[1]);

                        parent[2] = population.Get(n);

                        // Apply DE crossover 
                        child = (Solution)crossover[0].Execute(new object[] { population.Get(n), parent });

                        // Apply mutation
                        mutation.Execute(child);

                        // Evaluation
                        problem.Evaluate(child);

                        evaluations++;

                        // STEP 2.3. Repair. Not necessary

                        // STEP 2.4. Update z_
                        UpdateReference(child);

                        // STEP 2.5. Update of solutions
                        UpdateProblem(child, n, type);
                    }
                }
                else
                {
                    //SBX
                    //Create the offSpring solutionSet
                    SolutionSet offspringPopulation = new SolutionSet(populationSize);
                    for (int i = 0; i < (populationSize / 2); i++)
                    {
                        int n = permutation[i]; // or int n = i;

                        int type;
                        double rnd = JMetalRandom.NextDouble();

                        // STEP 2.1. Mating selection based on probability
                        if (rnd < delta) // if (rnd < realb)    
                        {
                            type = 1;   // neighborhood
                        }
                        else
                        {
                            type = 2;   // whole population
                        }
                        List<int> p = new List<int>();
                        MatingSelection(p, n, 2, type);

                        parents[0] = population.Get(p[0]);
                        parents[1] = population.Get(p[1]);

                        //obtain parents
                        Solution[] offSpring = (Solution[])crossover[1].Execute(parents);
                        //Solution child;
                        mutation.Execute(offSpring[0]);
                        mutation.Execute(offSpring[1]);
                        /*if(rnd < 0.5)
                        {
                            child = offSpring[0];
                        }
                        else
                        {
                            child = offSpring[1];
                        }*/
                        problem.Evaluate(offSpring[0]);
                        problem.Evaluate(offSpring[1]);
                        problem.EvaluateConstraints(offSpring[0]);
                        problem.EvaluateConstraints(offSpring[1]);
                        offspringPopulation.Add(offSpring[0]);
                        offspringPopulation.Add(offSpring[1]);
                        evaluations += 2;

                        // STEP 2.3. Repair. Not necessary

                        // STEP 2.4. Update z_
                        UpdateReference(offSpring[0]);
                        UpdateReference(offSpring[1]);

                        // STEP 2.5. Update of solutions
                        UpdateProblem(offSpring[0], n, type);
                        UpdateProblem(offSpring[1], n, type);
                    }
                }

                /**/string filevar = dir + "/VAR" + iteration;
                string filefun = dir + "/FUN" + iteration;
                population.PrintVariablesToFile(filevar);
                population.PrintObjectivesToFile(filefun);

                iteration++;

                if ((indicators != null) && (requiredEvaluations == 0))
                {
                    double HV = indicators.GetHypervolume(population);
                    if (HV >= (0.98 * indicators.TrueParetoFrontHypervolume))
                    {
                        requiredEvaluations = evaluations;
                    }
                }

            }

            Logger.Log.Info("ITERATION: " + iteration);
            Console.WriteLine("ITERATION: " + iteration);

            SolutionSet Result = population;

            //return population;

            // Return the first non-dominated front
            Ranking rank = new Ranking(population);

            //SolutionSet Result = rank.GetSubfront(0);

            return Result;
        }

        #endregion

        #region Private Methods

        /// <summary>
        /// InitUniformWeight
        /// </summary>
        private void InitUniformWeight()
        {
            if ((problem.NumberOfObjectives == 2) && (populationSize <= 300))
            {
                for (int n = 0; n < populationSize; n++)
                {
                    double a = 1.0 * n / (populationSize - 1);
                    lambda[n][0] = a;
                    lambda[n][1] = 1 - a;
                }
            }
            else
            {
                string dataFileName;
                dataFileName = "W" + problem.NumberOfObjectives + "D_" +
                  populationSize + ".dat";

                try
                {
                    // Open the file
                    using (StreamReader reader = new StreamReader(dataDirectory + "/" + dataFileName))
                    {

                        int numberOfObjectives = 0;
                        int i = 0;
                        int j = 0;
                        string aux = reader.ReadLine();
                        while (aux != null)
                        {
                            string[] st = aux.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            j = 0;
                            numberOfObjectives = st.Length;

                            foreach (string s in st)
                            {
                                double value = JMetalCSharp.Utils.Utils.ParseDoubleInvariant(s);
                                lambda[i][j] = value;
                                j++;
                            }
                            aux = reader.ReadLine();
                            i++;
                        }
                    }
                }
                catch (Exception ex)
                {
                    Logger.Log.Error("InitUniformWeight: failed when reading for file: " + dataDirectory + "/" + dataFileName, ex);
                    Console.WriteLine("InitUniformWeight: failed when reading for file: " + dataDirectory + "/" + dataFileName);
                }
            }
        }

        private void InitNeighborhood()
        {
            double[] x = new double[populationSize];
            int[] idx = new int[populationSize];

            for (int i = 0; i < populationSize; i++)
            {
                // calculate the distances based on weight vectors
                for (int j = 0; j < populationSize; j++)
                {
                    x[j] = JMetalCSharp.Metaheuristics.MOEAD.Utils.DistVector(lambda[i], lambda[j]);
                    idx[j] = j;
                } // for

                // find 'niche' nearest neighboring subproblems
                JMetalCSharp.Metaheuristics.MOEAD.Utils.MinFastSort(x, idx, populationSize, t);

                Array.Copy(idx, 0, neighborhood[i], 0, t);
            } // for
        }

        private void InitPoputalion()
        {
            for (int i = 0; i < populationSize; i++)
            {
                Solution newSolution = new Solution(problem);

                problem.Evaluate(newSolution);
                evaluations++;
                population.Add(newSolution);
            }
        }

        private void InitIdealPoint()
        {
            for (int i = 0; i < problem.NumberOfObjectives; i++)
            {
                z[i] = 1.0e+30;
                //znad[i] = 1.0e-30;
                indArray[i] = new Solution(problem);
                problem.Evaluate(indArray[i]);
                evaluations++;
            } // for

            for (int i = 0; i < populationSize; i++)
            {
                UpdateReference(population.Get(i));
            }
        }

        private void UpdateReference(Solution individual)
        {
            for (int n = 0; n < problem.NumberOfObjectives; n++)
            {
                if (individual.Objective[n] < z[n])
                {
                    z[n] = individual.Objective[n];

                    indArray[n] = individual;
                }
                /*if (individual.Objective[n] > znad[n])
                {
                    znad[n] = individual.Objective[n];
                }*/
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="list">the set of the indexes of selected mating parents</param>
        /// <param name="cid">the id of current subproblem</param>
        /// <param name="size">the number of selected mating parents</param>
        /// <param name="type">1 - neighborhood; otherwise - whole population</param>
        private void MatingSelection(List<int> list, int cid, int size, int type)
        {
            int ss;
            int r;
            int p;

            ss = neighborhood[cid].Length;
            while (list.Count < size)
            {
                if (type == 1)
                {
                    r = JMetalRandom.Next(0, ss - 1);
                    p = neighborhood[cid][r];
                }
                else
                {
                    p = JMetalRandom.Next(0, populationSize - 1);
                }
                bool flag = true;
                for (int i = 0; i < list.Count; i++)
                {
                    if (list[i] == p) // p is in the list
                    {
                        flag = false;
                        break;
                    }
                }

                if (flag)
                {
                    list.Add(p);
                }
            }
        }

        //double[] pro = new double[2];

        /// <summary>
		/// 
		/// </summary>
		/// <param name="cid">the id of current subproblem</param>
		/// <param name="type">1 - minimum; otherwise - whole neighborhood probability</param>
        public int ACOrSelection(int cid, int type)
        {
            int ss;
            ss = neighborhood[cid].Length;
            double r;
            int p = neighborhood[cid][0];
            Solution[] parents = new Solution[ss];
            double[] fitness = new double[ss];
            int indexOfmin = 0;
            double sum = 0;
            double[] pro = new double[ss];
            double a1 = 0;
            double a2 = 0;

            if (type == 1)
            {
                for (int i = 0; i < ss; i++)
                {
                    parents[i] = population.Get(neighborhood[cid][i]);
                    fitness[i] = FitnessFunction(parents[i], lambda[cid]);
                    if (fitness[i] < FitnessFunction(population.Get(p), lambda[cid]))
                    {
                        indexOfmin = i;
                    }
                    p = neighborhood[cid][indexOfmin];
                }
            }
            else
            {
                for (int i = 0; i < ss; i++)
                {
                    parents[i] = population.Get(neighborhood[cid][i]);
                    fitness[i] = 1 / FitnessFunction(parents[i], lambda[cid]);
                    sum = sum + fitness[i];
                }
                for (int j = 0; j < ss; j++)
                {
                    pro[j] = fitness[j] / sum;
                }
                r = JMetalRandom.NextDouble();
                for (int k = 0; k < pro.Length; k++)
                {
                    a2 = a2 + pro[k];
                    if (r < a2 && r >= a1)
                    {
                        p = neighborhood[cid][k];
                        break;
                    }
                    a1 = a1 + pro[k];
                }
            }

            return p;
        }

        /// <summary>
		/// 
		/// </summary>
		/// <param name="list">the set of the indexes of selected mating parents</param>
		/// <param name="cid">the id of current subproblem</param>
		/// <param name="type">1 - neighborhood; otherwise - whole  population</param>
        public int ACOrSelection2(int cid, int type, double[] pro)
        {
            double r;
            int b = neighborhood[cid][0];
            Solution[] parents = new Solution[t];
            Solution[] parent = new Solution[populationSize];
            double[] fitness = new double[t];
            double[] fit = new double[populationSize];
            Dictionary<int, double> tmp = new Dictionary<int, double>();
            Dictionary<int, double> resultSort = new Dictionary<int, double>();
            double a1 = 0;
            double a2 = 0;

            r = JMetalRandom.NextDouble();
            for (int k = 0; k < pro.Length; k++)
            {
                a2 = a2 + pro[k];
                if (r < a2 && r >= a1)
                {
                    b = k;
                    break;
                }
                a1 = a1 + pro[k];
            }

            if (type == 1)
            {
                for (int i = 0; i < t; i++)
                {
                    parents[i] = population.Get(neighborhood[cid][i]);
                    fitness[i] = FitnessFunction(parents[i], lambda[cid]);
                    tmp.Add(neighborhood[cid][i], fitness[i]);
                }
            }
            else
            {
                for (int i = 0; i < populationSize; i++)
                {
                    parent[i] = population.Get(i);
                    fit[i] = FitnessFunction(parent[i], lambda[cid]);
                    tmp.Add(i, fit[i]);
                }
            }
            resultSort = tmp.OrderBy(Data => Data.Value).ToDictionary(keyvalue => keyvalue.Key, keyvalue => keyvalue.Value);
            int[] person = resultSort.Keys.ToArray();

            return person[b];
        }

        /// <summary>
		/// 
		/// </summary>
		/// <param name="cid">the id of current subproblem</param>
		/// <param name="type">1 - neighborhood; otherwise - whole  population</param>
        public int ACOrSelection3(int cid, int type)
        {
            int ss;
            ss = neighborhood[cid].Length;
            double r;
            int p = neighborhood[cid][0];
            Solution[] parents = new Solution[ss];
            Solution[] parent = new Solution[populationSize];
            double[] fitness = new double[ss];
            double[] fit = new double[populationSize];
            int indexOfmin = 0;
            double sum = 0;
            double[] pro = new double[populationSize];
            double a1 = 0;
            double a2 = 0;

            if (type == 1)
            {
                indexOfmin = JMetalRandom.Next(0, ss - 1);
                p = neighborhood[cid][indexOfmin];
            }
            else
            {
                for (int i = 0; i < populationSize; i++)
                {
                    parent[i] = population.Get(i);
                    fit[i] = 1 / FitnessFunction(parent[i], lambda[cid]);
                    sum = sum + fit[i];
                }
                for (int j = 0; j < populationSize; j++)
                {
                    pro[j] = fit[j] / sum;
                }
                r = JMetalRandom.NextDouble();
                for (int k = 0; k < pro.Length; k++)
                {
                    a2 = a2 + pro[k];
                    if (r < a2 && r >= a1)
                    {
                        p = k;
                        break;
                    }
                    a1 = a1 + pro[k];
                }
            }

            return p;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="indiv">child solution</param>
        /// <param name="id">the id of current subproblem</param>
        /// <param name="type">update solutions in - neighborhood (1) or whole population (otherwise)</param>
        private void UpdateProblem(Solution indiv, int id, int type)
        {
            int size;
            int time;

            time = 0;

            if (type == 1)
            {
                size = neighborhood[id].Length;
            }
            else
            {
                size = population.Size();
            }
            int[] perm = new int[size];

            JMetalCSharp.Metaheuristics.MOEAD.Utils.RandomPermutation(perm, size);

            for (int i = 0; i < size; i++)
            {
                int k;
                if (type == 1)
                {
                    k = neighborhood[id][perm[i]];
                }
                else
                {
                    k = perm[i];      // calculate the values of objective function regarding the current subproblem
                }
                double f1, f2;

                f1 = FitnessFunction(population.Get(k), lambda[k]);
                f2 = FitnessFunction(indiv, lambda[k]);

                if (f2 < f1)
                {
                    population.Replace(k, new Solution(indiv));
                    time++;
                }
                // the maximal number of solutions updated is not allowed to exceed 'limit'
                if (time >= nr)
                {
                    return;
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="indiv">child solution</param>
        /// <param name="id">the id of current subproblem</param>
        /// <param name="type">update solutions in - neighborhood (1) or whole population (otherwise)</param>
        private int UpdateProblemWithReplace(Solution indiv, int id, int type)
        {
            int size;
            int time;

            time = 0;

            if (type == 1)
            {
                size = neighborhood[id].Length;
            }
            else
            {
                size = population.Size();
            }
            int[] perm = new int[size];

            JMetalCSharp.Metaheuristics.MOEAD.Utils.RandomPermutation(perm, size);

            for (int i = 0; i < size; i++)
            {
                int k;
                if (type == 1)
                {
                    k = neighborhood[id][perm[i]];
                }
                else
                {
                    k = perm[i];      // calculate the values of objective function regarding the current subproblem
                }
                double f1, f2;

                f1 = FitnessFunction(population.Get(k), lambda[k]);
                f2 = FitnessFunction(indiv, lambda[k]);

                if (f2 < f1)
                {
                    population.Replace(k, new Solution(indiv));
                    time++;
                }
                // the maximal number of solutions updated is not allowed to exceed 'limit'
                if (time >= nr)
                {
                    break;
                }
            }
            return time;
        }

        private double FitnessFunction(Solution individual, double[] lambda)
        {
            double fitness;
            fitness = 0.0;

            //if (functionType == "_TCHE1")
            //{
                double maxFun = -1.0e+30;

                for (int n = 0; n < problem.NumberOfObjectives; n++)
                {
                    //double diff = Math.Abs((individual.Objective[n] - z[n]) / (znad[n] - z[n]));
                    double diff = Math.Abs(individual.Objective[n] - z[n]);

                    double feval;
                    if (lambda[n] == 0)
                    {
                        feval = 0.0001 * diff;
                    }
                    else
                    {
                        feval = diff * lambda[n];
                    }
                    if (feval > maxFun)
                    {
                        maxFun = feval;
                    }
                }

                fitness = maxFun;
            /*}
            else
            {
                Logger.Log.Error("MOEAD.FitnessFunction: unknown type " + functionType);
                Console.WriteLine("MOEAD.FitnessFunction: unknown type " + functionType);
                Environment.Exit(-1);
                throw new Exception("MOEAD.FitnessFunction: unknown type " + functionType);

            }*/
            return fitness;
        }

        private void GetStdDev(int[][] n)
        {
            //GetFS();
            for (int i = 0; i < populationSize; i++)
            {
                XReal xTmp = new XReal(population.Get(i));
                for (int k = 0; k < population.Get(i).NumberOfVariables(); k++)
                {
                    double r = 0;
                    for (int l = 0; l < n[i].Length; l++)
                    {
                        XReal xTmp1 = new XReal(population.Get(n[i][l]));
                        double abs = Math.Abs(xTmp1.GetValue(k) - xTmp.GetValue(k));
                        //if (abs < Math.Abs(xTmp.GetUpperBound(k) - xTmp.GetLowerBound(k)) / n[i].Length)
                        r = r + (abs / (n[i].Length - 1));
                    }
                    xTmp.SetstdDev(k, r);
                }
            }
        }

        private void GetStdDev1(int[][] n, int type)
        {
            if (type == 1)
            {
                for (int i = 0; i < populationSize; i++)
                {
                    XReal xTmp = new XReal(population.Get(i));
                    for (int k = 0; k < population.Get(i).NumberOfVariables(); k++)
                    {
                        double r = 0;
                        /*for (int l = 0; l < n[i].Length; l++)
                        {
                            XReal xTmp1 = new XReal(population.Get(n[i][l]));
                            double abs = Math.Abs(xTmp1.GetValue(k) - xTmp.GetValue(k));
                            r = r + (abs / (n[i].Length - 1));
                        }*/
                        for (int l = 0; l < populationSize; l++)
                        {
                            XReal xTmp1 = new XReal(population.Get(l));
                            double abs = Math.Abs(xTmp1.GetValue(k) - xTmp.GetValue(k));
                            r = r + (abs / (populationSize - 1));
                        }
                        xTmp.SetstdDev(k, r);
                    }
                }
            }
            else
            {
                for (int i = 0; i < populationSize; i++)
                {
                    XReal xTmp = new XReal(population.Get(i));
                    for (int k = 0; k < population.Get(i).NumberOfVariables(); k++)
                    {
                        double r = 0;
                        for (int l = 0; l < n[i].Length; l++)
                        {
                            XReal xTmp1 = new XReal(population.Get(n[i][l]));
                            double abs = Math.Abs(xTmp1.GetValue(k) - xTmp.GetValue(k));
                            r = r + (abs / (n[i].Length - 1));
                        }
                        /*for (int l = 0; l < populationSize; l++)
                        {
                            XReal xTmp1 = new XReal(population.Get(l));
                            double abs = Math.Abs(xTmp1.GetValue(k) - xTmp.GetValue(k));
                            r = r + (abs / (populationSize - 1));
                        }*/
                        xTmp.SetstdDev(k, r);
                    }
                }
            }
        }

        public int GetEvaluations()
        {
            return evaluations;
        }

        #endregion
    }
}
