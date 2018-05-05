using CrackProblem.Contracts;
using CrackProblem.Helpers;
using CrackProblem.Tests;
using System;

namespace CrackProblem
{
    class Program
    {
        static void Main(string[] args)
        {
            Printer.Append = false;
            Printer.Mode = WriteMode.File;
            //SolveDirectProblem();
            SolveInverseProblem();
            //IInversProblemTestData testData = new PlanarInverseProblemTastData();
            //Func<double, double> func = (p) => 3.0 * Math.Cos(p)*Math.Cos(p) + 2;
            //Func<double, double> derivative = (p) => 6.0 * Math.Cos(p);

            //Func<double, double> func = (p) => testData.CorrectInnerCurve.GetY(p);
            //Func<double, double> derivative = (p) => testData.CorrectInnerCurve.GetDerivetiveY(p);
            //var pol = new ChebishevPolinom(func, 5);
            //foreach (var coefficient in pol.Coefficients)
            //{
            //    Console.Write($"{coefficient};   ");
            //}
            //double param = 0;
            //int N = 10;
            //double step = 2.0 * Math.PI / N;
            //while (param < 2.0 * Math.PI)
            //{
            //    //var deviation = Math.Abs(func(param) - pol.Value(param));
            //    Console.WriteLine($"Real = {derivative(param)}");
            //    Console.WriteLine($"Approx = {pol.Derivative(param)}");
            //    param += step;
            //}
            //Console.ReadKey();
        }

        public static void SolveInverseProblem()
        {
            double radius = 2;
            int pointsNumber = 16;
            int chebishevpolinomPower = 5;

            IInversProblemTestData testData = new PlanarInverseProblemTastData();
            InverseCrackProblemState state = new InverseCrackProblemState(radius, pointsNumber, chebishevpolinomPower, testData);

            // set initial curve
            state.InnerCurve = new ApproxParametrizedCurve(
                new ChebishevPolinom(InnerXFuntion, chebishevpolinomPower),
                new ChebishevPolinom(InnerYFunction, chebishevpolinomPower));
            InverseProblemSolver solver = new InverseProblemSolver(state);
            solver.CalculateCurve();
        }

        private static double InnerXFuntion(double t)
        {
            //return Math.Cos(t) + 0.2; // TODO
            return 0.5*Math.Cos(t);
        }

        private static double InnerXFuntionDerivetive(double t)
        {
            //return 1;
            //return (Math.Abs(t%Math.PI) > 0.0001)? -0.5*Math.Sin(t): 0.001; // better but wrong
            return 0.5;
        }

        private static double InnerYFunction(double t)
        {
            //return Math.Cos(t) + 0.2; // TODO
            return 0;
        }

        private static double InnerYFunctionDerivetive(double t)
        {
            //return 1;
            return 0;
        }

        public static void SolveDirectProblem()
        {
            double radius = 2;
            int pointsNumber = 32;

            //IDirectProblemTestData testData = new PlanarTestData();
            //IDirectProblemTestData testData = new FundamentalSolutionTestData();
            var outerCurve = new StarCurve((t) => radius);
            IDirectProblemTestData testData = new FundamentalSolutionDevidedTastData(outerCurve);
            //IDirectProblemTestData testData = new CoordinatesDeviationDoubleTestData();
            DirectProblemState state = new DirectProblemState(radius, pointsNumber, testData);
            DirectProblemSolver solver = new DirectProblemSolver(state);

            Printer.Append = false;
            Printer.Mode = WriteMode.Console;
            double[] density = solver.CalculateDensity();
            Printer.WriteLine("Density :");
            Printer.Write(density);
            double[] solution = solver.BuildSolutionOn(density,
                xFunc: SolutionCurveX,
                yFunc: SolutionCurveY,
                paramStart: 0.001,
                paramEnd: Math.PI - 0.001,
                solutionPointsNumber: 16);

            DeviationCalculator calculator = new DeviationCalculator();
            double deviation = calculator.MaxDeviation(
                solution,
                accurate: (double t) =>
                {
                    Point x = new Point(SolutionCurveX(t), SolutionCurveY(t));
                    return testData.OnTestCurveValue(x);
                },
                paramStart: 0.001,
                paramEnd: Math.PI - 0.001);

            Printer.Append = true;
            Printer.WriteLine("Solution :");
            Printer.Write(solution);
            Printer.WriteLine($"Deviation : {deviation}");
            Console.ReadLine();
        }

        public static double SolutionCurveX(double param)
        {
            //return param / (2.0 * Math.PI) - 1.5;
            // return param / (2.0 * Math.PI) * 1.5 + 0.05;
            //return param / (2.0 * Math.PI) * 1.5 + 0.2;
            return Math.Cos(param);
        }

        public static double SolutionCurveY(double param)
        {
            //return param / (2.0 * Math.PI) - 0.5;
            //return 0;
            return Math.Cos(param);
        }
    }
}
