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
        }

        public static void SolveInverseProblem()
        {
            double radius = 2;
            int pointsNumber = 16;
            int chebishevpolinomPower = 5;

            IInversProblemTestData testData = new PlanarInverseProblemTastData();
            InverseCrackProblemState state = new InverseCrackProblemState(radius, pointsNumber, chebishevpolinomPower, testData);
            state.InnerCurve = new ParametrizedCurve(InnerXFuntion,InnerYFunction, InnerXFuntionDerivetive, InnerYFunctionDerivetive);
            InverseProblemSolver solver = new InverseProblemSolver(state);
            solver.CalculateCurve();
        }

        private static double InnerXFuntion(double t)
        {
            return Math.Cos(t) + 0.2; // TODO
        }

        private static double InnerXFuntionDerivetive(double t)
        {
            return 1;
            //return - Math.Sin(t);
        }

        private static double InnerYFunction(double t)
        {
            return Math.Cos(t) + 0.2; // TODO
        }

        private static double InnerYFunctionDerivetive(double t)
        {
            return 1;
            //return - Math.Sin(t);
        }

        public static void SolveDirectProblem()
        {
            double radius = 2;
            int pointsNumber = 32;

            //IDirectProblemTestData testData = new PlanarTestData();
            //IDirectProblemTestData testData = new FundamentalSolutionTestData();
            IDirectProblemTestData testData = new FundamentalSolutionDevidedTastData();
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
