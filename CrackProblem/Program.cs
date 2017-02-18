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
            double radius = 2;
            int pointsNumber = 16;

            ITestData testData = new PlanarTestData();
            //ITestData testData = new FundamentalSolutionTestData();
            CrackProblemState state = new CrackProblemState(radius, pointsNumber, testData);
            DirectProblemSolver solver = new DirectProblemSolver(state);

            Printer.Append = false;
            Printer.Mode = WriteMode.Console;
            double[] density = solver.CalculateDensity();
            Printer.WriteLine("Density :");
            Printer.Write(density);
            double[] solution = solver.BuildSolutionOn(density,
                xFunc: SolutionCurveX,
                yFunc: (t) => t / (2.0 * Math.PI) - 0.5,
                paramStart: 0,
                paramEnd: 2 * Math.PI,
                solutionPointsNumber: 6);

            DeviationCalculator calculator = new DeviationCalculator();
            double deviation = calculator.MaxDeviation(
                solution,
                accurate: (double t) =>
                {
                    Point x = new Point(SolutionCurveX(t), SolutionCurveY(t));
                    return testData.OnTestCurveValue(x);
                },
                paramStart: 0,
                paramEnd: 2 * Math.PI);

            Printer.Append = true;
            Printer.WriteLine("Solution :");
            Printer.Write(solution);
            Printer.WriteLine($"Deviation : {deviation}");
            Console.ReadKey();
        }

        public static double SolutionCurveX(double param)
        {
            return param / (2.0 * Math.PI) - 1.5;
        }

        public static double SolutionCurveY(double param)
        {
            return param / (2.0 * Math.PI) - 0.5;
        }
    }
}
