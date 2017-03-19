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

            //ITestData testData = new PlanarTestData();
            ITestData testData = new FundamentalSolutionTestData();
            //ITestData testData = new FundamentalSolutionDevidedTastData();
            //ITestData testData = new CoordinatesDeviationDoubleTestData();
            CrackProblemState state = new CrackProblemState(radius, pointsNumber, testData);
            DirectProblemSolver solver = new DirectProblemSolver(state);

            Printer.Append = false;
            Printer.Mode = WriteMode.Console;
            double[] density = solver.CalculateDensity();
            Printer.WriteLine("Density :");
            Printer.Write(density);
            double[] solution = solver.BuildSolutionOn(density,
                xFunc: SolutionCurveX,
                yFunc: SolutionCurveY,
                paramStart: 0,
                paramEnd: 2 * Math.PI,
                solutionPointsNumber: 8);

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
            //return param / (2.0 * Math.PI) - 1.5;
            return param / (2.0 * Math.PI) * 1.5 + 0.05;
           
        }

        public static double SolutionCurveY(double param)
        {
            //return param / (2.0 * Math.PI) - 0.5;
            return 0;
        }
    }
}
