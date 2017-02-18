using CrackProblem.Contracts;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CrackProblem.Helpers;

namespace CrackProblem.Tests
{
    public class FundamentalSolutionTestData : ITestData
    {
        public Func<Point, double> OnCrackCurveValue { get; set; } = (Point p) => FundamentalSolution(p);
        public Func<Point, double> OnOuterCurveValue { get; set; } = (Point p) => FundamentalSolution(p);
        public Func<Point, double> OnTestCurveValue { get; set; } = (Point p) => FundamentalSolution(p);

        private static double FundamentalSolution(Point x)
        {
            Point outerPoint = new Point(10, 0);

            double deviationAbs = Math.Sqrt(
                Math.Pow(x.X - outerPoint.X, 2)
                +
                Math.Pow(x.Y - outerPoint.Y, 2));

            double result = Math.Log(deviationAbs);

            return result;
        }
    }
}
