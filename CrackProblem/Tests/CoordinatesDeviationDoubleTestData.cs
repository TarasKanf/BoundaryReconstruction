using CrackProblem.Contracts;
using CrackProblem.Helpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Tests
{
    public class CoordinatesDeviationDoubleTestData : ITestData
    {
        public Func<Point, double> OnCrackCurveValue { get; set; } = (Point p) => Solution(p);
        public Func<Point, double> OnOuterCurveValue { get; set; } = (Point p) => Solution(p);
        public Func<Point, double> OnTestCurveValue { get; set; } = (Point p) => Solution(p);

        private static double Solution(Point x)
        {
            double result = x.X * x.X - x.Y * x.Y;

            return result;
        }
    }
}
