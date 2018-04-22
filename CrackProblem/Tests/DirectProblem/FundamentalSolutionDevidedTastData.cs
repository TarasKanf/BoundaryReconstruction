using CrackProblem.Contracts;
using CrackProblem.Helpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Tests
{
    //public class FundamentalSolutionDevidedTastData : IDirectProblemTestData
    //{
    //    public Func<Point, double> OnCrackCurveValue { get; set; } = (Point p) => FundamentalSolution(p);
    //    public Func<Point, double> OnOuterCurveValue { get; set; } = (Point p) => FundamentalSolution(p);
    //    public Func<Point, double> OnTestCurveValue { get; set; } = (Point p) => FundamentalSolution(p);

    //    private static double FundamentalSolution(Point x)
    //    {
    //        Point outerPoint = new Point(10, 0);

    //        double deviationAbs = Math.Sqrt(
    //            Math.Pow(x.X - outerPoint.X, 2)
    //            +
    //            Math.Pow(x.Y - outerPoint.Y, 2));

    //        double result = Math.Log(1.0/deviationAbs) / (2.0 * Math.PI);

    //        return result;
    //    }
    //}
}
