using CrackProblem.Contracts;
using CrackProblem.Helpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Tests
{
    public class FundamentalSolutionDevidedTastData : IDirectProblemTestData
    {
        private IParametrizedCurve _outerCurve;
        public FundamentalSolutionDevidedTastData(IParametrizedCurve outerCurve)
        {
            _outerCurve = outerCurve;
            OnCrackCurveValue = (Point p) => FundamentalSolution(p);
            OnOuterCurveValue = (p) => FundamentalSolution(p);
            OnTestCurveValue = (Point p) => FundamentalSolution(p);
        }

        public Func<Point, double> OnCrackCurveValue { get; set; }
        public Func<double, double> OnOuterCurveValue { get; set; }
        public Func<Point, double> OnTestCurveValue { get; set; }

        private double FundamentalSolution(Point x)
        {
            Point outerPoint = new Point(10, 0);

            double deviationAbs = Math.Sqrt(
                Math.Pow(x.X - outerPoint.X, 2)
                +
                Math.Pow(x.Y - outerPoint.Y, 2));

            double result = Math.Log(1.0 / deviationAbs) / (2.0 * Math.PI);

            return result;
        }

        private double FundamentalSolution(double t)
        {
            Point x = new Point(_outerCurve.GetX(t), _outerCurve.GetY(t));
            return FundamentalSolution(x);
        }
    }
}
