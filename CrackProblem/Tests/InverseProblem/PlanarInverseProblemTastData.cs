using CrackProblem.Contracts;
using CrackProblem.Helpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Tests
{
    public class PlanarInverseProblemTastData : IInversProblemTestData
    {
        public ParametrizedCurve CorrectInnerCurve { get; set; } = new ParametrizedCurve(InnerXFuntion, InnerYFunction,
                InnerXFuntionDerivetive, InnerYFunctionDerivetive);
        public Func<Point, double> OnCrackCurveValue { get; set; } = (Point p) => 0;
        public Func<double, double> OnOuterCurveValue { get; set; } = (double t) => 1 + Math.Cos(t)*Math.Cos(t);
        public IDirectProblemTestData GetDirectProblemData()
        {
            return new PlanarTestData()
            {
                OnCrackCurveValue = this.OnCrackCurveValue,
                OnOuterCurveValue = this.OnOuterCurveValue,
                OnTestCurveValue = null
            };
        }

        private static double InnerXFuntion(double t)
        {
            //return Math.Cos(t); // TODO
            return 0.5 * Math.Cos(t);
        }

        private static double InnerXFuntionDerivetive(double t)
        {
            //return - Math.Sin(t);
            return 0.5;
        }

        private static double InnerYFunction(double t)
        {
            //return Math.Cos(t); // TODO
            return 0.5 * (Math.Cos(t) * Math.Cos(t) - 0.5);
        }

        private static double InnerYFunctionDerivetive(double t)
        {
            //return 1;
            return Math.Cos(t);
        }
    }
}
