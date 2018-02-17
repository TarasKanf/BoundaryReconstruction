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
        public Func<Point, double> OnCrackCurveValue { get; set; } = (Point p) => 1;
        public Func<Point, double> OnOuterCurveValue { get; set; } = (Point p) => 1;
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
            return Math.Cos(t); // TODO
        }

        private static double InnerXFuntionDerivetive(double t)
        {
            return 1;
            //return - Math.Sin(t);
        }

        private static double InnerYFunction(double t)
        {
            return Math.Cos(t); // TODO
        }

        private static double InnerYFunctionDerivetive(double t)
        {
            return 1;
            //return - Math.Sin(t);
        }
    }
}
