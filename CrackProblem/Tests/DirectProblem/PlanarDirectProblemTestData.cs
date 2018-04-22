using CrackProblem.Contracts;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CrackProblem.Helpers;

namespace CrackProblem.Tests
{
    public class PlanarTestData : IDirectProblemTestData
    {
        public Func<Point, double> OnCrackCurveValue { get; set; }

        public Func<double, double> OnOuterCurveValue { get; set; } 

        public Func<Point, double> OnTestCurveValue { get; set; }

        public PlanarTestData()
        {
            OnCrackCurveValue =  (Point p) => 1;
            OnOuterCurveValue = (double p) => 1;
            OnTestCurveValue = (Point p) => 1;
        }
    }
}
