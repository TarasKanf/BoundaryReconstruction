using CrackProblem.Contracts;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CrackProblem.Helpers;

namespace CrackProblem.Tests
{
    public class PlanarTestData : ITestData
    {
        public Func<Point, double> OnCrackCurveValue { get; set; } = (Point p) => 1; 

        public Func<Point, double> OnOuterCurveValue { get; set; } = (Point p) => 1;

        public Func<Point, double> OnTestCurveValue { get; set; } = (Point p) => 1;
    }
}
