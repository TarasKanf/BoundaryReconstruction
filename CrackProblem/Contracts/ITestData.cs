using CrackProblem.Helpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Contracts
{
    public interface ITestData
    {
        Func<Point, double> OnCrackCurveValue { get; set; }
        Func<Point,double> OnOuterCurveValue { get; set; }
        Func<Point, double> OnTestCurveValue { get; set; }
    }
}
