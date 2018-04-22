using CrackProblem.Helpers;
using System;

namespace CrackProblem.Contracts
{
    public interface IDirectProblemTestData
    {
        Func<Point, double> OnCrackCurveValue { get; set; }
        Func<double,double> OnOuterCurveValue { get; set; }
        Func<Point, double> OnTestCurveValue { get; set; }
    }
}
