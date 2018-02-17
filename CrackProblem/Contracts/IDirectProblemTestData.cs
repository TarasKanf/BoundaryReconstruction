using CrackProblem.Helpers;
using System;

namespace CrackProblem.Contracts
{
    public interface IDirectProblemTestData
    {
        Func<Point, double> OnCrackCurveValue { get; set; }
        Func<Point,double> OnOuterCurveValue { get; set; }
        Func<Point, double> OnTestCurveValue { get; set; }
    }
}
