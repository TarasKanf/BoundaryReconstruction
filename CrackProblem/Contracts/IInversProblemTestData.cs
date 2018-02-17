using CrackProblem.Helpers;
using System;

namespace CrackProblem.Contracts
{
    public interface IInversProblemTestData
    {
        ParametrizedCurve CorrectInnerCurve { get; set; }
        Func<Point, double> OnCrackCurveValue { get; set; }
        Func<Point, double> OnOuterCurveValue { get; set; }

        IDirectProblemTestData GetDirectProblemData();
    }
}
