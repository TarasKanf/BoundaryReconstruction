using CrackProblem.Contracts;
using System;

namespace CrackProblem.Helpers
{
    /// <summary>
    /// Continious parametrized curve with derivetive
    /// </summary>
    public class ApproxParametrizedCurve : IParametrizedCurve
    {
        private readonly TrigonPolinom _xApproxPolinom;
        private readonly TrigonPolinom _yApproxPolinom;

        public ApproxParametrizedCurve(TrigonPolinom xApproxPolinom, TrigonPolinom yApproxPolinom)
        {
            _xApproxPolinom = xApproxPolinom;
            _yApproxPolinom = yApproxPolinom;
        }

        public double GetX(double t)
        {
            return _xApproxPolinom.Value(t);
        }

        public double GetY(double t)
        {
            return _yApproxPolinom.Value(t);
        }

        public double GetDerivetiveX(double t)
        {
            return _xApproxPolinom.Derivative(t);
        }

        public double GetDerivetiveY(double t)
        {
            return _yApproxPolinom.Derivative(t);
        }
    }
}
