using CrackProblem.Contracts;
using System;

namespace CrackProblem.Helpers
{
    public class StarCurve : IParametrizedCurve
    {
        private ISource<double> valueSource;

        public StarCurve(ISource<double> radialFunction)
        {
            valueSource = radialFunction;
        }

        public StarCurve(Func<double, double> radialFunction)
        {
            valueSource = new ParametrizedSource(radialFunction);
        }

        public double GetX(double t)
        {
            return valueSource.Value(t) * Math.Cos(t);
        }

        public double GetY(double t)
        {
            return valueSource.Value(t) * Math.Sin(t);
        }

        public double GetDerivetiveX(double t)
        {
            throw new NotImplementedException();
        }

        public double GetDerivetiveY(double t)
        {
            throw new NotImplementedException();
        }
    }
}
