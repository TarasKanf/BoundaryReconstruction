using CrackProblem.Contracts;
using System;

namespace CrackProblem.Helpers
{
    public class ParametrizedSource : ISource
    {
        private Func<double, double> _function;

        public ParametrizedSource(Func<double, double> function)
        {
            _function = function;
        }

        public double Value(double x)
        {
            return _function(x);
        }
    }
}
