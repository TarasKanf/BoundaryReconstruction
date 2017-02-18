using CrackProblem.Contracts;
using System;

namespace CrackProblem.Helpers
{
    /// <summary>
    /// Continious parametrized curve with derivetive
    /// </summary>
    public class ParametrizedCurve : IParametrizedCurve
    {
        private Func<double, double> _xFunction;
        private Func<double, double> _yFunction;
        private Func<double, double> xDerivetive;
        private Func<double, double> yDerivetive;


        public ParametrizedCurve(
            Func<double, double> xFunction, Func<double, double> yFunction)
        {
            _xFunction = xFunction;
            _yFunction = yFunction;
            // TODO Build polinomials to get derivetives
        }

        public ParametrizedCurve(
            Func<double, double> xFunction,
            Func<double, double> yFunction,
            Func<double, double> xFunctionDerivetive,
            Func<double, double> yFunctionDerivetive)
        {
            _xFunction = xFunction;
            _yFunction = yFunction;
            xDerivetive = xFunctionDerivetive;
            yDerivetive = yFunctionDerivetive;
        }

        public double GetX(double t)
        {
            return _xFunction(t);
        }

        public double GetY(double t)
        {
            return _yFunction(t);
        }

        public double GetDerivetiveX(double t)
        {
            return xDerivetive(t);
        }

        public double GetDerivetiveY(double t)
        {
            return yDerivetive(t);
        }
    }
}
