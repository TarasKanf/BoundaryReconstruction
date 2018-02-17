using CrackProblem.Contracts;
using CrackProblem.Integrals;
using System;
using CrackProblem.Helpers;
using CrackProblem.LinearSystem;
using CrackProblem.Tests;

namespace CrackProblem
{
    public class InverseProblemSolver
    {
        private readonly InverseCrackProblemState _state;

        public InverseProblemSolver(InverseCrackProblemState state)
        {
            _state = state;
        }

        public void CalculateCurve()
        {
            if (_state.DerivativeOnOuterCurve == null)
            {
                // Solver is in testing mode 
                // Find derivative
                var initialInnerCurve = _state.InnerCurve;
                _state.InnerCurve = _state.InversProblemTestData.CorrectInnerCurve;
                DirectProblemSolver solver = new DirectProblemSolver(_state);
                double[] density = solver.CalculateDensity();
                _state.DerivativeOnOuterCurve = _state.BuildSolutionDerivativeOnOuterCurve(density);

                // change initial state back
                _state.InnerCurve = initialInnerCurve;
            }

            bool innnerCurveFound = false;
            while (!innnerCurveFound)
            {
                DirectProblemSolver solver = new DirectProblemSolver(_state);
                _state.Density = solver.CalculateDensity();

                var curveCorrection = CalculateCurveCorrection();
                // TODO calculate correction to crack
                // TODO calculate exit conditiona
                // TODO calculate new approximation
            }
        }

        private double[] CalculateCurveCorrection()
        {
            IntegralEquationDiscretezer discretezer = new IntegralEquationDiscretezer(
                _state.PointsNumber,
                _state.FreshetCore,
                _state.LinDataEquationRightPart);

            double[,] matrix;
            double[] rightPart;
            discretezer.FormDiscreteEquation(out matrix, out rightPart,
                (value, j) => value * _state.Density[j],
                (value, j) => _state.DerivativeOnOuterCurve[j] + value);

            var tihanovRegularization = new TihanovRegularization(matrix, rightPart, rightPart.Length, rightPart.Length);
            var lambda = 0.01;
            var curveCorrection = tihanovRegularization.Solve(lambda);

            return curveCorrection;
        }


    }
}
