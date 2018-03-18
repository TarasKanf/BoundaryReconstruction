using CrackProblem.Contracts;
using CrackProblem.Integrals;
using System;
using System.Collections.Generic;
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
                Printer.WriteLine("Correct density:");
                Printer.Write(density);
                _state.DerivativeOnOuterCurve = _state.BuildSolutionDerivativeOnOuterCurve(density);
                Printer.WriteLine("Derivetive on outer curve:");
                Printer.Write(_state.DerivativeOnOuterCurve);

                // change initial state back
                _state.InnerCurve = initialInnerCurve;
            }

            double[] points = IntegralEquationDiscretezer.GetDiscretePoints(_state.PointsNumber);
            Printer.WriteLine("Innitial inner curve X:");
            Printer.Write(GetCurveValues(_state.InnerCurve, points, onX: true));
            Printer.WriteLine("Innitial inner curve Y:");
            Printer.Write(GetCurveValues(_state.InnerCurve, points, onX: false));

            bool innerCurveFound = false;
            int iteration = 0;
            while (iteration < 6 )
            {
                DirectProblemSolver solver = new DirectProblemSolver(_state);
                _state.Density = solver.CalculateDensity();
                Printer.WriteLine("Calculated density");
                Printer.Write(_state.Density);

                var curveCorrection = CalculateCurveCorrection();
                var correctedCurve = BuildNewCurve(_state.InnerCurve,
                    new ArraySegment<double>(curveCorrection, 0, _state.CorrectionPolinomPower + 1),
                    new ArraySegment<double>(curveCorrection, _state.CorrectionPolinomPower + 1,
                        _state.CorrectionPolinomPower + 1));
                // TODO calculate exit conditiona
               
                Printer.WriteLine("Correct curve X");
                Printer.Write(GetCurveValues(_state.InversProblemTestData.CorrectInnerCurve, points, onX: true));
                Printer.WriteLine("Solution curve x");
                Printer.Write(GetCurveValues(correctedCurve, points, onX: true));
                Printer.WriteLine("Correct curve y");
                Printer.Write(GetCurveValues(_state.InversProblemTestData.CorrectInnerCurve, points, onX: false));
                Printer.WriteLine("Solution curve y");
                Printer.Write(GetCurveValues(correctedCurve, points, onX: false));

                _state.InnerCurve = correctedCurve;
                iteration++;
            }
        }

        private double[] CalculateCurveCorrection()
        {
            IntegralTwoCoreEquationDiscretezer discretezer = new IntegralTwoCoreEquationDiscretezer(
                _state.PointsNumber,
                _state.FreshetCoreByX,
                _state.FreshetCoreByY,
                _state.LinDataEquationRightPart);

            double[,] matrix;
            double[] rightPart;
            int columnsNumber = _state.CorrectionPolinomPower * 2 + 2;
            discretezer.FormDiscreteEquation(out matrix, out rightPart, columnsNumber,
                (value, i) => _state.DerivativeOnOuterCurve[i] + value);

            var tihanovRegularization = new TihanovRegularization(matrix, rightPart, rightPart.Length, columnsNumber);
            var lambda = 0.1;
            var curveCorrection = tihanovRegularization.Solve(lambda);
            return curveCorrection;
        }

        private IParametrizedCurve BuildNewCurve(IParametrizedCurve currentCurve, IList<double> xCurveCorrection, IList<double> yCurveCorrection)
        {
            var xChebishevpolinom = new ChebishevPolinom(xCurveCorrection, _state.CorrectionPolinomPower);
            var yChebishevpolinom = new ChebishevPolinom(yCurveCorrection, _state.CorrectionPolinomPower);
            double[] points = IntegralEquationDiscretezer.GetDiscretePoints(_state.PointsNumber);
            double[] newXValues = new double[points.Length];
            double[] newYValues = new double[points.Length];
            for (int i = 0; i < points.Length; i++)
            {
                newXValues[i] = _state.InnerCurve.GetX(points[i]) + xChebishevpolinom.Value(points[i]);
                newYValues[i] = _state.InnerCurve.GetY(points[i]) + yChebishevpolinom.Value(points[i]);
            }
            return new ApproxParametrizedCurve(
                new TrigonPolinom(newXValues, newXValues.Length / 2), 
                new TrigonPolinom(newYValues, newYValues.Length/2));
        }

        private double[] GetCurveValues(IParametrizedCurve curve, double[] points, bool onX)
        {
            double[] values = new double[_state.PointsNumber];
            for (int i = 0; i < values.Length; i++)
            {
                values[i] = onX ? curve.GetX(points[i]) : curve.GetY(points[i]);
            }
            return values;
        }
    }
}
