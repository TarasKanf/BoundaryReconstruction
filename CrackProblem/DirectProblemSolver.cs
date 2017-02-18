using CrackProblem.Contracts;
using CrackProblem.Integrals;
using System;

namespace CrackProblem
{
    using Helpers;
    using LinearSystem;
    using Tests;

    /// <summary>
    /// Forms matrix and right part for direct problem
    /// </summary>
    public class DirectProblemSolver
    {
        private CrackProblemState _state;

        public DirectProblemSolver(double radius, int pointsNumber)
        {
            _state = new CrackProblemState(radius, pointsNumber, new PlanarTestData());
        }

        public DirectProblemSolver(CrackProblemState state)
        {
            _state = state;
        }

        public double[] CalculateDensity()
        {
            IntegralEquationDiscretezer discretezer = new IntegralEquationDiscretezer(
                _state.PointsNumber,
                _state.DirectProblemCore,
                _state.RightPartFunction);

            double[,] matrix;
            double[] rightPart;
            discretezer.FormDiscreteEquation(out matrix, out rightPart);

            SystemOfLinearEquations systemSolver = new SystemOfLinearEquations();
            double[] density = systemSolver.LU_methodSolving(matrix, rightPart, rightPart.Length);

            return density;
        }

        public double[] BuildSolutionOn(double[] density, 
            Func<double,double> xFunc, 
            Func<double,double> yFunc, 
            double paramStart, 
            double paramEnd,
            int solutionPointsNumber)
        {
            double t = 0;
            double h = (paramEnd - paramStart) / solutionPointsNumber;
            double[] result = new double[solutionPointsNumber];

            for (int i = 0; i < solutionPointsNumber; i++)
            {
                Point point = new Point(xFunc(t), yFunc(t));
                result[i] = _state.BuildSolution(density, point);
                t += h;
            }

            return result;
        }
    }
}
