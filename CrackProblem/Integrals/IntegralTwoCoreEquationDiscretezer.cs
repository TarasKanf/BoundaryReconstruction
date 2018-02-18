using CrackProblem.Contracts;
using CrackProblem.Helpers;
using System;

namespace CrackProblem.Integrals
{
    public class IntegralTwoCoreEquationDiscretezer
    {
        private int pointsNumber = 0;
        public int PointsNumber
        {
            get
            {
                return pointsNumber;
            }
            private set
            {
                if (value % 2 != 0)
                {
                    throw new ArgumentException("Count of point must be even");
                }
                else
                {
                    pointsNumber = value;
                }
            }
        }

        private readonly Func<double,double> RightPartFunction;
        private readonly DoubleCore<double> CoreFirst;
        private readonly DoubleCore<double> CoreSecond;

        public IntegralTwoCoreEquationDiscretezer(int pointsNumber, DoubleCore<double> core1, DoubleCore<double> core2, Func<double, double> rightPartFunction)
        {
            PointsNumber = pointsNumber;
            RightPartFunction = rightPartFunction;
            CoreFirst = core1;
            CoreSecond = core2;
        }

        public static double CalculateDiscreteStep(int pointsNumber)
        {
            return 2.0 * Math.PI / pointsNumber;
        }

        public static double[] GetDiscretePoints(int pointsNumber)
        {
            double H = CalculateDiscreteStep(pointsNumber);
            double[] points = new double[pointsNumber];
            double ti = 0;
            for (int i = 0; i < pointsNumber; i++)
            {
                points[i] = ti;
                ti += H;
            }
            return points;
        }

        public void FormDiscreteEquation(out double[,] matrix, out double[] rightPart, int columnsNumber, Func<double, int, double> rightPartModifier)
        {
            double H = CalculateDiscreteStep(PointsNumber);

            // calculate right part vector
            rightPart = new double[this.PointsNumber];
            double ti = 0;
            for (int i = 0; i < PointsNumber; i++)
            {
                rightPart[i] = rightPartModifier?.Invoke(RightPartFunction(ti), i) ?? RightPartFunction(ti);
                ti += H;
            }

            // calculate matrix
            int columnsPerCore = columnsNumber  / 2; // chebichev polinom power
            matrix = new double[PointsNumber,columnsNumber];
            double sj = 0;
            ti = 0;

            for (int i = 0; i < PointsNumber; i++)
            {
                CoreFirst.Prepare(ti);
                CoreSecond.Prepare(ti);
                for (int k = 0; k < columnsPerCore; k++)
                {
                    matrix[i, k] = CoreFirst.GetValue(k);
                    matrix[i, k + columnsPerCore] = CoreSecond.GetValue(k);
                }
                ti += H;
            }
        }
    }    
}
