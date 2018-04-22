using CrackProblem.Contracts;
using CrackProblem.Helpers;
using System;

namespace CrackProblem.Integrals
{
    public class IntegralEquationDiscretezer
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
        private readonly DoubleCore<double> Core;
        
        public IntegralEquationDiscretezer(int pointsNumber, DoubleCore<double> core, Func<double,double> rightPartFunction)
        {
            PointsNumber = pointsNumber;
            RightPartFunction = rightPartFunction;            
            Core = core;
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

        public void FormDiscreteEquation(out double[,] matrix, out double[] rightPart, Func<double,int, double> matrixModifier = null, Func<double, int, double> rightPartModifier = null )
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
            matrix = new double[PointsNumber, PointsNumber];
            double sj = 0;
            ti = 0;
            for (int i = 0; i < PointsNumber; i++)
            {
                
                Core.Prepare(ti);             
                sj = 0;
                for (int j = 0; j < PointsNumber; j++)
                {
                    matrix[i, j] = matrixModifier?.Invoke(Core.GetValue(sj), j) ?? Core.GetValue(sj);
                    sj += H;
                }
                ti += H;
            }
        }
    }    
}
