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

        public void FormDiscreteEquation(out double[,] matrix, out double[] rightPart)
        {
            double H = 2.0 * Math.PI / PointsNumber;

            // calculate right part vector
            rightPart = new double[this.PointsNumber];
            double ti = 0;
            for (int i = 0; i < PointsNumber; i++)
            {
                rightPart[i] = RightPartFunction(ti);
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
                    matrix[i, j] = Core.GetValue(sj);
                    sj += H;
                }
                ti += H;
            }

            Printer.WriteLine("Right part :");
            Printer.Write(rightPart);
            Printer.WriteLine("Matrix :");
            Printer.Write(matrix);
        }
    }    
}
