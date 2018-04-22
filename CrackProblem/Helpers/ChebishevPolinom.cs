using System;
using System.Collections.Generic;

namespace CrackProblem.Helpers
{
    public class ChebishevPolinom
    {
        public int N { get; private set; }
        private IList<double> Coefficients { get; set; }

        private const double Eps = 1E-3;

        public ChebishevPolinom(IList<double> coefficients)
        {
            Coefficients = coefficients;
        }
        public double Value(double t)
        {
            double sum = 0;
            for (int i = 0; i < Coefficients.Count; i++)
            {
                sum += Coefficients[i] * Math.Cos(i * t);
            }
            return sum;
        }
        public bool Add(IList<double> coefficients)
        {
            for (int i = 0; i < Coefficients.Count; i++)
            {
                Coefficients[i] += coefficients[i];
            }
            return true;
        }

        public double Derivative(double t)
        {
            double sum = 0;
            //for (int i = 0; i < Coefficients.Count; i++)
            //{
            //    sum += -Coefficients[i] * i * Math.Sin(i * t);
            //}
            var isSinZeroPoint = Math.Abs(t % Math.PI) < Eps;
            for (int i = 0; i < Coefficients.Count; i++)
            {
                sum += Coefficients[i] * (isSinZeroPoint 
                    ? i*i 
                    : i * Math.Sin(i * t)/ Math.Abs(Math.Sin(t)));
            }
            return sum;
        }
    }
}
