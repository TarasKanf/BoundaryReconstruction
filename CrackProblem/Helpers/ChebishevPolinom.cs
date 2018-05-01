using System;
using System.Collections.Generic;

namespace CrackProblem.Helpers
{
    #region old polinom
    public class ChebishevPolinomOld
    {
        public int N { get; private set; }
        private IList<double> Coefficients { get; set; }

        private const double Eps = 1E-3;

        public ChebishevPolinomOld(IList<double> coefficients)
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
    #endregion

    public class ChebishevPolinom
    {
        public int N { get; private set; }
        public IList<double> Coefficients { get; set; }
        public IList<double> DerivativeCoefficients { get; set; } 

        private const double Eps = 1E-3;

        public ChebishevPolinom(Func<double, double> funcFromMinusOneToOne, int polinomPower )
        {
            N = polinomPower;
            Coefficients = new double[N+1];
            for (int j = 0; j <= N; j++)
            {
                double sum = 0;
                for (int k = 1; k <= N + 1; k++)
                {
                    sum += funcFromMinusOneToOne(Math.PI*(k - 0.5)/(N + 1))*
                           Math.Cos(Math.PI*j*(k - 0.5)/(N + 1));
                }
                Coefficients[j] = sum*2.0/(N + 1);
            }
            Coefficients[0] = Coefficients[0]/2.0;

            DerivativeCoefficients = new double[N+1];
            DerivativeCoefficients[N] = 0;
            DerivativeCoefficients[N - 1] = 2.0*N*Coefficients[N];
            for (int j = N - 2; j >= 0; j--)
            {
                DerivativeCoefficients[j] = DerivativeCoefficients[j + 2] + 2.0*(j + 1)*Coefficients[j + 1];
            }
            DerivativeCoefficients[0] = DerivativeCoefficients[0] / 2.0;

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

            DerivativeCoefficients[N] = 0;
            DerivativeCoefficients[N - 1] = 2.0 * N * Coefficients[N];
            for (int j = N - 2; j >= 0; j--)
            {
                DerivativeCoefficients[j] = DerivativeCoefficients[j + 2] + 2.0 * (j + 1) * Coefficients[j + 1];
            }
            DerivativeCoefficients[0] = DerivativeCoefficients[0] / 2.0;

            return true;
        }

        public double Derivative(double t)
        {
            //double sum = 0;
            //for (int i = 0; i < Coefficients.Count; i++)
            //{
            //    sum += - Coefficients[i]*i*Math.Sin(i*t);
            //}

            ////for (int i = 0; i < Coefficients.Count; i++)
            ////{
            ////    var x = Math.Cos(t);
            ////    sum += Math.Sin(t) * Coefficients[i] * i * Math.Sin(i * Math.Acos(x))/Math.Sqrt((1 - x*x));
            ////}
            //return sum;

            double sum = 0;
            for (int i = 0; i < DerivativeCoefficients.Count; i++)
            {
                sum += DerivativeCoefficients[i] * Math.Cos(i * t);
            }
            return sum;
        }
    }
}
