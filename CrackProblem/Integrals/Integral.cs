using CrackProblem.Contracts;
using System;

namespace CrackProblem.Integrals
{
    public static class Integral
    {
        // for periodik function on [0, 2*PI]
        public static double CalculateWithTrapeziumMethod(Func<double, double> f, int n)
        {
            n = DevideByTwo(n);
            int N = 2 * n;
            double temp = 0;
            double h = Math.PI / n;
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += f(temp);
                temp += h;
            }
            return sum * Math.PI / n;
        }

        public static double CalculateWithTrapeziumMethod(double[] density, Func<double, double> f)
        {
            if (density.Length % 2 != 0) throw new ArgumentException("Density must have even length");
            int N = density.Length;
            int n = density.Length / 2;
            double temp = 0;
            double h = Math.PI / n;
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += density[i] * f(temp);
                temp += h;
            }
            return sum * Math.PI / n;
        }

        // for periodik function f on [0, 2*PI]
        public static double CalculateWithTrapeziumMethod(ICore f,int n) 
        {
            n = DevideByTwo(n);
            int N = 2*n;
            double temp = 0;
            double h = Math.PI/n;
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += f.GetValue(temp);
                temp += h;
            }
            return sum*Math.PI / n;
        }

        // for periodic hyper singular function f on [0, 2*PI] and f = g*/(2*(sin[(t-k)/2])^2)
        public static double CalculateWithHyperSingularCore(DoubleCore<double> f, int n)
        {
            n = DevideByTwo(n);
            int N = 2 * n;
            double temp = 0;
            double h = Math.PI / n;
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += f.GetValue(temp)*CoefficientForHyperSingular(f.Param,n,temp);
                temp += h;
            }
            return sum *2.0* Math.PI ;
        }

        // for periodiс cores with smooth part f on [0, 2*PI]x[0, 2*PI]
        public static double CalculateWithWeakSingularCore(DoubleCore<double> f, int n)
        {
            n = DevideByTwo(n);
            int N = 2 * n;
            double temp = 0;
            double h = Math.PI / n;
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += f.GetValue(temp) * CoefficientForWeakSingular(f.Param, n, temp);
                temp += h;
            }
            return sum * 2.0 * Math.PI;
        }

        public static double CoefficientForWeakSingular(double t, int n, double ti)
        {
            n = DevideByTwo(n);
            //якщо ln((4/e)*(...))
            double sum = 0;
            double delta_t = (Math.Abs(t - ti) < 1e-10) ? 0 : t - ti;
            for (int i = 1; i < n; i++)
            {
                sum += Math.Cos(i * (delta_t)) / i;
            }
            sum *= 2.0;
            sum += 1.0;
            sum += Math.Cos(n * (delta_t)) / n;
            sum *= -1.0 / (2.0 * n);
            return sum;

            // якщо ln((4)*(...))
            //double sum = 0;
            //for (int i = 1; i < n; i++)
            //{
            //    sum += Math.Cos(i * (t - ti)) / i;
            //}
            //sum *= -1.0 / n;
            //sum -= Math.Cos(n * (t - ti)) / (2.0 * n * n);
            //return sum;           
        }
        private static double CoefficientForHyperSingular(double t, int n, double ti)
        {
            n = DevideByTwo(n);
            double sum = 0;
            double delta_t = (Math.Abs(t - ti) < 1e-10) ? 0 : t - ti;
            for (int i = 1; i < n; i++)
            {
                sum += i * Math.Cos(i * (delta_t));
            }
            sum *= -1.0 / (double)n;
            sum -= Math.Cos(n * (delta_t)) / 2.0;
            return sum;
        }

        private static int DevideByTwo(int n)
        {
            if (n % 2 != 0) throw new Exception("Two calculate integral points number must be even number");
            return n / 2;
        }
    }
}
