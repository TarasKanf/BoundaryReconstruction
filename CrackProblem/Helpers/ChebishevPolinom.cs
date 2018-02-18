using System;
using System.Collections.Generic;

namespace CrackProblem.Helpers
{
    public class ChebishevPolinom
    {
        public int Power { get; set; }
        private IList<double> index;
        public ChebishevPolinom(IList<double> indexes, int power) // y - array with length 2xn
        {
            Power = power;
            if(indexes.Count < power + 1) throw new Exception("Not anough indexes for specified polinom power.");
            index = indexes;
        }
        public double Value(double x)
        {
            double sum = 0;
            for (int i = 0; i <= Power; i++)
            {
                sum += index[i]*Math.Cos(i*x);
            }
            return sum;
        }
    }
}
