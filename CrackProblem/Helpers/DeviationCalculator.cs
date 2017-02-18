using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Helpers
{
    public class DeviationCalculator
    {
        public double MaxDeviation(double[] aproximateValues, 
            Func<double, double> accurate,
            double paramStart, 
            double paramEnd)
        {
            double step = (paramEnd - paramStart) / aproximateValues.Length;
            if(aproximateValues.Length == 0)
            {
                throw new Exception("Array of values could not be empty");
            }

            double maxDeviation = Math.Abs(aproximateValues[0] - accurate(paramStart));

            int i = 0;
            List<double> accValues = new List<double>();
            for(double t = paramStart; t < paramEnd && i< aproximateValues.Length; t+= step, i++)
            {
                double deviation = Math.Abs(aproximateValues[i] - accurate(t));
                accValues.Add(accurate(t));

                maxDeviation = maxDeviation < deviation ? deviation : maxDeviation;
            }

            Printer.WriteLine("Accurate solution :");
            Printer.Write(accValues.ToArray());

            return maxDeviation;
        }
    }
}
