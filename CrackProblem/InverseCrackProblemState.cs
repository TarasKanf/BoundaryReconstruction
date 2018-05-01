using System;
using CrackProblem.Helpers;
using CrackProblem.Integrals;
using CrackProblem.Contracts;

namespace CrackProblem
{
    public class InverseCrackProblemState : DirectProblemState
    {
        public int CorrectionPolinomPower { get; set; }
        public double[] DerivativeOnOuterCurve { get; set; }
        public double[] Density { get; set; }

        public DoubleCore<double> FreshetCoreByX { get; }
        public DoubleCore<double> FreshetCoreByY { get; }
        public Func<double, double> LinDataEquationRightPart { get; } 

        public IInversProblemTestData InversProblemTestData { get; set; }

        public InverseCrackProblemState(double radius, int pointsNumber, int correctionPolinomPower, IInversProblemTestData inverseProblemTestData) 
            : base (radius, pointsNumber, inverseProblemTestData.GetDirectProblemData())
        {
            CorrectionPolinomPower = correctionPolinomPower;
            InversProblemTestData = inverseProblemTestData;
            FreshetCoreByX = new DoubleCore<double>(FreshetByXCoreFunction);
            FreshetCoreByY = new DoubleCore<double>(FreshetByYCoreFunction);
            LinDataEquationRightPart = new Func<double, double>(LinDataEquationRightPartFucntion);
        }

        private double FreshetByXCoreFunction(double t, double k)
        {
            double H = IntegralEquationDiscretezer.CalculateDiscreteStep(PointsNumber);
            double coreSum = 0;
            double sj = 0;
            for (int j = 0; j < PointsNumber; j++)
            {
                double x = InnerCurve.GetX(sj), y = InnerCurve.GetY(sj),
                        dx = InnerCurve.GetDerivetiveX(sj), dy = InnerCurve.GetDerivetiveY(sj);

                double denominator = x * x + y * y
                    + Radius * Radius - 2.0 * x * Radius * Math.Cos(t)
                    - 2.0 * y * Radius * Math.Sin(t);
                double firstTerm = -(x * x + y * y - Radius * Radius)
                    * (2.0 * x - 2.0 * Radius * Math.Cos(t))
                    / (2.0 * Math.PI * Radius * Math.Pow(denominator, 2));
                double secondTerm = x / (Math.PI * Radius * denominator);
                double result = (firstTerm + secondTerm)/2.0;
                double derivetiveNorm = dx*dx + dy*dy;
                double dataEquationCore = DataEquationOperatorCore(new Point(Radius*Math.Cos(t), Radius*Math.Sin(t)),
                    new Point(x, y));
                coreSum += (2.0 * Math.PI / PointsNumber) * Density[j]
                           * (result * Math.Cos(k * sj) + 2.0 * dataEquationCore * dx * ChebishevPolDerivative(k, sj) / derivetiveNorm); // 2.0 * dataEquationCore

                //coreSum += (2.0 * Math.PI / PointsNumber) * Density[j]
                //          * (result * Math.Cos(k * sj));
                sj += H;
            }
            return coreSum;
        }

        private double FreshetByYCoreFunction(double t, double k)
        {
            double H = IntegralEquationDiscretezer.CalculateDiscreteStep(PointsNumber);
            double coreSum = 0;
            double sj = 0;
            for (int j = 0; j < PointsNumber; j++)
            {
                double x = InnerCurve.GetX(sj), y = InnerCurve.GetY(sj),
                          dx = InnerCurve.GetDerivetiveX(sj), dy = InnerCurve.GetDerivetiveY(sj);

                double denominator = x * x + y * y
                    + Radius * Radius - 2.0 * x * Radius * Math.Cos(t)
                    - 2.0 * y * Radius * Math.Sin(t);
                double firstTerm = -(x * x + y * y - Radius * Radius)
                    * (2.0 * y - 2.0 * Radius * Math.Sin(t))
                    / (2.0 * Math.PI * Radius * Math.Pow(denominator, 2));
                double secondTerm = y / (Math.PI * Radius * denominator);
                double result = (firstTerm + secondTerm)/2.0;
                double derivetiveNorm = dx * dx + dy * dy;
                double dataEquationCore = DataEquationOperatorCore(new Point(Radius * Math.Cos(t), Radius * Math.Sin(t)),
                   new Point(x, y));
                coreSum += (2.0 * Math.PI / PointsNumber) * Density[j]
                    * (result * Math.Cos(k * sj) + 2.0 * dataEquationCore * dy * ChebishevPolDerivative(k, sj) / derivetiveNorm);// додано двійку без формули 2.0 * dataEquationCore

                //coreSum += (2.0 * Math.PI / PointsNumber) * Density[j]
                //    * (result * Math.Cos(k * sj));
                sj += H;
            }
            return coreSum;
        }

        private double ChebishevPolDerivative(double power, double t)
        {
            //return -power*Math.Sin(power*t);
            //var result = (t%Math.PI > 0.001)
            //    ? Math.Sin(power*t)*power/(Math.Abs(Math.Sin(t)))
            //    : power*power;
            //var result = -Math.Sin(power*t)*power;
            //var result = Math.Cos(power*t);
            //var result = -power*Math.Sign(Math.Sin(t))*Math.Sin(power*Math.Acos(Math.Cos(t)));
            var x = Math.Cos(t);
            int powerInt = Convert.ToInt32(power);
            double result = 0;
            switch (powerInt)
            {
                case 0: result = 0; break;
                case 1: result = 1; break;
                case 2: result = 4 *x;
                    break;
                case 3: result = 12 *x*x - 3;
                    break;
                case 4: result = 32 *Math.Pow(x, 3) - 16*x;
                    break;
                case 5:result = 80*Math.Pow(x, 4) - 60*x*x + 5;
                    break;
                default:result = double.NaN;
                    break;
            }

            return result;
        }

        private double LinDataEquationRightPartFucntion(double t)
        {
            var core = new DoubleCore<Point>(DataEquationOperatorCore);
            core.Prepare(new Point(
                OuterCurve.GetX(t),
                OuterCurve.GetY(t)));

            var result = - Integral.CalculateWithTrapeziumMethod(Density, core)
            - Omega1(t);
            return result;
        }
    }
}
