using System;
using CrackProblem.Helpers;
using CrackProblem.Integrals;

namespace CrackProblem
{
    using Contracts;

    public class InverseCrackProblemState : DirectProblemState
    {
        public double[] DerivativeOnOuterCurve { get; set; }
        public double[] Density { get; set; }

        public DoubleCore<double> FreshetCore { get; }
        public Func<double, double> LinDataEquationRightPart { get; } 

        public IInversProblemTestData InversProblemTestData { get; set; }

        public InverseCrackProblemState(double radius, int pointsNumber, IInversProblemTestData inverseProblemTestData) 
            : base (radius, pointsNumber, inverseProblemTestData.GetDirectProblemData())
        {
            InversProblemTestData = inverseProblemTestData;
            FreshetCore = new DoubleCore<double>(FreshetCoreFunction);
            LinDataEquationRightPart = new Func<double, double>(LinDataEquationRightPartFucntion);
        }

        /// <summary>
        /// Core function of linealized data equation
        /// </summary>
        /// <param name="t"></param>
        /// <param name="tau"></param>
        /// <returns></returns>
        private double FreshetCoreFunction(double t, double tau)
        {
            return 0;
        }

        private double FreshetByXCoreFunction(double t, double tau)
        {
            double x = InnerCurve.GetX(tau), y = InnerCurve.GetY(tau);

            double denominator = x * x + y * y
                + Radius * Radius - 2.0 * x * Radius * Math.Cos(t)
                - 2.0 * y * Radius * Math.Sin(t);
            double firstTerm = -(x * x + y * y - Radius * Radius)
                               * (2.0 * x - 2.0 * Radius * Math.Cos(t))
                               / (2.0 * Math.PI * Radius * Math.Pow(denominator, 2));
            double secondTerm = x / (Math.PI * Radius * denominator);
            double result = firstTerm + secondTerm;
            return result;
        }

        private double FreshetByYCoreFunction(double t, double tau)
        {
            double x = InnerCurve.GetX(tau), y = InnerCurve.GetY(tau);

            double denominator = x * x + y * y
                + Radius * Radius - 2.0 * x * Radius * Math.Cos(t)
                - 2.0 * y * Radius * Math.Sin(t);
            double firstTerm = -(x * x + y * y - Radius * Radius)
                               * (2.0 * y - 2.0 * Radius * Math.Sin(t))
                               / (2.0 * Math.PI * Radius * Math.Pow(denominator, 2));
            double secondTerm = y / (Math.PI * Radius * denominator);
            double result = firstTerm + secondTerm;
            return result;
        }

        private double LinDataEquationRightPartFucntion(double t)
        {
            var core = new DoubleCore<Point>(DataEquationOperatorCore);
            core.Prepare(new Point(
                OuterCurve.GetX(t),
                OuterCurve.GetY(t)));

           var result = - Integral.CalculateWithTrapeziumMethod(Density, core) / 2.0
            - Omega1(t);
            return result;
        }
    }
}
