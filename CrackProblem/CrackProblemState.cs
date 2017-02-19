namespace CrackProblem
{
    using Contracts;
    using Integrals;
    using CrackProblem.Helpers;
    using System;

    public class CrackProblemState
    {
        private double DeviationEps = 1e-7;
        public double Radius;
        public StarCurve OuterCurve;
        public ParametrizedCurve InnerCurve { get; set; }
        public DoubleCore<double> DirectProblemCore;
        public Func<double,double> RightPartFunction { get; }
        public int PointsNumber { get; set; }
        private ITestData _testData { get; set; }

        public CrackProblemState(double radius, int pointsNumber, ITestData testData)
        {
            Radius = radius;
            OuterCurve = new StarCurve(OuterRadialFunction);
            InnerCurve = new ParametrizedCurve(InnerXFuntion, InnerYFunction, 
                InnerXFuntionDerivetive, InnerYFunctionDerivetive);
            PointsNumber = pointsNumber;
            DirectProblemCore = new DoubleCore<double>(CoreFunction);
            RightPartFunction = DirectProblemRightPartFunction;
            _testData = testData;
        }

        #region Initails
        private double OnEdgeValueFunction(double t)
        {
            Point x = new Point(OuterCurve.GetX(t), OuterCurve.GetY(t));
            return _testData.OnOuterCurveValue(x);
        }

        private double OnCrackValueFunction(double t)
        {
            Point x = new Point(InnerCurve.GetX(t), InnerCurve.GetY(t));
            return _testData.OnCrackCurveValue(x);
        }

        private double OuterRadialFunction(double t)
        {
            return Radius;
        }

        private double InnerXFuntion(double t)
        {
            return Math.Cos(t); // TODO
        }

        private double InnerXFuntionDerivetive(double t)
        {
            return 1;
            //return - Math.Sin(t);
        }

        private double InnerYFunction(double t)
        {
            return Math.Cos(t); // TODO
        }

        private double InnerYFunctionDerivetive(double t)
        {
            return 1;
            //return - Math.Sin(t);
        }

        
        #endregion

        private double CoreFunction(double t, double tau)
        {
            return - Integral.CoefficientForWeakSingular(t, PointsNumber, tau) 
                + H(t, tau)* 2.0 * Math.PI 
                / (PointsNumber);
        }

        private double H(double t, double tau)
        {
            if(Math.Abs(t - tau) > DeviationEps)
            {
                double curveDeviation = Math.Pow(InnerCurve.GetX(t) - InnerCurve.GetX(tau), 2) 
                    + Math.Pow(InnerCurve.GetY(t) - InnerCurve.GetY(tau), 2);

                double denominator = Math.E * Math.E * curveDeviation;

                double numerator = 4.0 * Math.Pow(Math.Cos(t) - Math.Cos(tau), 2);

                return Math.Log(numerator / denominator)/(4.0 * Math.PI) 
                    + GrinFunction(t, tau);
            }
            else
            {
                double curveAbs = Math.Pow(InnerCurve.GetDerivetiveX(t), 2)
                    + Math.Pow(InnerCurve.GetDerivetiveY(t), 2);
                double logParam = 4.0 / (Math.E * Math.E * curveAbs);

                double result = Math.Log(logParam)/(4.0 * Math.PI)
                    + GrinFunction(t, t);

                return result;
            }            
        }

        private double GrinFunction(double t, double tau)
        {
            Point pointX = new Point(InnerCurve.GetX(t), InnerCurve.GetY(t));
            Point pointY = new Point(InnerCurve.GetX(tau), InnerCurve.GetY(tau));

            return GrinFunction(pointX, pointY);
        }

        private double GrinFunction(Point x, Point y)
        {
            double xAbs = Math.Pow(x.X, 2) + Math.Pow(x.Y, 2);
            double yAbs = Math.Pow(y.X, 2) + Math.Pow(y.Y, 2);
            double scalarProduct = x.X * y.X + y.Y * x.Y;
            double logParam = (Math.Pow(Radius, 4) + xAbs * yAbs - 2.0 * Radius * Radius * scalarProduct)
                / Math.Pow(Radius, 2);

            return Math.Log(logParam) / (4.0 * Math.PI);
        }

        private double DirectProblemRightPartFunction(double tParam)
        {
            DoubleCore<double> core = new DoubleCore<double>((t, tau) => {
                    Point pointX = new Point(InnerCurve.GetX(t), InnerCurve.GetY(t));
                    return RightPartCore(pointX, tau);
                });

            core.Prepare(tParam);
            return 2.0 * Integral.CalculateWithTrapeziumMethod(core, PointsNumber) + OnCrackValueFunction(tParam);
        }

        /// <summary>
        /// Is used to build solution
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double PartialSolution(Point x)
        {
            DoubleCore<Point> core = new DoubleCore<Point>(RightPartCore);
            core.Prepare(x);
            return - 2.0 * Integral.CalculateWithTrapeziumMethod(core, PointsNumber);
        }

        private double RightPartCore(Point pointX, double tau)
        {
            Point pointY = new Point(Radius * Math.Cos(tau), Radius * Math.Sin(tau));
            double scalarProduct = pointX.X * pointY.X + pointX.Y * pointY.Y;
            double yDoubleAbs = Radius * Radius;
            double xDoubleAbs = pointX.X * pointX.X + pointX.Y * pointX.Y;
            double deviationAbs = Math.Pow(pointX.X - pointY.X, 2) + Math.Pow(pointX.Y - pointY.Y, 2);

            double firstTerm = (scalarProduct - yDoubleAbs) / deviationAbs;
            double secondTerm = (yDoubleAbs * xDoubleAbs - Radius * Radius * scalarProduct) 
                / (Math.Pow(Radius, 4) + xDoubleAbs * yDoubleAbs - 2.0 * Radius * Radius * scalarProduct);

            // Radius = |y|
            double yAbs = Radius;
            return (firstTerm + secondTerm) / (2.0 * Math.PI * yAbs);
        }

        private double SolutionBuilderCoreFunction(Point toFindSolutionOn, Point yPoint)
        {
            double deviationAbs = Math.Sqrt(
                Math.Pow(toFindSolutionOn.X - yPoint.X, 2)
                + Math.Pow(toFindSolutionOn.Y - yPoint.Y, 2));

            double result = Math.Log(1.0 / deviationAbs) / (2.0 * Math.PI) 
                + GrinFunction(toFindSolutionOn, yPoint);
            return result;
        }

        public double BuildSolution(double[] density, Point toFindSolutionOn)
        {
            Func<double, double> coreFunction = (tau) =>
             {
                 Point yPoint = new Point(InnerCurve.GetX(tau), InnerCurve.GetY(tau));
                 return SolutionBuilderCoreFunction(toFindSolutionOn, yPoint);
             };

            return Integral.CalculateWithTrapeziumMethod(density, coreFunction)
                + PartialSolution(toFindSolutionOn);            
        }
    }
}
