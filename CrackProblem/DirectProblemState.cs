using System.Runtime.Remoting.Messaging;

namespace CrackProblem
{
    using Contracts;
    using Integrals;
    using CrackProblem.Helpers;
    using System;

    public class DirectProblemState
    {
        protected double DeviationEps = 1e-7;
        protected double Radius;
        protected StarCurve OuterCurve;
        public IParametrizedCurve InnerCurve { get; set; }
        public DoubleCore<double> DirectProblemCore;
        public Func<double,double> RightPartFunction { get; }

        /// <summary>
        /// Real number of points/ Notice <see cref="PointsNumber"/> is twice bigger than n in trigonometric quadratures
        /// </summary>
        public int PointsNumber { get; set; }
        protected IDirectProblemTestData DirectProblemTestData { get; set; }

        public DirectProblemState(double radius, int pointsNumber, IDirectProblemTestData directProblemTestData)
        {
            Radius = radius;
            OuterCurve = new StarCurve(OuterRadialFunction);
            InnerCurve = new ParametrizedCurve(InnerXFuntion, InnerYFunction, 
                InnerXFuntionDerivetive, InnerYFunctionDerivetive);
            PointsNumber = pointsNumber;
            DirectProblemCore = new DoubleCore<double>(CoreFunction);
            RightPartFunction = DirectProblemRightPartFunction;
            DirectProblemTestData = directProblemTestData;
        }

        #region Initails
        private double OnEdgeValueFunction(double t)
        {
            //Point x = new Point(OuterCurve.GetX(t), OuterCurve.GetY(t));
            return DirectProblemTestData.OnOuterCurveValue(t);
        }

        private double OnCrackValueFunction(double t)
        {
            Point x = new Point(InnerCurve.GetX(t), InnerCurve.GetY(t));
            return DirectProblemTestData.OnCrackCurveValue(x);
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
                + H(t, tau) * 2.0 * Math.PI
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
            } else
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
            return 2.0 * Integral.CalculateWithTrapeziumMethod(core, PointsNumber) + 2.0 * OnCrackValueFunction(tParam);
        }

        private double RightPartCore(Point pointX, double tau)
        {
            //Printer.WriteLine(pointX.X.ToString() + " /// " + pointX.Y.ToString());
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
            double result = Radius * (firstTerm + secondTerm) / (2.0 * Math.PI * yAbs);
            return result * OnEdgeValueFunction(tau);
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
            return -Integral.CalculateWithTrapeziumMethod(core, PointsNumber);
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

            return Integral.CalculateWithTrapeziumMethod(density, coreFunction) / 2.0
                + PartialSolution(toFindSolutionOn);
        }
        
        // Functions to calculate solution derivative on outer curve
        protected double Omega1(double tx)
        {
            DoubleCore<double> core = new DoubleCore<double>(Omega1CoreNotSingular);
            core.Prepare(tx);
            return - Integral.CalculateWithHyperSingularCore(core, PointsNumber); // перевірити обчисленн гіперсинугялрного інтегралу
        }

        protected double DataEquationOperatorCore(Point pointX, double tau)
        {
            Point pointY = new Point(InnerCurve.GetX(tau), InnerCurve.GetY(tau));
            return DataEquationOperatorCore(pointX, pointY);
        }

        protected double DataEquationOperatorCore(Point pointX, Point pointY)
        {
            double scalarProduct = pointX.X * pointY.X + pointX.Y * pointY.Y;

            double xDoubleAbs = Radius * Radius;
            double yDoubleAbs = pointY.X * pointY.X + pointY.Y * pointY.Y;

            double deviationAbs = Math.Pow(pointX.X - pointY.X, 2) + Math.Pow(pointX.Y - pointY.Y, 2);

            double firstTerm = (scalarProduct - xDoubleAbs) / deviationAbs;

            double secondTerm = (yDoubleAbs * xDoubleAbs - Radius * Radius * scalarProduct)
                / (Math.Pow(Radius, 4) + xDoubleAbs * yDoubleAbs - 2.0 * Radius * Radius * scalarProduct);
            // Radius = |y|
            double xAbs = Radius;
            double result = (firstTerm + secondTerm) / (4.0 * Math.PI * xAbs); // похідну по нормалі поділено на два
            return result;
        }

        private double Omega1CoreNotSingular(double tx, double ty)
        {
            return RightPartFunction(ty) / (2.0 * Math.PI * Radius);
        }

        public double[] BuildSolutionDerivativeOnOuterCurve(double[] density)
        {
            var descretePoints = IntegralEquationDiscretezer.GetDiscretePoints(PointsNumber);
            var solutionDerivative = new double[descretePoints.Length];
            for (int i = 0; i < descretePoints.Length; i++)
            {

                var core = new DoubleCore<Point>(DataEquationOperatorCore);
                core.Prepare(new Point(
                    OuterCurve.GetX(descretePoints[i]),
                    OuterCurve.GetY(descretePoints[i])));

                solutionDerivative[i] = Integral.CalculateWithTrapeziumMethod(density, core)
                + Omega1(descretePoints[i]);
            }

            return solutionDerivative;
        }
    }
}
