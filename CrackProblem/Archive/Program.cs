//using System;
//using System.IO;

//namespace InverseDPForLE
//{
//    // Метод нелінійних інтегральних рівнянь для обернених задач теорії потенціалу
//    public class Program
//    {
//        private const double R = 2, InsideCurveRadius = 0.5, ApproxInsideCurveRadius =0.8, RForSolution = 1,ParameterOfRegul = 1e-2;        
//        private const double ExRadius = 10, ExAgle = 0; // exterior point   
//        private static void Main2(string[] args)
//        {
//            Console.WriteLine("Курсова робота \n Enter N (N*2 = number of points):");
//            int N = int.Parse(Console.ReadLine());
//            //DIRECT PROBLEM
//            //Problem pr = new Problem();
//            //pr.SolveDirectProblem(N, R, u_On_D2, InsideCurveRadius, u_On_D1, RForSolution);
//            //int n = 2 * N;
//            //Console.WriteLine("\n Solution on some curve: \n Accurate\t Received\t Acc Derivetive\t Derivetive ");
//            //for (int i = 0; i < n; i++)
//            //{
//            //    Console.WriteLine("{0:F8} \t {1:F8} \t {2:F8} \t {3:F8}", u_ForFoundedSolution(i * Math.PI / N), pr.ui[i],
//            //                                                                       g_Accurate(i * Math.PI / N), pr.gi[i]);
//            //}

//            // Запис в файл
//            //using (StreamWriter stream  = new StreamWriter("Derivative.txt",false))
//            //{
//            //    stream.WriteLine("Inner circle radius = {0}",InsideCurveRadius);
//            //    stream.WriteLine("Outer circle radius = {0}", R);
//            //    stream.WriteLine("u = {0}  on Г0",0);
//            //    stream.WriteLine("u = {0} on Г1", u_On_D2(0));
//            //    stream.WriteLine("n = {0}", N);
//            //    stream.WriteLine("   t   \t        g(t) on Г1   ");
//            //    for (int i = 0; i < n; i++)
//            //    {
//            //        stream.WriteLine("{0:F8} \t {1:F8}",  i * Math.PI / N, pr.gi[i]);
//            //    }
//            //}

//            //double maxDeviationDerivetive = Math.Abs(pr.gi[0] - g_Accurate(0));
//            //double maxDeviation = Math.Abs(pr.ui[0] - u_ForFoundedSolution(0));
//            //for (int i = 1; i < 2 * N; i++)
//            //{
//            //    maxDeviation = Math.Max(maxDeviation, Math.Abs(pr.ui[i] - u_ForFoundedSolution(i * Math.PI / N)));
//            //    maxDeviationDerivetive = Math.Max(maxDeviationDerivetive, Math.Abs(pr.gi[i] - g_Accurate(i * Math.PI / N)));
//            //}
//            //Console.WriteLine("\n Max deviation of solution = {0}", maxDeviation);
//            //Console.WriteLine("\n Max deviation of derivetive = {0}", maxDeviationDerivetive);

//            //// INVERSE PROBLEM
//            //double paramBest = 0;
//            //double maxPrev = 100;
//            //for (double par = 1e-5; par > 1e-16; par /= 10)
//            //{
//            //    Problem pr = new Problem();
//            //    pr.SolveInverseProblem(N, R, u_On_D2, InsideCurve, ApproxInsideCurveRadius, par);
//            //    int n = 2 * N;
//            //    double maxDev = 0;
//            //    //Console.WriteLine("\nr accurate \t r aproximate");
//            //    for (int i = 0; i < n; i++)
//            //    {
//            //        maxDev = Math.Max(maxDev, Math.Abs(InsideCurve(i * Math.PI / N) - pr.RadiusOfInsideCurve.Value(i * Math.PI / N)));
//            //        //Console.WriteLine("{0:F8} \t {1:F8} ", InsideCurve(i * Math.PI / N), pr.RadiusOfInsideCurve.Value(i * Math.PI / N));
//            //    }
//            //    if (maxPrev > maxDev)
//            //    {
//            //        paramBest = par;
//            //        maxPrev = maxDev;
//            //    }
//            //    //Console.WriteLine("Max deviation = {0}", maxDev);
//            //}
//            Problem pr = new Problem();
//            pr.SolveInverseProblem(N, R, u_On_D2, InsideCurve, ApproxInsideCurveRadius, ParameterOfRegul);
//            int n = 2 * N;
//            double maxDev = 0;
//            Console.WriteLine("\nr accurate \t r aproximate");
//            for (int i = 0; i < n; i++)
//            {
//                //maxDev = Math.Max(maxDev, Math.Abs(InsideCurve(i * Math.PI / N) - pr.RadiusOfInsideCurve.Value(i * Math.PI / N)));
//                //Console.WriteLine("{0:F8} \t {1:F8} ", InsideCurve(i * Math.PI / N), pr.RadiusOfInsideCurve.Value(i * Math.PI / N));
//            }

//            //Func<double, double> upFunc = (x) => Math.Pow(InsideCurve(x) - pr.RadiusOfInsideCurve.Value(x),2);

//            //Console.WriteLine("Best param = {0}",paramBest);
//            Console.WriteLine("Deviation  = {0}", maxDev);

//            Console.ReadKey();
//        }



//        public static double u_On_D2(double tx)
//        {
//            //Direct Problem
//            // Приклад 1 фундаментальний розвязок ln(1/|x-y|)/(2Pi)
//            //return FundamentalSolution(R,tx);

//            //Приклад 2  константа
//            return  1;

//            // Приклад 3  розвязок ln(|x-y|)
//            //return FundamentalSolutionFirst(R, tx);

//            //Приклад 4  x1^2 - x2^2
//            //return R * R * Math.Cos(2*tx);

//            // Inverse problem
//           // return 2;
//            //return 3 * Math.Sqrt(Math.Pow(Math.Cos(tx),2) + 0.25 * Math.Pow(Math.Sin(tx),2)) / 4;
//            //return Math.Exp(-Math.Pow(Math.Cos(tx),2));
//        }
//        static public double u_On_D1(double tx)
//        {
//            //Direct Problem
//            // Приклад 1 фундаментальний розвязок ln(1/|x-y|)/(2Pi)
//            //return FundamentalSolution(r0, tx);

//            //Приклад 2  константа
//            return 0; 

//            //Приклад 3  розвязок ln(|x-y|)
//            //return FundamentalSolutionFirst(r0,tx);

//            //Приклад 4  x1^2 - x2^2
//            //return r0 * r0 * Math.Cos(2*tx);

//            // Inverse Problem
//           // return 0;
//        }
//        public static double u_ForFoundedSolution(double tx)
//        {            
//            //Direct Problem
//            // Приклад 1 фундаментальний розвязок ln(1/|x-y|)/(2Pi)
//           // return FundamentalSolution(RForSolution,tx);

//            //Приклад 2  константа
//            return 1;

//            //Приклад 3  розвязок ln(|x-y|)
//            //return FundamentalSolutionFirst(r_for_solution, tx);

//            //Приклад 4  x1^2 - x2^2
//            //return r_f * r_f * Math.Cos(2*tx);             
//        } 
//        public static double g_Accurate( double tx )
//        {
//            //Direct Problem
//            // Приклад 1 фундаментальний розвязок ln(1/|x-y|)/(2Pi)
//            //return -(R - ExRadius * Math.Cos(tx - ExAgle)) /((2.0*Math.PI)* (R * R + ExRadius * ExRadius - 2 * R * ExRadius * Math.Cos(tx - ExAgle)));

//            //Приклад 2  константа
//            return 0;

//            //Приклад 3  розвязок ln(|x-y|)
//            //return (R - ExRadius * Math.Cos(tx - ExAgle)) / (R * R + ExRadius * ExRadius - 2 * R * ExRadius * Math.Cos(tx - ExAgle));

//            //Приклад 4  x1^2 - x2^2
//            //return 2.0 * R * Math.Cos(2.0*tx);
//        }
//        public static double InsideCurve(double t)
//        {
//            //Direct problem
//            //return InsideCurveRadius;
//            return Math.Sqrt(Math.Pow(Math.Cos(t),2) + 0.25* Math.Pow(Math.Sin(t),2));

//            // Inverse problem
//            //return InsideCurveRadius;
//            //return 3.0 * Math.Sqrt(Math.Pow(Math.Cos(t), 2) + 0.25 * Math.Pow(Math.Sin(t), 2)) / 4.0; 
//        }
//       private static double FundamentalSolution(double rx, double tx)
//        {
//            return Math.Log(1.0/Math.Sqrt(Math.Pow(rx * Math.Cos(tx) - ExRadius * Math.Cos(ExAgle), 2) + Math.Pow(rx * Math.Sin(tx) - ExRadius * Math.Sin(ExAgle), 2)))/(2.0*Math.PI);
//        }
//        private static double FundamentalSolutionFirst(double rx, double tx)
//        {
//            return Math.Log(Math.Sqrt(Math.Pow(rx * Math.Cos(tx) - ExRadius * Math.Cos(ExAgle), 2) + Math.Pow(rx * Math.Sin(tx) - ExRadius * Math.Sin(ExAgle), 2)));
//        }
//    }
//}
