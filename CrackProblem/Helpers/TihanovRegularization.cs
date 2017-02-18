using CrackProblem.LinearSystem;
using System;
using System.Linq;

namespace CrackProblem.Helpers
{
    public class TihanovRegularization
    {          
        private double[,] At;
        private int m, n; // number of rows and columns
        private double[] Atb;       
        //матриця А - розмірності m Х n, b - vector of values, lamda - parameter of regularisation
        public TihanovRegularization(double[,] A, double[] b, int _m, int _n) 
        {
            n = _n;
            m = _m;
            if (_m < _n) throw new Exception("Number of rows in matrix must be not lesser then number of columns.");
            if (_m != b.Count()) throw new Exception("Number of values must be the same as number of rows in matrix");            
            // будуємо матрицю Ат * А  де Ат - транспонована матриця А
            At = new double[n, n];
            double temp;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    temp = 0;
                    for (int k = 0; k < m; k++)
                    {
                        temp += A[k, i] * A[k, j];
                    }
                    At[i, j] = temp;
                }
            }
            Atb = new double[n];
            for (int i = 0; i < n; i++)
            {
                temp = 0;
                for (int k = 0; k < m; k++)
                    temp += A[k, i] * b[k];
                Atb[i] = temp;
                // Since Y0[i] = 0 in our case, we comment next row
                //Atb[i] += lamda * Y0[i];
            }
        }
        // матриця А - розмірності m Х n, b - vector of values, lamda - parameter of regularisation
        public double[] Solve(double lamda) 
        {
            // solve equation (At*A + lamda*I)*Y = At*b;           

            for (int i = 0; i < n; i++)
            {
                At[i, i] += lamda;
            }
            //будуємо вектор Ат * b + lamda*Y0 , де Y0 - початкове наближення до розвязку Y

            // розв`язуємо систему методом Якобі  
            //solution = SystemOfEquations.SolveWithIakobiMethod(At,Atb,0.001);
            SystemOfLinearEquations systemSolver = new SystemOfLinearEquations();

            return systemSolver.SolveWithQRmethod(At, Atb, Atb.Length);         
        }
       
    }
}
