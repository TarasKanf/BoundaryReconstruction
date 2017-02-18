using System;

namespace CrackProblem.LinearSystem
{
    public class SystemOfLinearEquations
    {
        private const double epsilon = 0.000001;
        public double[] JacobiMethodSolving(double[,] matrix, double[] rightPart, double eps)
        {
            double[] xk1 = new double[rightPart.Length];
            double[] xk2 = new double[rightPart.Length];
            for (int i = 0; i < rightPart.Length; i++)
            {
                xk1[i] = 0;
                xk2[i] = 1;
            }

            double temp = 0;
            int iteration = 0;
            while ((LoopContinue(xk1, xk2, eps)) && (iteration < 10000))
            {
                iteration++;
                for (int i = 0; i < xk1.Length; i++)
                    xk1[i] = xk2[i];
                for (int i = 0; i < rightPart.Length; i++)
                {
                    temp = 0;
                    for (int j = 0; j < rightPart.Length; j++)
                    {
                        if (j == i) continue;
                        temp += matrix[i, j] * xk1[j];
                    }
                    temp -= rightPart[i];
                    xk2[i] = -temp / matrix[i, i];
                }
            }
            if (iteration > 9999) throw new Exception("Jacobi method does not converge to solution.");
            return xk2;
        }
        private bool LoopContinue(double[] x1, double[] x2, double eps)
        {
            bool ans = true;
            double maxdif = Math.Abs(x1[0] - x2[0]);
            for (int i = 1; i < x1.Length; i++)
                if (maxdif < Math.Abs(x1[i] - x2[i])) maxdif = Math.Abs(x1[i] - x2[i]);
            if (maxdif < eps) ans = false;
            return ans;
        }

        public double[] LU_methodSolving(double[,] matrix, double[] rightPart, int n)
        {
            // розклад матриці A
            double[,] lu = new double[n, n];
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    sum = 0;
                    for (int k = 0; k < i; k++)
                        sum += lu[i, k] * lu[k, j];
                    lu[i, j] = matrix[i, j] - sum;
                }
                for (int j = i + 1; j < n; j++)
                {
                    sum = 0;
                    for (int k = 0; k < i; k++)
                        sum += lu[j, k] * lu[k, i];
                    lu[j, i] = (1 / lu[i, i]) * (matrix[j, i] - sum);
                }
            }
            double[] y = new double[n];
            for (int i = 0; i < n; i++)
            {
                sum = 0;
                for (int k = 0; k < i; k++)
                    sum += lu[i, k] * y[k];
                y[i] = rightPart[i] - sum;
            }
            double[] x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                sum = 0;
                for (int k = i + 1; k < n; k++)
                    sum += lu[i, k] * x[k];
                x[i] = (1 / lu[i, i]) * (y[i] - sum);
            }
            return x;
        }
        public double[] SolveWithQRmethod(double[,] matrix, double[] rightPart, int n)
        {
            double[,] M = (double[,])matrix.Clone();
            double[] M1;
            double[] M2;
            //bool IsDegenerate;            
            M1 = new double[n];
            M2 = new double[n];          
            bool sing = Factoring(ref M, ref M1, ref M2, n); // is deganarated => sing == true
            QRSolve(ref M, ref M1, ref M2, n, ref rightPart);
            return rightPart;
        }        
        private bool Factoring(ref double[,] M, ref double[] M1, ref double[] M2, int N)
        {
            bool sing = false; // матриця невироджена
            for (int k = 0; k < N - 1; k++)
            {
                double eta = Math.Abs(M[N - 1, k]);

                for (int i = k; i < N - 1; i++)
                {
                    if (Math.Abs(M[i, k]) > eta)
                    {
                        eta = Math.Abs(M[i, k]);
                    }
                }

                if (Math.Abs(eta) == 0)
                {
                    // вироджена                    
                    M1[k] = 0;
                    M2[k] = 0;
                    sing = true;
                }
                else
                {
                    for (int i = k; i < N; i++)
                    {
                        M[i, k] /= eta;
                    }
                    double sum = 0;
                    for (int j = k; j < N; j++)
                    {
                        sum += M[j, k] * M[j, k];
                    }
                    double sigma = ((M[k, k] == 0) ? Math.Sqrt(sum) : Math.Sign(M[k, k])) * Math.Sqrt(sum);
                    //double sigma = Math.Sign(M[k, k]) * Math.Sqrt(sum);
                    M[k, k] += sigma;
                    M1[k] = sigma * M[k, k];
                    M2[k] = -1 * eta * sigma;
                    for (int j = k + 1; j < N; j++)
                    {
                        sum = 0;
                        for (int l = k; l < N; l++)
                        {
                            sum += M[l, k] * M[l, j];
                        }
                        double tao = sum / M1[k];
                        for (int i = k; i < N; i++)
                        {
                            M[i, j] = M[i, j] - tao * M[i, k];
                        }
                    }
                }
            }
            if (Math.Abs(M[N - 1, N - 1]) < epsilon) sing = true;
            M2[N - 1] = M[N - 1, N - 1];
            return sing;
        }
        private void QRSolve(ref double[,] M, ref double[] M1, ref double[] M2, int N, ref double[] b)
        {
            // b <= QT*b
            double t;
            for (int j = 0; j < N - 1; j++)
            {
                t = 0;
                for (int i = j; i < N; i++)
                {
                    t += M[i, j] * b[i];
                }
                t /= M1[j];
                for (int i = j; i < N; i++)
                {
                    b[i] -= t * M[i, j];
                }
            }
            RSolve(ref M, ref M1, ref M2, N, ref b);
        }

        private void RSolve(ref double[,] M, ref double[] M1, ref double[] M2, int N, ref double[] b) //  розвязок зберігають в векторі b
        {
            // Rx = b  -> x
            b[N - 1] /= M2[N - 1];
            double sum = 0;
            for (int i = N - 2; i >= 0; i--)
            {
                sum = 0;
                for (int j = i + 1; j < N; j++)
                {
                    sum += M[i, j] * b[j];
                }
                b[i] = (b[i] - sum) / M2[i];
            }
        }
    }
}
