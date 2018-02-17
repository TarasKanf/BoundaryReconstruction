//using System;
//using InverseDPForLE.Helpers;
//using InverseDPForLE.Integrals;
//using Курсова__5_семестр_.Helpers;

//namespace InverseDPForLE
//{
//    public class Problem
//    {
        //#region Fields
        //// метод нелінійних інтегральних рівнянь для обернених задач теорії потенціалу
        //// delta(U) =0 in D
        //// U =0 on Г1
        //// d(U)/dV on Г2 = g 
        ////U on Г2 = f 
        //// Г2 - коло радіуса R
        //// знайти  Г1.  Г1 - ?
        //// Позначатимемо в програмі D1 == Г1,D2 == Г2   
        //private double R; // радіус зовнішнього кола
        //private Func<double,double> f; // U on Г2 = f (bound function);  f(t) t є [0,2*PI]
        //private Func<double, double> f0; // U on Г1 = f0;   f(t) є [0,2*PI]        
        //public TrigonPolinom RadiusOfInsideCurve => r;
        //private TrigonPolinom r; // задає радіальну функцію для внутрінього кола
        //private TrigonPolinom x_radius; // завжди задаватиме радіал. функцію тої кривої на якій лежить точка x 
        //private int N ; // кількість точок за якими будується наближення кривої Г1 і Г2 буде рівна 2*N
        //public double[] fi,gi; // густина, похідна, вектор для радіусів    
        
        ////Поля тільки для ПРЯМОЇ ЗАДАЧІ    
        //public double[] ui;
        //public double RForSolution, R0;// радіус кривої на якій шукатимемо розвязок; радіус внутрішнього кола
        
        //// Поля тільки для ЗВОРОТНЬОЇ ЗАДАЧІ
        //private Func<double, double> rD0; // початкове точне задання внутр кривої щоб знайти похідну по нормалі
        //public double[] correction; // корекція радільної функції для вн. кривої
        //public double parOfRegular; // параметр регуляризаціїї
        //private int correctionN = 8;
        //#endregion

        //#region Direct Problem
        //public void SolveDirectProblem(int _N, double _R, Func<double,double> _f, double _r0, Func<double,double> _f0, double r_f_s)
        //{
        //    RForSolution = r_f_s;
        //    N = _N;
        //    R = _R;
        //    R0 = _r0;
        //    f = _f;
        //    f0 = _f0;                      
        //    double[] radius = new double[2 * N];
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        radius[i] = R0;
        //    }
        //    r = new TrigonPolinom(radius, N);
        //    x_radius = r; // в рівнянні поля x на Г1 (вн. крива)                      
        //    IntegralEquation equation = new IntegralEquation( RightPartOfFieldEquation);
        //    SmoothCore s1 = new SmoothCore(GrinFuncWeakSingularPart);
        //    SmoothCore s2 = new SmoothCore(GrinFuncSmoothPart);
        //    fi = equation.SolveIE_WithWeakAndSmoothCore(s1,s2,N); // наближений розв`язок інтегрального рівняння поля в точках t[j] = j*PI/N ,  j = 0, 2*N -1   
        //    ui = FindSolution(fi);            
        //    gi = FindDerivetiveOnD2(fi);
        //}
        //#endregion

        //#region Inverse Problem
        //public void SolveInverseProblem(int _N, double _R, Func<double, double> _f, Func<double, double> _rD0,double rInitial,double parTihan)
        //{
        //    // Find normal derivetive on D2
        //    N = _N;
        //    R = _R;            
        //    f = _f;
        //    rD0 = _rD0;
        //    parOfRegular = parTihan;
        //    f0 = ( i => 0);
        //    double h = Math.PI / N;
        //    double temp = 0;
        //    double[] radius = new double[2 * N];
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        radius[i] = rD0(temp);                
        //        temp += h;
        //    }

        //    Printer.WriteLine($"R = {R}");
        //    Printer.WriteLine($"r_0 = {rInitial}");
        //    Printer.WriteLine($"f = {f(0)} on G2");
        //    Printer.WriteLine($"n = {N}");
        //    Printer.WriteLine($"alfa = {parTihan}");

        //    Printer.WriteLine(" a) ***************************************");
        //    Printer.WriteLine("r accurate");
        //    Printer.Write(radius);

        //    r = new TrigonPolinom(radius, N);
        //    x_radius = r;
        //    IntegralEquation equation = new IntegralEquation(RightPartOfFieldEquation);
        //    SmoothCore s1 = new SmoothCore(GrinFuncWeakSingularPart);
        //    SmoothCore s2 = new SmoothCore(GrinFuncSmoothPart);

        //    fi = equation.SolveIE_WithWeakAndSmoothCore(s1, s2, N);

        //    temp = 0;
        //    double[] ftemp = (double[])fi.Clone();
        //    for (int j = 0; j < fi.Length; j++)
        //    {
        //        ftemp[j] /= Jacobi(r, temp);
        //        temp += h;
        //    }
        //    Printer.WriteLine("Density on G1");
        //    Printer.Write(ftemp);  

        //    gi = FindDerivetiveOnD2(fi);// чому дорівнює -fi   і чи так має бути     
        //    Printer.WriteLine("Flux on G2");
        //    Printer.Write(gi);
        //    Printer.WriteLine("b) *************************************");

        //    // перевірка похідної Фреше
        //    Examine(fi, N);
             
        //    //ініціалізуємо початкове наближення до вн. кривої
        //    temp = 0;
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        radius[i] = rInitial;
        //        temp += h;
        //    }
        //    r = new TrigonPolinom(radius, N);
        //    // зовнішня крива
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        radius[i] = R;                
        //    }
        //    TrigonPolinom exteriorCurve = new TrigonPolinom(radius, N);
        //    int iter = 0;
        //    bool exitLoop = false;
        //    double myconst = 1, delta = 0.0001;
        //    correction = new double[2 * N];
        //    double maxPrevCorection = 100;
        //    double[] coefficients = new double[2*correctionN +1];
        //    for (int i = 0; (i < 50) && (!exitLoop); i++)
        //    {
        //        x_radius = r; // x на Г1
        //        equation = new IntegralEquation(RightPartOfFieldEquation);
        //        s1 = new SmoothCore(GrinFuncWeakSingularPart);
        //        s2 = new SmoothCore(GrinFuncSmoothPart);

        //        fi = equation.SolveIE_WithWeakAndSmoothCore(s1, s2, N);
        //        for (int j = 0; j < fi.Length; j++)
        //        {
        //            fi[j] /= Jacobi(r, j*h);
        //        }
        //        Printer.WriteLine($"interation № {i}");
        //        Printer.WriteLine("Density on G1");
        //        Printer.Write(fi);

        //        x_radius = exteriorCurve; // х на Г2
        //        var tr = GetTihanovRegularization(gi, fi);
        //        double[] coefficientsNext = tr.Solve(parOfRegular);

        //        for (int k = 0; k < 2 * correctionN + 1; k++)
        //        {
        //           coefficients[k] = coefficientsNext[k];
        //        }

        //        Printer.WriteLine("Розв'язки системи : ");
        //        Printer.Write(coefficients);

        //        double max = 0;
        //        for (int j = 0; j < correction.Length; j++)
        //        {
        //            correction[j] = 0;
        //        }
        //        for (int j = 0; j < correction.Length; j++)
        //        {
        //            for (int k = 0; k < 2*correctionN + 1; k++)
        //            {
        //                correction[j] = correction[j] + coefficients[k]*L(k, j*h);
        //            }
        //            max = Math.Max(Math.Abs(correction[j]), max);
        //        }

        //        Printer.WriteLine("Correction");
        //        Printer.Write(correction);
             
        //        Console.WriteLine("Correction {0} iter {1}", correction[0],iter);
        //        if (!r.Add(correction))
        //        {
        //            Console.WriteLine("Wrong adding of polinoms");
        //            Console.ReadKey();
        //            return;
        //        }

        //        var numerator = correction;
        //        var denominator = radius;
        //        //double EC = ExitCondition(up, down);
        //        exitLoop = (ExitCondition(numerator, denominator) < (myconst * delta));
        //        iter++;
        //        //
        //        temp = 0;
        //        for (int j = 0; j < 2 * N; j++)
        //        {
        //            radius[j] = r.Value(temp);
        //            //if (max > maxPrevCorection) exitLoop = true;
        //            maxPrevCorection = max;
        //            if (radius[j] < 0) exitLoop = true;
        //            temp += h;
        //        }
        //        //Console.WriteLine("Radius  {0}", radius[0]);
        //        //if(exitLoop) break;
                
        //    }
        //    Console.WriteLine("Number of iterations {0}",iter);           

        //}
        //public double ExitCondition( double[] up, double[] down)
        //{
        //    double sum = 0;
        //    for (int i = 0; i < 2*N; i++)
        //    {
        //        sum += Math.Pow(up[i],2);
        //    }
        //    sum *= Math.PI / N;
        //    double sumdown = 0;
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        sumdown += Math.Pow(down[i],2);
        //    }
        //    sumdown *= Math.PI / N;
        //    return Math.Sqrt(sum) / Math.Sqrt(sumdown);
        //}
        //#endregion

        //#region Functions

        //private double Jacobi(TrigonPolinom radialFunc ,double t)
        //{
        //    return Math.Sqrt(Math.Pow(radialFunc.Value(t),2) + Math.Pow(radialFunc.Derivative(t),2));
        //}
        //private double GrinFuncFirstSummand(double tx, double ty)
        //{
        //    double rx = x_radius.Value(tx),ry= r.Value(ty);                  
        //    return Math.Log((Math.Pow(R, 4) + Math.Pow(rx, 2) * Math.Pow(ry, 2) - 2 * R * R * rx * ry * Math.Cos(tx - ty))/(R*R)) / (4 * Math.PI);
        //}
        //private double GrinFuncWeakSingularPart(double tx, double ty)
        //{
        //    return -1.0 / (4.0 * Math.PI); // це множник, що стоїть перед особливістю ln((4/e) * sin^2((tx-ty)/2))
        //}
        //private double GrinFuncSmoothPart(double tx, double ty)
        //{
        //    if (Math.Abs(tx - ty) > 1e-7)
        //    {
        //        double rx = x_radius.Value(tx);
        //        double ry = r.Value(ty);
        //        return GrinFuncFirstSummand(tx, ty) 
        //            + Math.Log((4.0 * Math.Pow(Math.Sin((tx - ty) / 2.0), 2)) 
        //            / (Math.E * (rx * rx + ry * ry - 2 * rx * ry * Math.Cos(tx - ty)))) 
        //            / (4.0 * Math.PI);
        //    }
        //    return GrinFuncFirstSummand(tx, tx) 
        //        + Math.Log(1.0 / (Math.E * (Math.Pow(r.Derivative(tx), 2) 
        //        + Math.Pow(r.Value(tx), 2)))) / (4.0 * Math.PI);
        //}

        //private double RightPartOfFieldEquation(double tx)
        //{
        //    return f0(tx) - Omega(tx);
        //}        
        //private double Omega(double tx)
        //{            
        //    SmoothCore core = new SmoothCore(OmegaCore);
        //    core.Prepare(tx);
        //    return -Integral.CalculateWithTrapeziumMethod(core,N);
        //}             
        //private double OmegaCore(double tx, double ty)
        //{
        //    double rx = x_radius.Value(tx);          
        //    return f(ty)*((rx * rx - R * R) / (R * R + rx * rx - 2 * R * rx * Math.Cos(tx - ty))) / (2 * Math.PI);
        //}        
        //private double Omega1(double tx)
        //{            
        //    SmoothCore core = new SmoothCore(Omega1CoreNotSingular);
        //    core.Prepare(tx);            
        //    return -Integral.CalculateWithHyperSingularCore(core, N); // перевірити обчисленн гіперсинугялрного інтегралу
        //}
        //private double Omega1CoreNotSingular(double tx, double ty)
        //{
        //    return f(ty) / (2.0*Math.PI * R);
        //}
        
        //private double[] FindDerivetiveOnD2(double[] y)
        //{
        //    double h = Math.PI / N;
        //    double[] g = new double[2*N];
        //    int n = 2 * N;
        //    double xt = 0;
        //    for (int i = 0; i < n; i++)
        //    {
        //        double temp = 0;
        //        double sum = 0;
        //        for (int j = 0; j < n; j++)
        //        {
        //            sum += y[j] * Core_dGvx_OnD2(xt,temp);
        //            temp += h;
        //        }
        //        sum *= Math.PI / N;
        //        g[i] = sum + Omega1(xt);
        //        xt += h;
        //    }
        //    return g;
        //}
        //private double Core_dGvx_OnD2(double tx, double ty)
        //{
        //    double ry = r.Value(ty);
        //    return (( ry*ry - R * R ) / (Math.Pow(R,2) +  ry*ry - 2.0*R * ry * Math.Cos(tx - ty)))/(2.0*Math.PI*R); // перевірити в математиці
        //}        
        //private double[] FindSolution(double[] y)
        //{
        //    double h = Math.PI / N;
        //    double[] u = new double[2 * N];
        //    int n = 2 * N;
        //    double xt = 0;
        //    double[] radius = new double[n];
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        radius[i] = RForSolution;
        //    }
        //    TrigonPolinom prev = x_radius;
        //    x_radius = new TrigonPolinom(radius, N);
        //    for (int i = 0; i < n; i++)
        //    {
        //        double temp = 0;
        //        double sum = 0;
        //        for (int j = 0; j < n; j++)
        //        {
        //            sum += y[j] * CoreFullGrinFunc(xt, temp);
        //            temp += h;
        //        }
        //        sum *= Math.PI / N;                
        //        u[i] = sum + Omega(xt);
        //        xt += h;
        //    }
        //    x_radius = prev;
        //    return u;           
        //}
        //private double CoreFullGrinFunc(double tx,double ty)//
        //{
        //    double ry = r.Value(ty);
        //    double rx = x_radius.Value(tx);             
        //    return Math.Log(1.0 / Math.Sqrt(rx*rx + ry*ry - 2.0 * rx * ry * Math.Cos(tx - ty)))/(2.0*Math.PI) + GrinFuncFirstSummand(tx, ty);
        //}        
        //// find frechet derivetive, form the equation Ax = b from integral equation (2) and return an object to regularize with different parameters of regulazation
        //private TihanovRegularization GetTihanovRegularization(double[] g, double[] yi)
        //{           
        //    double h = Math.PI / N;
        //    int n = 2 * N;
        //    // обчислюємо праву частину
        //    double [] b = new double[n];
        //    double ti = 0, sj;
        //    for (int i = 0; i < n; i++)
        //    {               
        //        b[i] = g[i] - Omega1(ti);
        //        sj = 0;
        //        double sum =0;
        //        for (int j = 0; j < n; j++)
        //        {
        //            sum += yi[j] * Core_dGvx_OnD2(ti, sj) * Jacobi(r,sj);
        //            sj += h;
        //        }
        //        sum *= Math.PI / N;
        //        b[i] -= sum;
        //        ti += h;
        //    }
        //    //            
        //    // заповнюємо матрицю А[,]
        //    double[,] A = new double[n, 2*correctionN + 1];
        //    ti = 0;            
        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j <= 2*correctionN; j++)
        //        {
        //            double sum = 0;
        //            sj = 0;
        //            for (int k = 0; k < n; k++)
        //            {
        //                double radialFunc = r.Value(sj);
        //                double devRadialFunc = r.Derivative(sj);
        //                sum += fi[k]*
        //                       (FreshetCore(ti, sj)*Jacobi(r, sj)*L(j, sj) 
        //                       + Core_dGvx_OnD2(ti, sj)*(radialFunc*L(j, sj) 
        //                       + devRadialFunc*LDev(j, sj))/Jacobi(r, sj));
        //                sj += h;
        //            }
        //            A[i, j] = sum * Math.PI / N;
        //        }
        //        ti += h;
        //    }
        //    return new TihanovRegularization(A,b,n,2*correctionN + 1);                       
        //}

        //private double L(int k, double t)
        //{
        //    if (k < 0 || k > 2 * correctionN)
        //        throw new Exception("LDev index overflow");
        //    if (k <= correctionN)
        //        return Math.Cos(k*t);
        //    return Math.Sin((k - correctionN)*t);
        //}

        //private double LDev(int k, double t)
        //{
        //    if(k < 0 || k > 2*correctionN) 
        //        throw  new Exception("LDev index overflow");
        //    if (k <= correctionN + 1)
        //        return  - k*Math.Sin(k * t);
        //    return (k - correctionN)*Math.Cos((k - correctionN) * t);
        //}

        //private double FreshetCore(double tx,double ty)
        //{
        //    double ry = r.Value(ty);
        //    return  (2.0 * R * ry - (R * R + ry * ry) * Math.Cos(tx - ty)) / 
        //        (Math.PI  * Math.Pow(R * R + ry * ry - 2 * R * ry * Math.Cos(tx - ty), 2));            
        //}
        //#endregion

        //#region Examine freshet core
        //private double FreshetCore(double ry,double tx, double ty)// +
        //{            
        //    return (2.0 * R * ry - (R * R + ry * ry) * Math.Cos(tx - ty)) /
        //        (Math.PI * Math.Pow(R * R + ry * ry - 2 * R * ry * Math.Cos(tx - ty), 2));
        //}
        //private double Core_dGvx_OnD2(double ry,double tx, double ty)
        //{            
        //    return ((ry * ry - R * R) / (Math.Pow(R, 2) + ry * ry - 2.0 * R * ry * Math.Cos(tx - ty))) / (2.0 * Math.PI * R); // перевірити в математиці
        //}
        //private double[] Ar(double[] hVec,int N, double[] mu)
        //{
        //    double h = Math.PI / N;
        //    double[] result = new double[2 * N];
        //    double sum = 0;
        //    double ti = 0, sj = 0 ;           
        //    for (int i = 0; i < 2*N; i++)
        //    {
        //        sum = 0;
        //        sj = 0;
        //        for (int j = 0; j < 2*N; j++)
        //        {
        //            sum += FreshetCore(r.Value(sj),ti, sj) * hVec[j] * mu[j];
        //            sj += h;
        //        }
        //        result[i] = sum * Math.PI / (double)N;
        //        ti += h;
        //    }
        //    return result;
        //}
        //private double[] Arizn(double[] hVec, int N, double[] mu, double t)
        //{
        //    double h = Math.PI / N;
        //    double[] result = new double[2 * N];
        //    double sum = 0;
        //    double ti = 0, sj = 0;
        //    for (int i = 0; i < 2 * N; i++)
        //    {
        //        sum = 0;
        //        sj = 0;
        //        for (int j = 0; j < 2 * N; j++)
        //        {
        //            sum += (Core_dGvx_OnD2(r.Value(sj)+ hVec[j]*t, ti, sj)- Core_dGvx_OnD2(r.Value(sj), ti, sj)) * mu[j];
        //            sj += h;
        //        }
        //        result[i] = sum * Math.PI / (N*t);
        //        ti += h;
        //    }
        //    return result;
        //}
        //private double MaxNorm(double[] vec)
        //{
        //    double result = Math.Abs(vec[0]);
        //    for(int i = 0; i < vec.Length; i++)
        //    {
        //        result = Math.Max(Math.Abs(vec[i]),result);
        //    }
        //    return result;
        //}
        //private void Examine(double[] mu, int N)
        //{
        //    double[] h = new double[2*N];
        //    // будь-який вектор h
        //    Random random = new Random();
        //    for (int i = 0; i < 2 * N; i++)
        //        h[i] = random.NextDouble();
        //    double t;

        //    for(int i = -1; i > -5; i--)
        //    {
        //        t = Math.Pow(10, i);
        //        double[] dev1 = Ar(h, N, mu);
        //        double[] dev2 = Arizn(h, N, mu, t);
        //        for(int k = 0;k< dev1.Length; k++)
        //        {
        //            dev1[k] -= dev2[k];
        //        }
        //        double maxnorm = MaxNorm(dev1);
        //        Console.WriteLine("Maxnorm divation of Freshet core = {0}", maxnorm);
        //    }
           
        //}
        //#endregion
//    }
//}
