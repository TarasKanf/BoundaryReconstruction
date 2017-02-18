using CrackProblem.Contracts;
using System;

namespace CrackProblem.Integrals
{ 
    public class DoubleCore<ParamT> : ICore
    {
        public ParamT Param { get; private set; }
        private Func<ParamT, double,double> F;
        public DoubleCore(Func<ParamT, double, double> f)
        {
            F = f;
        }
        public void Prepare(ParamT par)
        {
            Param = par;
        }
        public double GetValue(double t)
        {
            return F(Param,t);
        }
    }
}
