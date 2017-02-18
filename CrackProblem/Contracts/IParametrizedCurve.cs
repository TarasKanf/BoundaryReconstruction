using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Contracts
{
    public interface IParametrizedCurve
    {
        double GetX(double t);
        double GetY(double t);
    }
}
