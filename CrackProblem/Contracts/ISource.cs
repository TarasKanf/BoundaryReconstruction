using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrackProblem.Contracts
{
    public interface ISource<T>
    {
        T Value(T x);
    }
}
