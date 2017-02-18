using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using  System.IO;

namespace CrackProblem.Helpers
{
    public enum WriteMode
    {
        Console,
        File
    }

    public static class Printer
    {
        public static bool Append = false;
        private static string FileName = "Results.txt";

        public static WriteMode Mode = WriteMode.File;

        public static void WriteLine(string text)
        {
            if (Mode == WriteMode.File)
            {
                using (StreamWriter stream = new StreamWriter(FileName, Append))
                {
                    stream.WriteLine(text);
                }

                Append = true;
            }

            if (Mode == WriteMode.Console)
            {
                Console.WriteLine(text);
            }
           
        }

        public static void Write(double[] array)
        {
            if(Mode == WriteMode.File)
            {
                using (StreamWriter stream = new StreamWriter(FileName, Append))
                {
                    for (int i = 0; i < array.Length; i++)
                    {
                        stream.Write("{0:F8}  ", array[i]);
                        if ((i + 1) % 4 == 0) stream.WriteLine();
                    }

                    stream.WriteLine();
                }

                Append = true;
            }
            
            if( Mode == WriteMode.Console)
            {
                for (int i = 0; i < array.Length; i++)
                {
                    Console.Write("{0:F8}  ", array[i]);
                    if ((i + 1) % 4 == 0) Console.WriteLine();
                }

                Console.WriteLine();
            }
        }
    }
}
