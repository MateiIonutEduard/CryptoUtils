using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CryptoUtils
{
    public class Core
    {
        static bool atkin = false;
        static bool search = false;
        static bool escape = false;

        static int gcd(int a, int b)
        {
            if (a == 1 || b == 1) return 1;

            while (b != 0)
            {
                int r = a % b;
                a = b;
                b = r;
            }

            return a;
        }

        static int phi(int n)
        {
            int r = 0;

            for (int c = 1; c < n; c++)
                if (gcd(c, n) == 1) r++;

            return r;
        }
    }
}
