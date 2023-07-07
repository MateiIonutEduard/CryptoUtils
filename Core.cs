using Eduard;
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
        static EllipticCurve curve;

        public static void GenDomain(int bits, bool Search, bool Atkin=false)
        {
            curve = new EllipticCurve(bits);
            search = Search;
            atkin = Atkin;
        }

        static void mul(int p, int qnr, int x, int y, ref int a, ref int b)
        {
            int olda = a;
            a = (a * x + b * y * qnr) % p;
            b = (olda * y + b * x) % p;
        }

        static void quad(int p, int qnr, int x, int y, int e, ref int a, ref int b)
        {
            int k = e;
            a = 1;
            b = 0;
            if (k == 0) return;

            for (; ; )
            {
                if (k % 2 != 0)
                    mul(p, qnr, x, y, ref a, ref b);

                k >>= 1;

                if (k == 0) return;
                mul(p, qnr, x, y, ref x, ref y);
            }
        }

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
