using Eduard;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
#pragma warning disable

namespace CryptoUtils
{
    public class Core
    {
        static bool atkin = false;
        static bool search = false;

        static List<ModularPolynomial> Gl;
        static List<int> prime;

        static bool escape = false;
        static EllipticCurve curve;
        static StreamReader sr;

        public static void GenDomain(int bits, bool Search, bool Atkin=false)
        {
            sr = new StreamReader("mueller.raw");
            curve = new EllipticCurve(bits);
            ReadModularPolynomialsByLimit(curve.field);

            search = Search;
            atkin = Atkin;
        }

        public static EllipticCurve GetEllipticCurve()
        {
            if(search)
            {
                BigInteger lim = 2;
                BigInteger a = curve.a;
                BigInteger b = curve.b;

                BigInteger field = curve.field;
                BigInteger Q4 = (field >> 64).Sqrt();

                if (field.GetBits() > 256) 
                    Q4 = (field >> 72).Sqrt();

                List<int> u = new List<int>();
                List<int> v = new List<int>();

                bool found = false;
                int index = 0;

                while(!found)
                {
                    Util.modulo(field);
                    Field.modulo(field);
                    Polynomial.SetField(field);

                    Polynomial X = new Polynomial(1, 0);
                    Polynomial Y2 = new Polynomial(1, 0, a, b);

                    Polynomial XP = Polynomial.Pow(X, field, Y2);
                    Polynomial pp = Polynomial.Gcd(XP - X, Y2);

                    int t = (pp == 1) ? 1 : 0;

                    if(t == 0)
                    {
                        b++;
                        continue;
                    }

                    u.Add(t);
                    v.Add(2);

                    Field delta = -16 * (4 * a * a * a + 27 * b * b);
                    Field j = (-1728 * 64 * a * a * a) / delta;
                    Field deltal = 0;

                    BigInteger tau = 0;

                    ModularPolynomial dGx, dGy, dGxx, dGxy, dGyy;
                    Field Eg, Ej, Exy, Dg, Dj, E4b, E6b, p1;
                    Field atilde, btilde, jl, E4bl, E2bs, gd, jd, E0b;
                    Field Dgd, Djd, E0bd, f, fd, Dgs, Djs, jld, E6bl;
                    int s, el, count;

                    E4b = (Field)(-a) / 3;
                    E6b = (Field)(-b) / 2;
                    delta = (E4b * E4b * E4b - E6b * E6b) / 1728;

                    while(index < Gl.Count && lim < Q4)
                    {
                        int p = prime[index];

                        for (s = 1; ; s++)
                            if (s * (p - 1) % 12 == 0) break;

                        Polynomial mod = Gl[index].F((BigInteger)j);
                        XP = Polynomial.Pow(X, field, mod);

                        Polynomial G = Polynomial.Gcd(XP - X, mod);
                        if (G.Degree == p + 1) continue;

                        if(G.Degree == 0)
                        {
                            // Atkin case
                            PolyMod.SetModulus(mod);
                            PolyMod[] m = new PolyMod[p + 1];
                            int max_r = p + 1, l = 0, r = 2;

                            m[0] = XP;
                            bool useful = false;

                            int k = (int)(field % p);
                            m[1] = PolyMod.Compose(m[0], m[0]);
                            if (!atkin) max_r = 2;

                            // Compute Frobenius map
                            for (r = 2; r <= max_r; r++)
                            {
                                PolyMod C = 0;
                                int js = (p + 1) / r;
                                int sym = BigInteger.Jacobi(k, p);

                                if ((p + 1) % r != 0) continue;
                                if (js % 2 == 0 && sym == -1) continue;
                                if (js % 2 == 1 && sym == 1) continue;

                                int kk = r, i = 0;
                                bool first = true;

                                while (kk > 0)
                                {
                                    if ((kk & 1) == 1)
                                    {
                                        if (first) C = m[i];
                                        else C = PolyMod.Compose(m[i], C);
                                        first = false;
                                    }

                                    kk >>= 1;
                                    if (kk == 0) break;
                                    i++;

                                    if (i > l)
                                    {
                                        m[i] = PolyMod.Compose(m[i - 1], m[i - 1]);
                                        l = i;
                                    }
                                }

                                PolyMod ok = C - X;

                                if (ok == 0)
                                {
                                    useful = true;
                                    break;
                                }
                            }

                            if (!useful) continue;
                            int ax = 0, by = 0, candidates, gx, gy, ord, qnr = 2;
                            bool gen = false;

                            while (BigInteger.Jacobi(qnr, p) != -1)
                                qnr++;

                            ord = p * p - 1;
                            gy = 1;

                            for (gx = 1; gx < p; gx++)
                            {
                                gen = true;

                                for (int jj = 2; jj <= ord / 2; jj++)
                                {
                                    if (ord % jj != 0) continue;
                                    quad(p, qnr, gx, gy, ord / jj, ref ax, ref by);

                                    if (ax == 1 && by == 0)
                                    {
                                        gen = false;
                                        break;
                                    }
                                }

                                if (gen) break;
                            }

                            candidates = 0;
                            lim *= p;

                            for (int jj = 1; jj < r; jj++)
                            {
                                if (jj > 1 && gcd(jj, r) != 1) continue;
                                quad(p, qnr, gx, gy, jj * ord / r, ref ax, ref by);

                                BigInteger two = 2;
                                tau = ((ax + 1) * k * (int)two.Inverse(p)) % p;

                                if (tau == 0)
                                {
                                    if ((field + 1 - tau) % p == 0)
                                    {
                                        if (search)
                                            escape = true;
                                    }

                                    u.Add((int)tau);
                                    v.Add(p);
                                    break;
                                }
                                else if (BigInteger.Jacobi(tau, p) == 1)
                                {
                                    tau = Eduard.Core.sqrmp((uint)tau, (uint)p);
                                    tau = (2 * tau) % p;
                                    if (candidates == phi(r)) break;
                                }
                            }

                            if (escape) break;
                            continue;
                        }
                        else
                        {

                        }
                    }
                }
            }
            else
            {

            }
            return null;
        }

        public static BigInteger Kangaroo(BigInteger a, BigInteger b, BigInteger p, BigInteger order, BigInteger ordermod)
        {
            ECPoint[] wild = new ECPoint[80], tame = new ECPoint[80];
            BigInteger[] wdist = new BigInteger[80], tdist = new BigInteger[80];
            int[] wname = new int[80], tname = new int[80];

            EllipticCurve curve = new EllipticCurve(a, b, p, order);
            ECPoint ZERO = ECPoint.POINT_INFINITY;
            ECPoint[] K = new ECPoint[10], TE = new ECPoint[10];
            ECPoint X, P, G, trap;

            ECPoint[] table = new ECPoint[128];
            BigInteger[] start = new BigInteger[10];
            BigInteger txc, wxc, mean, leaps, upper, lower, middle, x, y, n, w, t, nrp;
            int i, jj, j, m, sp = 1, nw, nt, cw, ct, k, distinguished, nbits;
            BigInteger[] D = new BigInteger[10];

            BigInteger s = 0, real_order;
            BigInteger[] distance = new BigInteger[128];
            bool bad, collision, abort;

            Console.WriteLine("\nReleasing 5 Tame and 5 Wild Kangaroos!");

            for (; ; )
            {
                // find a random point on the curve
                P = curve.BasePoint;
                lower = p + 1 - 2 * p.Sqrt() - 3; // lower limit of search
                upper = p + 1 + 2 * p.Sqrt() + 3; // upper limit of search

                w = 1 + (upper - lower) / ordermod;
                leaps = w.Sqrt();

                mean = (5 * leaps) >> 1;
                nbits = (leaps >> 4).GetBits();
                if (nbits > 30) nbits = 30;
                distinguished = 1 << nbits;

                for (s = 1, m = 1; ; m++)
                { 
                    /* find table size */
                    distance[m - 1] = s * ordermod;
                    s <<= 1;

                    if ((2 * s / m) > mean)
                        break;
                }

                table[0] = ECMath.Multiply(curve, ordermod, P);

                for (i = 1; i < m; i++)
                { // double last entry
                    table[i] = table[i - 1];
                    table[i] = ECMath.Add(curve, table[i], table[i - 1]);
                }

                middle = (upper + lower) >> 1;

                if (ordermod > 1)
                    middle += (ordermod + order - middle % ordermod);

                for (i = 0; i < 5; i++) start[i] = middle + 13 * ordermod * i;
                for (i = 0; i < 5; i++) start[i + 5] = 13 * ordermod * i;

                ECPoint U = ECMath.Multiply(curve, middle, P);
                ECPoint V = ECMath.Multiply(curve, start[6], P);
                K[0] = U; K[5] = ECPoint.POINT_INFINITY;
                for (i = 0; i < 10; i++) D[i] = 0;

                for (i = 1; i < 5; i++)
                {
                    K[i] = ECMath.Add(curve, K[i - 1], V);
                    K[i + 5] = ECMath.Add(curve, K[i + 4], V);
                }

                nt = 0; nw = 0; cw = 0; ct = 0;
                collision = false; abort = false;

                for (; ; )
                {
                    for (jj = 0; jj < 5; jj++)
                    {
                        // random function...
                        txc = object.ReferenceEquals(K[jj].GetAffineX(), null) ? 0 : K[jj].GetAffineX();
                        i = (int)(txc % m);

                        if (txc % distinguished == 0)
                        {
                            if (nt >= 80)
                            {
                                abort = true;
                                break;
                            }

                            Console.Write(".");

                            tame[nt] = K[jj];
                            tdist[nt] = D[jj];
                            tname[nt] = jj;

                            for (k = 0; k < nw; k++)
                            {
                                if (wild[k] == tame[nt])
                                {
                                    ct = nt; cw = k;
                                    collision = true;
                                    break;
                                }
                            }

                            if (collision) break;
                            nt++;
                        }

                        D[jj] += distance[i];
                        TE[jj] = table[i];
                    }

                    if (collision || abort) break;

                    for (jj = 5; jj < 10; jj++)
                    {
                        // random function...
                        wxc = object.ReferenceEquals(K[jj].GetAffineX(), null) ? 0 : K[jj].GetAffineX();
                        j = (int)(wxc % m);

                        if (wxc % distinguished == 0)
                        {
                            if (nw >= 80)
                            {
                                abort = true;
                                break;
                            }

                            Console.Write(".");

                            wild[nw] = K[jj];
                            wdist[nw] = D[jj];
                            wname[nw] = jj;

                            for (k = 0; k < nt; k++)
                            {
                                if (tame[k] == wild[nw])
                                {
                                    ct = k; cw = nw;
                                    collision = true;
                                    break;
                                }
                            }

                            if (collision) break;
                            nw++;
                        }

                        D[jj] += distance[j];
                        TE[jj] = table[j];
                    }

                    if (collision || abort) break;

                    for (int l = 0; l < 10; l++)
                        K[l] = ECMath.Add(curve, K[l], TE[l]);
                }

                if (abort) continue;
                nrp = start[tname[ct]] - start[wname[cw]] + tdist[ct] - wdist[cw];
                G = P;

                G = ECMath.Multiply(curve, nrp, G);
                if (G != ZERO) continue;
                if (BigInteger.IsProbablePrime(nrp)) break;
                real_order = nrp; i = 0;
                Sieve sieve = new Sieve(256);

                for (; ; )
                {
                    if (i >= sieve.Count) break;
                    sp = sieve[i];

                    if (real_order % sp == 0)
                    {
                        G = P;
                        G = ECMath.Multiply(curve, real_order / sp, G);

                        if (G == ZERO)
                        {
                            real_order /= sp;
                            continue;
                        }
                    }

                    i++;
                }

                if (real_order <= 4 * p.Sqrt()) continue;
                real_order = nrp;

                for (i = 0; i < sieve.Count; i++)
                {
                    while (real_order % sp == 0)
                        real_order /= sp;
                }

                if (real_order == 1) break;

                if (BigInteger.IsProbablePrime(real_order))
                {
                    G = P;
                    G = ECMath.Multiply(curve, nrp / real_order, G);

                    if (G == ZERO) continue;
                    else break;
                }

                bad = false;

                for (i = 0; i < 20; i++)
                {
                    P = curve.BasePoint;
                    G = P;
                    G = ECMath.Multiply(curve, nrp, G);

                    if (G != ZERO)
                    {
                        bad = true;
                        break;
                    }
                }

                if (bad) continue;
                break;
            }

            return nrp;
        }

        static void ReadModularPolynomialsByLimit(BigInteger field)
        {
            BigInteger lim = 2;
            BigInteger Q4 = field.GetBits() > 256 ?
                (field >> 72).Sqrt() : (field >> 64).Sqrt();

            prime = new List<int>();
            Gl = new List<ModularPolynomial>();

            while (lim < Q4)
            {
                int p = int.Parse(sr.ReadLine());
                ModularPolynomial mp = new ModularPolynomial();

                while (true)
                {
                    BigInteger c = new BigInteger(sr.ReadLine());
                    int dx = int.Parse(sr.ReadLine());

                    int dy = int.Parse(sr.ReadLine());
                    mp.AddTerm(c, dx, dy);

                    if (dx == 0 && dy == 0)
                        break;
                }

                Gl.Add(mp);
                prime.Add(p);
            }
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
