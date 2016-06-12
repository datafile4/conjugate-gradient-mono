using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConjugateGradient
{
    class Program
    {
        public struct interval
        {
            public double a;
            public double b;
        }
        const double eps = 0.00001;
        const int N = 2;

        static void Main(string[] args)
        {
            double[] X0 = new double[N];
            double[] result = new double[N];
            Console.Write("Enter X0: ");
            X0[0] = double.Parse(Console.ReadLine());
            Console.Write("Enter X1: ");
            X0[1] = double.Parse(Console.ReadLine());

            Console.Write("Enter tau: ");
            double tau = double.Parse(Console.ReadLine());
            /*X0[0] = 2;
			X0[1] = -2;*/
            ConjugateGradientMethod(X0, ref result, tau);
            //SimpleGradient(X0, ref result);
            //gradient(X0, ref result);

            Console.ReadKey();
        }

        public static double Sqr(double x)
        {
            return x * x;
        }

        static double f(double[] X)
        {

            double result = Sqr(Sqr(X[0]) - X[1]) + Sqr(1 - X[0]);
            return result;
        }

        static double f(double[] Xk, double[] dk, double t)
        {
            double[] X = new double[N];
            for (int i = 0; i < N; i++)
            {
                X[i] = Xk[i] + t * dk[i];
            }
            double result = f(X);
            return result;
        }

        static void gradient(double[] X, ref double[] res)
        {
            /*for (int i = 0; i < N; i++)
                Console.WriteLine(" in gr X {0} ", X[i]);
            Console.WriteLine();*/
            double h = 0.001;
            double result;
            X[0] = X[0] + h;
            double f1 = f(X);
            X[0] = X[0] - 2 * h;
            double f2 = f(X);
            result = 0.5 * (f1 - f2) / h;
            res[0] = result;
            X[0] = X[0] + h;

            X[1] = X[1] + h;
            f1 = f(X);
            X[1] = X[1] - 2 * h;
            f2 = f(X);
            result = 0.5 * (f1 - f2) / h;
            res[1] = result;
            X[1] = X[1] + h;
            /*
			res [0] = 2 * X [0];
			res [1] = 2 * X [1]; */
        }

        //одномерная оптимизация

        public static interval IntervalSearch(double x0, double h, double[] Xk, double[] direction)
        {
            double delta;
            if (f(Xk, direction, (x0 + h)) > f(Xk, direction, x0))
            {
                delta = -h;
                Console.WriteLine("-delta");
            }
            else
            {
                Console.WriteLine("+delta");
                delta = h;
            }

            double f1;
            double f2;
            double x1 = x0;
            double x2 = x0 + delta;
            do
            {
                f1 = f(Xk, direction, x1);
                f2 = f(Xk, direction, x2);
                x2 = x2 + delta;
            } while (f1 > f2);

            interval result;
            if (h > 0)
            {
                result.a = x1;
                result.b = x2;
            }
            else
            {
                result.a = x2;
                result.b = x1;
            }
            return result;
        }
        static interval Swan(double x0, double h, double[] Xk, double[] direction)
        {
            double a, b;
            double xl = x0 - h;
            double xr = x0 + h;
            double fl = f(Xk, direction, xl);
            double fc = f(Xk, direction, x0);
            double fr = f(Xk, direction, xr);
            interval segment;
            segment.a = 0;
            segment.b = 0;

            if (fl >= fc && fc <= fr)
            {
                a = xl;
                b = xr;
                segment.a = a;
                segment.b = b;
                return segment;
            }

            if (fl <= fc && fc >= fr)
            {
                segment.a = 0;
                segment.b = 0;
                return segment;
            }

            if (fl >= fc && fc >= fr)
            {
                double pow = 2;
                a = x0;
                xl = xr;
                fl = fr;
                xr = x0 + pow * h;
                fr = f(Xk, direction, xr);

                while (fl >= fr)
                {
                    pow *= 2;
                    if (pow > 2454542)
                    {
                        segment.a = 0;
                        segment.b = 0;
                        return segment;
                    }
                    a = xl;
                    xl = xr;
                    fl = fr;
                    xr = x0 + pow * h;
                    fr = f(Xk, direction, xr);
                }
                b = xr;
                segment.a = a;
                segment.b = b;
                return segment;
            }
            if (fl <= fc && fc <= fr)
            {
                double pow = 2;
                b = x0;
                xr = xl;
                fr = fl;
                xl = x0 - pow * h;
                fl = f(Xk, direction, xl);
                while (fl <= fr)
                {
                    pow *= 2;
                    if (pow > 2452458)
                    {
                        segment.a = 0;
                        segment.b = 0;
                        return segment;
                    }
                    b = xr;
                    xr = xl;
                    fr = fl;
                    xl = x0 - pow * h;
                    fl = f(Xk, direction, xl);
                }
                a = xl;
                segment.a = a;
                segment.b = b;
                return segment;
            }
            return segment;
        }



        static double GoldenMethod(double a, double b, double[] Xk, double[] direction)
        {
            double x = a + (3 - Math.Sqrt(5)) * (b - a) / 2;
            double y = a + (Math.Sqrt(5) - 1) * (b - a) / 2;
            double fx = f(Xk, direction, x);
            double fy = f(Xk, direction, y);

            while (b - a > eps)
            {
                if (fx <= fy)
                {
                    b = y;
                    y = x;
                    fy = fx;
                    x = a + b - y;
                    fx = f(Xk, direction, x);
                }
                else
                {
                    a = x;
                    x = y;
                    fx = fy;
                    y = a + b - x;
                    fy = f(Xk, direction, y);
                }

            }
            double point = (b + a) / 2;
            return point;
        }


        //a - left, b - center, c - right
        public static double GoldenSearch(double a, double b, double c, double tau, double[] xNow, double[] grDirection)
        {
            double x;
            double phi = (1 + Math.Sqrt(5)) / 2;
            double resphi = 2 - phi;
            if (b < c)
                x = b + resphi * (c - b);
            else
                x = b - resphi * (b - a);
            if (Math.Abs(c - a) < tau * (Math.Abs(b) + Math.Abs(x)))
                return (c + a) / 2;
            if (f(xNow, grDirection, x) < f(xNow, grDirection, b))
                return (b < c) ? GoldenSearch(b, x, c, tau, xNow, grDirection) : GoldenSearch(a, x, b, tau, xNow, grDirection);
            else
                return (b < c) ? GoldenSearch(a, b, x, tau, xNow, grDirection) : GoldenSearch(x, b, c, tau, xNow, grDirection);
        }


        static double skalar(double[] X, double[] Y)
        {
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += X[i] * Y[i];
            }
            return sum;
        }

        static double norm(double[] X)
        {
            double sum = 0;
            for (int i = 0; i < N; i++)
                sum += X[i] * X[i];
            double result = Math.Sqrt(sum);
            return result;
        }

        //Polak-Ribiere method
        static double polak(double[] fk, double[] fkPrev)
        {
            double[] tmp = new double[N];
            for (int i = 0; i < N; i++)
            {
                tmp[i] = fk[i] - fkPrev[i];
            }
            double result = Math.Sqrt(skalar(fk, tmp)) / Math.Sqrt(skalar(fk, fk));
            return result;
        }

        //Flatcher-Rievse method
        static double fletcher(double[] gradNow, double[] gradPrev)
        {
            //double result = Math.Sqrt(skalar(gradNow, gradNow)) / Math.Sqrt(skalar(gradPrev, gradPrev));
            double result = skalar(gradNow, gradNow) / skalar(gradPrev, gradPrev);
            return result;
        }

        static void direction(double[] gradNow, double[] dirPrev, double bk, ref double[] dirNew)
        {
            for (int i = 0; i < N; i++)
                dirNew[i] = -gradNow[i] + bk * dirPrev[i];
        }

        static void SimpleGradient(double[] X0, ref double[] res)
        {
            double[] gradFNow = new double[N];
            double[] result = new double[N];
            double t = 0.1;
            double[] gradFPrev = new double[N];
            double[] xNow = new double[N];
            double[] xPrev = new double[N];
            double fPrev;
            double fNow;

            for (int i = 0; i < N; i++)
            {
                xPrev[i] = X0[i];
            }

            do
            {
                gradient(xPrev, ref gradFNow);
                for (int i = 0; i < N; i++)
                {
                    xNow[i] = xPrev[i] - t * gradFNow[i];
                    Console.Write(" xNow= {0}", xNow[i]);
                }
                Console.WriteLine();

                fPrev = f(xPrev);
                fNow = f(xNow);
                Console.WriteLine("fprev {0}", fPrev);
                Console.WriteLine("fnow {0} ", fNow);

                if (fNow > fPrev)
                    t /= 2;
                for (int i = 0; i < N; i++)
                    xPrev[i] = xNow[i];
            } while (Math.Abs(fNow - fPrev) > eps);

            for (int i = 0; i < N; i++)
                res[i] = xNow[i];
        }

        public static void ConjugateGradientMethod(double[] X0, ref double[] res, double intervalStep)
        {
            double[] destinationVect = new double[N];
            double skalarWeight;
            double[] xNow = new double[N];
            double[] xPrev = new double[N];
            double[] gradNow = new double[N];
            double[] gradPrev = new double[N];
            double t = 0.0;
            interval inter = new interval();
            for (int i = 0; i < N; i++)
            {
                X0[i] = xNow[i];
            }

            do
            {
                //for(int l = 0; l<2000;l++){
                for (int i = 0; i < N; i++)
                {
                    destinationVect[i] = -xNow[i];
                }

                for (int k = 0; k <= N; k++)
                {
                    inter = Swan(0, intervalStep, xNow, destinationVect);
                    t = GoldenMethod(inter.a, inter.b, xNow, destinationVect);
                    Console.WriteLine("t = {0}", t);

                    for (int i = 0; i < N; i++)
                    {
                        xPrev[i] = xNow[i];
                        xNow[i] = xPrev[i] + t * destinationVect[i];
                        Console.Write("X[{0}] = {1} ", i, xNow[i]);
                    }

                    Console.WriteLine();
                    gradient(xNow, ref gradNow);
                    gradient(xPrev, ref gradPrev);
                    skalarWeight = fletcher(gradNow, gradPrev);

                    for (int i = 0; i < N; i++)
                    {
                        destinationVect[i] = -gradNow[i] + skalarWeight * destinationVect[i];
                        Console.Write("dest[{0}] = {1} ", i, destinationVect[i]);
                    }
                    Console.WriteLine();
                }

            } while (Math.Abs(f(xNow) - f(xPrev)) / f(xNow) > eps);
        }

        static void ConjugateGradient(double[] X0, ref double[] res, double tau)
        {
            double[] gradFNow = new double[N];
            double[] result = new double[N];
            double[] xNow = new double[N];
            double[] xPrev = new double[N];
            double t = 0;
            double[] grDirection = new double[N];
            double[] grDirectionPrev = new double[N];
            double[] gradFPrev = new double[N];
            double grBk;
            double count = 0;
            int k = 0;
            for (int i = 0; i < N; i++)
                xNow[i] = X0[i];

            //for(int s = 0; s< 1000; s++)
            do
            {
                gradient(xNow, ref gradFNow);
                //******************************
                if (k == N || k == 0)
                {
                    for (int i = 0; i < N; i++)
                    {
                        grDirection[i] = -gradFNow[i];
                    }
                    k = 0;
                }
                else
                {
                    for (int i = 0; i < N; i++)
                        grDirectionPrev[i] = grDirection[i];
                    gradient(xPrev, ref gradFPrev);
                    grBk = fletcher(gradFNow, gradFPrev);
                    direction(gradFNow, grDirectionPrev, grBk, ref grDirection);
                }
                //****************************************


                interval segment = Swan(0, tau, xNow, grDirection);

                double ser = (segment.b + segment.a) / 2;
                t = GoldenMethod(segment.a, segment.b, xNow, grDirection);
                Console.WriteLine("a = {0}, b = {1}, t = {2}", segment.a, segment.b, t);
                for (int i = 0; i < N; i++)
                    xPrev[i] = xNow[i];
                for (int i = 0; i < N; i++)
                {
                    xNow[i] = xPrev[i] + t * grDirection[i];
                    Console.Write("xNow {0} ", xNow[i]);
                }
                Console.WriteLine();
                Console.WriteLine("F = {0} ", Math.Round(f(xNow), 6, MidpointRounding.AwayFromZero));
                k++;
                count++;

                //проверка на возрастание
                /*if(f(xNow) > f(xPrev)){
					Console.WriteLine("Функция возрастает");
					break;
				}*/
            } //while (Math.Abs(f(xNow) -f(xPrev))>=eps);
            while (Math.Abs(f(xNow) - f(xPrev)) / f(xNow) > eps);
            Console.WriteLine(count);
            for (int i = 0; i < N; i++)
                res[i] = xNow[i];
        }
    }
}
