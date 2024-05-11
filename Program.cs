using CryptoUtils.Data;
using CryptoUtils.Models;
using CryptoUtils.Services;
using Eduard;
using Microsoft.EntityFrameworkCore;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Options;
using System.Diagnostics;

namespace CryptoUtils
{
    public class Program
    {
        public static void Main(string[] args)
        {
            /* treat this application as a console tool */
            if (args.Length > 0)
            {
                /* help requested about this app */
                bool requestHelp = !string.IsNullOrEmpty(args.FirstOrDefault(arg => 
                    arg.StartsWith("-h") || arg.StartsWith("--help")));

                // if it is requested help, show how to use this app
                if (requestHelp)
                {
                    Console.WriteLine("usage: CryptoUtils [-bits] [size_in_bits] [-verb] [-S]\n");
                    Console.WriteLine("Represents a tool that generates the elliptic curve parameters in finite fields for cryptographic use.");
                    Console.WriteLine("This implementation favors strong elliptic curve parameter generation in finite fields, based on \nmodified SEA algorithm.\n");

                    Console.WriteLine("All options:");
                    Console.WriteLine("-h, --help\tshow this help message and exit");
                    Console.WriteLine("-bits, --bits\tSize in bits of the prime field of elliptic curve");
                    Console.WriteLine("-verb, --verbose  Flag that enable the console messages printing when compute cardinality of elliptic curve");
                    Console.WriteLine("-S, --search\tThis flag activates searching of the ideal strong cryptographic elliptic curves");
                    Environment.Exit(1);
                }
                else
                {
                    bool verbose = !string.IsNullOrEmpty(args.FirstOrDefault(arg =>
                        arg.StartsWith("-verb") || arg.StartsWith("--verbose")));

                    bool search = !string.IsNullOrEmpty(args.FirstOrDefault(arg =>
                        arg.StartsWith("-S") || arg.StartsWith("--search")));

                    int index = Array.IndexOf(args, "-bits");
                    if(index == -1) index = Array.IndexOf(args, "--bits");

                    int? bits = index >= 0 ? Convert.ToInt32(args[index + 1]) : null;
					Stopwatch sw = new Stopwatch();

					if (bits != null)
                    {
                        Core.SetDomain(bits.Value, search);
                        sw.Start();
						EllipticCurve curve = Core.GenParameters(verbose);
                        sw.Stop();

                        Console.WriteLine($"a = {curve.a}");
                        Console.WriteLine($"b = {curve.b}");

                        Console.WriteLine($"field = {curve.field}");
                        Console.WriteLine($"order = {curve.order}");

						Console.WriteLine($"Total Time: {sw.Elapsed.ToCustomString()}");
						Environment.Exit(0);
                    }
                    else
                    {
                        BigInteger a = new BigInteger(args[0]);
						BigInteger b = new BigInteger(args[1]);
						BigInteger field = new BigInteger(args[2]);

						Core.SetDomain(a, b, field, search);
                        sw.Start();
						EllipticCurve curve = Core.GenParameters(verbose);
                        sw.Stop();

						Console.WriteLine($"a = {curve.a}");
						Console.WriteLine($"b = {curve.b}");

						Console.WriteLine($"field = {curve.field}");
						Console.WriteLine($"order = {curve.order}");

                        Console.WriteLine($"Total Time: {sw.Elapsed.ToCustomString()}");
						Environment.Exit(0);
					}
                }
            }
            else
            {
                IHost host = Host.CreateDefaultBuilder(args)
                    .ConfigureServices((hostContext, services) =>
                    {
                        /* Add PostgreSQL database */
                        services.AddEntityFrameworkNpgsql().AddDbContext<TrustGuardContext>(opt =>
                            opt.UseNpgsql(hostContext.Configuration.GetConnectionString("TrustGuard")));

                        /* get AppSettings config section */
                        services.Configure<AppSettings>(
                            hostContext.Configuration.GetSection(nameof(AppSettings)));

                        /* register AppSettings as singleton service */
                        services.AddSingleton<IAppSettings>(sp =>
                            sp.GetRequiredService<IOptions<AppSettings>>().Value);

                        /* declares domain parameters service */
                        services.AddTransient<IDomainService, DomainService>();
                        services.AddHostedService<Worker>();
                    })
                    .Build();

                host.Run();
            }
        }
    }
}