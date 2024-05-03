using CryptoUtils.Data;
using CryptoUtils.Models;
using CryptoUtils.Services;
using Microsoft.EntityFrameworkCore;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Options;

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
                    Console.WriteLine("usage: BotanicTool [-bits] [size_in_bits] [-verb] [-S]\n");
                    Console.WriteLine("Represents a tool that generates the elliptic curve parameters in finite fields for cryptographic use.");
                    Console.WriteLine("This implementation favors strong elliptic curve parameter generation in finite fields, based on \nmodified SEA algorithm.\n");

                    Console.WriteLine("options:");
                    Console.WriteLine("-h, --help\tshow this help message and exit");
                    Console.WriteLine("-bits, --bits\tSize in bits of the prime field of elliptic curve");
                    Console.WriteLine("-verb, --verbose  Flag that enable the console messages printing when compute cardinality of elliptic curve");
                    Console.WriteLine("-S, --search\tThis flag activates searching of the ideal strong cryptographic elliptic curves");
                }
                else
                {

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