using CryptoUtils.Data;
using CryptoUtils.Models;
using Microsoft.EntityFrameworkCore;
using Microsoft.Extensions.Options;

namespace CryptoUtils
{
    public class Program
    {
        public static void Main(string[] args)
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
                    services.AddHostedService<Worker>();
                })
                .Build();

            host.Run();
        }
    }
}