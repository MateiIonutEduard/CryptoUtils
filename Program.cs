using CryptoUtils.Data;
using Microsoft.EntityFrameworkCore;

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
                    services.AddHostedService<Worker>();
                })
                .Build();

            host.Run();
        }
    }
}