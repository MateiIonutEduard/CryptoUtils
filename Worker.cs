using CryptoUtils.Data;
using CryptoUtils.Models;
using CryptoUtils.Services;
using Eduard;
using Microsoft.EntityFrameworkCore;
using System.Diagnostics;

namespace CryptoUtils
{
    public class Worker : BackgroundService
    {
        readonly IAppSettings appSettings;
        readonly IDomainService domainService;
        readonly ILogger<Worker> _logger;

        public Worker(IAppSettings appSettings, IDomainService domainService, ILogger<Worker> logger)
        {
            _logger = logger;
            this.appSettings = appSettings;
            this.domainService = domainService;
        }

        protected override async Task ExecuteAsync(CancellationToken stoppingToken)
        {
            while (!stoppingToken.IsCancellationRequested)
            {
                Demand[] demands = await domainService
                    .GetDemandsAsync();

                if (demands.Length == 0)
                    _logger.LogWarning("There are not exists demands! {time}", DateTime.UtcNow);
                else
                {
                    foreach(var demand in demands)
                    {
                        DateTime now = DateTime.UtcNow;
                        Stopwatch sw = new Stopwatch();
                        sw.Start();

                        _logger.LogInformation("Attempt to find strong elliptic curve domain parameters... {time}", now);
                        Core.SetDomain(appSettings.KeySize);
                        EllipticCurve curve = Core.GenParameters(appSettings.verbose);
                        sw.Stop();

                        double seconds = sw.Elapsed.TotalSeconds;
                        int demandId = demand.Id;
                        sw.Reset();

                        /* save elliptic curve parameters */
                        DateTime next = now.AddSeconds(seconds);
                        _logger.LogInformation("New elliptic curve domain parameters were generated at {time}", next);
                        await domainService.CreateDomainAsync(demandId, now, (int)seconds, curve);
                    }
                }

                await Task.Delay(appSettings.Delay, stoppingToken);
            }
        }
    }
}