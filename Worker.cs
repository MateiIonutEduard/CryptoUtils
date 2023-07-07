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
                    _logger.LogInformation("There are not exists demands! {time}", DateTime.UtcNow);
                else
                {
                    foreach(var demand in demands)
                    {
                        DateTime now = DateTime.UtcNow;
                        Stopwatch sw = new Stopwatch();
                        sw.Start();

                        Core.SetDomain(appSettings.KeySize);
                        EllipticCurve curve = Core.GenParameters(true);
                        sw.Stop();

                        int seconds = (int)sw.Elapsed.TotalSeconds;
                        int demandId = demand.Id;

                        /* save elliptic curve parameters */
                        await domainService.CreateDomainAsync(demandId, now, seconds, curve);
                    }
                }

                await Task.Delay(appSettings.Delay, stoppingToken);
            }
        }
    }
}