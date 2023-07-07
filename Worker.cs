using CryptoUtils.Data;
using CryptoUtils.Models;
using CryptoUtils.Services;
using Eduard;
using Microsoft.EntityFrameworkCore;

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
                _logger.LogInformation("Worker running at: {time}", DateTimeOffset.Now);
                await Task.Delay(appSettings.Delay, stoppingToken);
            }
        }
    }
}