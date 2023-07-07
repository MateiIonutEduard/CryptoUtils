using System;
using Eduard;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CryptoUtils.Data;
using Microsoft.EntityFrameworkCore;

namespace CryptoUtils.Services
{
    public class DomainService : IDomainService
    {
        readonly TrustGuardContext guardContext;

        public DomainService(IServiceScopeFactory scopeFactory)
        {
            guardContext = scopeFactory.CreateScope().ServiceProvider.GetRequiredService<TrustGuardContext>();
        }

        public async Task<Demand[]> GetDemandsAsync()
        {
            Demand[] demands = await guardContext.Demand
                .Where(e => (e.IsSeen == null || (e.IsSeen != null && !e.IsSeen.Value)))
                .ToArrayAsync();

            /* non solved demands */
            return demands;
        }

        public async Task<bool> SolveAsync(int demandId, DateTime solvedAt, int seconds)
        {
            Demand? demand = await guardContext.Demand
                .FirstOrDefaultAsync(e => e.Id == demandId && (e.IsSeen == null || (e.IsSeen != null && !e.IsSeen.Value)));

            /* demand exists, but was resolved */
            if (demand != null)
            {
                demand.SolvedAt = solvedAt;
                demand.TotalSeconds = seconds;

                await guardContext.SaveChangesAsync();
                return true;
            }


            return false;
        }

        public async Task CreateDomainAsync(int demandId, DateTime solvedAt, int seconds, EllipticCurve curve)
        {
            Guid guid = Guid.NewGuid();
            ECPoint basePoint = curve.BasePoint;

            Demand? demand = await guardContext.Demand
                .FirstOrDefaultAsync(e => e.Id == demandId);

            /* if not null, update demand */
            if(demand != null)
            {
                demand.SolvedAt = solvedAt;
                demand.TotalSeconds = seconds;

                demand.IsSeen = true;
                await guardContext.SaveChangesAsync();
            }

            /* creates new elliptic curve domain parameters */
            Domain domain = new Domain
            {
                a = curve.a.ToString(),
                b = curve.b.ToString(),
                p = curve.field.ToString(),
                N = curve.order.ToString(),
                x = basePoint.GetAffineX().ToString(),
                y = basePoint.GetAffineY().ToString(),
                webcode = guid.ToString()
            };

            /* save into database */
            guardContext.Domain.Add(domain);
            await guardContext.SaveChangesAsync();
        }
    }
}
