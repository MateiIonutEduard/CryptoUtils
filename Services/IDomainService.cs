using Eduard;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CryptoUtils.Services
{
    public interface IDomainService
    {
        Task CreateDomainAsync(int demandId, DateTime solvedAt, int seconds, EllipticCurve curve);
    }
}
