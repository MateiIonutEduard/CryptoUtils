using Eduard;
using CryptoUtils.Data;

namespace CryptoUtils.Services
{
    public interface IDomainService
    {
        Task<Demand[]> GetDemandsAsync();
        Task<bool> SolveAsync(int demandId, DateTime solvedAt, int seconds);
        Task CreateDomainAsync(int demandId, DateTime solvedAt, int seconds, EllipticCurve curve);
    }
}
