using Microsoft.EntityFrameworkCore;

namespace CryptoUtils.Data
{
    public class TrustGuardContext : DbContext
    {
        public TrustGuardContext(DbContextOptions<TrustGuardContext> options) : base(options)
        { }

        public DbSet<Domain> Domain { get; set; }
        public DbSet<Demand> Demand { get; set; }

        protected override void OnModelCreating(ModelBuilder modelBuilder)
        {
            modelBuilder.Entity<Domain>().ToTable(nameof(Domain));
            modelBuilder.Entity<Demand>().ToTable(nameof(Demand));
        }
    }
}
