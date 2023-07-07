using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CryptoUtils.Models
{
    public class AppSettings : IAppSettings
    {
        public int Delay { get; set; }
        public int KeySize { get; set; }
    }
}
