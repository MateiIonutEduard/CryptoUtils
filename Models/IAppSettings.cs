using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CryptoUtils.Models
{
    public interface IAppSettings
    {
        int Delay { get; set; }
        int KeySize { get; set; }
        bool verbose { get; set; }
    }
}
