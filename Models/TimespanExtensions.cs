using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CryptoUtils.Models
{
    public static class TimespanExtensions
    {
        /* Prints the timespan object in custom mode. */
        public static string ToCustomString(this TimeSpan span)
        {
            return string.Format("{0:00}:{1:00}:{2:00}", span.Hours, span.Minutes, span.Seconds);
        }
    }
}
