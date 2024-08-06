#include "statistics.hpp"
#include <cmath>
#include <stdexcept>
Sample::Sample()
    : n{0}
    , s1{0}
    , s2{0}
    , s3{0}
    , s4{0}
{}

void Sample::add(double x)
{
  ++n;
  s1 += x;
  s2 += x * x;
  s3 += std::pow(x, 3);
  s4 += std::pow(x, 4);
}
int Sample::size()
{
  return n;
}

Statistics Sample::statistics() const
{
  if (n < 4) {
    throw std::runtime_error{"Not enough entries to run a statistics"};
  }

  // https://mathworld.wolfram.com/SampleCentralMoment.html
  double M1 = s1 / n;
  double M2 = s2 / n - M1 * M1;
  double M3 = s3 / n - 3 * M1 * s3 + 2 * std::pow(M1, 3);
  double M4 =
      s4 / n - 4 * M1 * s3 / n + 6 * M1 * M1 * s2 / n + 3 * std::pow(M1, 4);

  double mean    = M1;
  double std_dev = std::sqrt(M2 * n / (n - 1));

  // http://brownmath.com/stat/shape.htm#SkewnessCompute
  //  adj fisher-pearson
  double skew = std::sqrt(static_cast<double>(n) * static_cast<double>(n - 1))
              / static_cast<double>(n - 2) * M3 / std::pow(M2, 1.5);
  // sample excess kurtosis: Joanes, D. N., and C. A. Gill. 1998.
  // “Comparing Measures of Sample Skewnessand Kurtosis”. The
  // Statistician 47(1): 183–189.
  double kurt = static_cast<double>(n - 1)
              / (static_cast<double>(n - 2) * static_cast<double>(n - 3))
              * (static_cast<double>(n + 1) * M4 / M2 * M2 + 6.);

  return {mean, std_dev, skew, kurt};
}