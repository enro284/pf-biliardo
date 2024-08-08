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

  // https://mathworld.wolfram.com/SampleCentralmoment.html
  double m1 = s1 / n;
  double m2 = s2 / n - m1 * m1;
  double m3 = s3 / n - 3 * m1 * s2 / n + 2 * std::pow(m1, 3);
  double m4 =
      s4 / n - 4 * m1 * s3 / n + 6 * m1 * m1 * s2 / n - 3 * std::pow(m1, 4);
  
  double N = static_cast<double>(n);

  double mean    = m1;

  double std_dev = std::sqrt(m2 * n / (n - 1));

  // http://brownmath.com/stat/shape.htm#SkewnessCompute
  double skew = std::sqrt(N * (N-1) )
              / (N - 2) * m3 / std::pow(m2, 1.5);

  double kurt = (N - 1)
              / ((N - 2) * (N - 3))
              * ((N + 1) * (m4 / (m2 * m2) - 3.) + 6.);

  return {mean, std_dev, skew, kurt};
}