#include <cmath>
#include <numeric>
#include <vector>

class Function
{
 public:
};

class Pol : Function
{
  std::vector<double> coeff_;

 public:
  double operator()(double x)
  {
    int i{-1};
    return std::accumulate(
        coeff_.begin(), coeff_.end(), 0,
        [x, &i](double acc, double coeff) { 
            ++i;
            return acc += coeff * std::pow(x, i); });
  };
  double der(double x_);
};

double intersect(Function t, Function b);