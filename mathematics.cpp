#include "mathematics.hpp"
#include <cmath>
#include <numeric>
#include <vector>

Pol::Pol(std::vector<double> coeff)
    : coeff_(coeff)
{}
double Pol::operator()(double x)
{
  int i{-1};
  return std::accumulate(coeff_.begin(), coeff_.end(), 0.,
                         [x, &i](double acc, double coeff) {
                           ++i;
                           return acc += coeff * std::pow(x, i);
                         });
};
