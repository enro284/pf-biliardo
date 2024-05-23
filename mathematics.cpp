#include "mathematics.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
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
int Pol::deg() const
{
  return coeff_.size() - 1;
}

std::vector<double> Pol::coeff() const
{
  return coeff_;
}

Pol Pol::operator-()
{
  Pol res{*this};
  std::transform(res.coeff_.begin(), res.coeff_.end(), res.coeff_.begin(),
                 [](double c) { return -c; });
  return res;
}

double eq_solve(Pol const& pol1, Pol const& pol2)
{
  std::vector<double> eq;
  std::vector<double> minus;
  int eq_deg;
  if (pol2.deg() > pol1.deg()) {
    eq_deg = pol2.deg();
    eq     = pol2.coeff();
    minus  = pol1.coeff();
  } else {
    eq_deg = pol1.deg();
    eq     = pol1.coeff();
    minus  = pol2.coeff();
  }
  std::transform(minus.begin(), minus.end(), eq.begin(), eq.begin(),
                 [](double m, double e) { return e - m; });
  switch (eq_deg) {
  case 1:
    return -eq[0] / eq[1];
    break;
  default:
    throw std::runtime_error("equation degree too high, I cannot solve it");
    break;
  }
}