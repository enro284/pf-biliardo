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
}

int Pol::deg() const
{
  return static_cast<int>(coeff_.size()) - 1;
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

std::vector<double> eq_solve(Pol const& pol1, Pol const& pol2)
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

  minus.resize(eq_deg, 0.0);

  std::transform(minus.begin(), minus.end(), eq.begin(), eq.begin(),
                 [](double m, double e) { return e - m; });

  std::vector<double> sol;
  switch (eq_deg) {
  case 1: { // ax + b = 0
    double a = eq[1];
    double b = eq[0];
    if (std::abs(a) < 1e-10) { // TODO: EPS
      break;
    }
    sol.push_back({-b / a});
  }
  case 2: { // ax^2 + bx + c = 0
    double a = eq[2];
    double b = eq[1];
    double c = eq[0];

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
      break;
    } else if (std::abs(discriminant) < 1e-10) {
      sol.push_back({-b / (2 * a)});
    } else {
      double sqrt_disc = std::sqrt(discriminant);
      sol.push_back((-b - sqrt_disc) / (2 * a));
      sol.push_back((-b + sqrt_disc) / (2 * a));
    }
  }
  default:
    throw std::runtime_error(
        "Equation degree too high, only 2nd degree is supported");
  }
  return sol;
}

bool Vec2::operator==(Vec2 const& rhs) const
{
  return (this->x_ == rhs.x_ && this->y_ == rhs.y_);
}

double dot(Vec2 const& rhs, Vec2 const& lhs)
{
  return rhs.x_ * lhs.x_ + rhs.y_ + rhs.y_;
}

double Vec2::norm() const
{
  return std::sqrt(x_ * x_ + y_ * y_);
}

Vec2 Vec2::ortho() const
{
  return {-y_, x_};
}

Vec2 operator+(Vec2 const& lhs, Vec2 const& rhs)
{
  return {lhs.x_ + rhs.x_, lhs.y_ + rhs.y_};
}

Vec2 Vec2::operator*(double rhs) const
{
  return {x_ * rhs, y_ * rhs};
}
