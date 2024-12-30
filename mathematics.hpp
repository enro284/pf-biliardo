#ifndef MATHEMATICS_HPP
#define MATHEMATICS_HPP

#include <vector>

class Pol
{
  std::vector<double> coeff_;

 public:
  Pol(std::vector<double> coeff);
  double operator()(double x);
  // double der(double x); //Operator to calculate derivative, needs to be
  // implemented for higher degree barriers

  int deg() const;
  std::vector<double> coeff() const;

  Pol operator-();
};

std::vector<double> eq_solve(Pol const& pol1, Pol const& pol2);

struct Vec2
{
  double x_;
  double y_;

  bool operator==(Vec2 const& rhs) const;
  Vec2 operator*(double rhs) const; // TODO: make free func
  double norm() const;
  Vec2 ortho() const;
};

Vec2 operator+(Vec2 const& lhs, Vec2 const& rhs);
double dot(Vec2 const& rhs, Vec2 const& lhs);

#endif