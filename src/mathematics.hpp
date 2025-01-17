#ifndef MATHEMATICS_HPP
#define MATHEMATICS_HPP

#include <vector>

class Pol
{
  // [0]x^0 + [1]x^1 + [2]x^2 ...
  std::vector<double> coeff_;

 public:
  Pol(std::vector<double> coeff);

  double operator()(double x) const;
  double der(double x) const;
  std::size_t deg() const;
  std::vector<double> coeff() const;

  Pol operator-();
};

std::vector<double> eq_solve(Pol const& pol1, Pol const& pol2);

struct Vec2
{
  double x_;
  double y_;

  bool operator==(Vec2 const& rhs) const;
  Vec2& operator*=(double rhs);
  Vec2& operator+=(Vec2 rhs);
  Vec2& operator-=(Vec2 rhs);

  double norm() const;
  Vec2 ortho() const;

  double dist2(Vec2 const& v) const;
};

Vec2 operator*(double rhs, Vec2 const& lhs);
Vec2 operator*(Vec2 const& rhs, double lhs);
Vec2 operator/(Vec2 const& rhs, double lhs);

Vec2 operator+(Vec2 const& lhs, Vec2 const& rhs);
Vec2 operator-(Vec2 const& lhs, Vec2 const& rhs);

double dot(Vec2 const& rhs, Vec2 const& lhs);

#endif