#include "kinematics.hpp"
#include <cassert>
#include <cmath>

bool Point::operator==(Point const& rhs) const
{
  return (this->x_ == rhs.x_ && this->y_ == rhs.y_);
}

Result::Result(double x, double y, double theta)
    : p_{x, y}
    , theta_{theta}
{}
Result::Result(Point const& p, double theta)
    : p_{p}
    , theta_{theta}
{}

Result::Result()
    : p_{-1., 0.}
{}

Result Trajectory::result() const
{
  return Result(p_, std::atan(m_));
}

void Trajectory::exit(double x)
{
  p_.y_ = m_ * (x - p_.x_) + p_.y_;
  p_.x_ = x;
}

bool Result::operator==(Result b) const
{
  return (this->p_ == b.p_ && this->theta_ == b.theta_);
}

double Result::get_x() const
{
  return p_.x_;
}
double Result::get_y() const
{
  return p_.y_;
}
double Result::get_theta() const
{
  return theta_;
}

bool Result::is_valid() const
{
  return p_.x_ == -1.;
}

std::ostream& operator<<(std::ostream& os, Result const& res)
{
  os << res.get_x() << ", " << res.get_y() << ", " << res.get_theta();
  return os;
}

Barrier::Barrier(Pol pol, double x_max)
    : max_{x_max, pol(x_max)}
    , pol_{pol}
{}

double Barrier::max() const
{
  return max_.x_;
}

Pol Barrier::pol() const
{
  return pol_;
}

Point intersect(Trajectory t, Barrier b)
{
  Pol t_pol{std::vector<double>{t.p_.y_ - t.m_ * t.p_.x_, t.m_}};
  double x = eq_solve(t_pol, b.pol());
  return Point{x, t_pol(x)};
}

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Point p0,
                                double m0)
{
  Trajectory t{{0., p0.y_}, m0, false};
  double l  = barrier_up.max();
  double r1 = barrier_up.pol()(0.);
  double r2 = barrier_up.pol()(l);

  assert(std::abs(p0.y_) < r1);

  {
    double up_x{intersect(t, barrier_up).x_};     // up intersection x
    double down_x{intersect(t, barrier_down).x_}; // down intersection x

    if (up_x > 0 && down_x > 0)
      t.up_ = (up_x < down_x);
    else if (up_x > 0. && down_x < 0.)
      t.up_ = true;
  }
  Point inter;

  for (int i{0}; i < 100; ++i) {
    if (t.up_) {
      t.up_ = false;
      inter = intersect(t, barrier_up);

      if (inter.x_ < l && inter.x_ > 0) {
        t.p_ = inter;
        t.m_ = std::tan((2. * std::atan(-l / (r2 - r1))) - std::atan(t.m_));
      } else if (inter.x_ < 0) {
        return Result();
      } else if (inter.x_ > l) {
        t.exit(l);
        return t.result();
      }
    } else {
      t.up_ = true;
      inter = intersect(t, barrier_down);

      if (inter.x_ < l && inter.x_ > 0) {
        t.p_ = inter;
        t.m_ = std::tan(2. * std::atan(l / (r2 - r1)) - std::atan(t.m_));
      } else if (inter.x_ < 0) {
        return Result();
      } else if (inter.x_ > l) {
        t.exit(l);
        return t.result();
      }
    }
  }
  return t.result();
}