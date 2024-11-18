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
  if (res.get_x() == -1) {
    std::cout << "Ball goes backwards and exits system from origin";
    return os;
  } else {
    os << res.get_x() << ", " << res.get_y() << ", " << res.get_theta();
    return os;
  }
}

Barrier::Barrier(Pol pol, double x_max)
    : max_{x_max, pol(x_max)}
    , pol_{pol}
{}

Barrier::Barrier(double l, double r1, double r2)
    : max_{l, r2}
    , pol_{std::vector<double>{r1, (r2 - r1) / l}}
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
  Trajectory t{{0., p0.y_}, m0};
  double l  = barrier_up.max();
  double r1 = barrier_up.pol()(0.);
  double r2 = barrier_up.pol()(l);

  assert(std::abs(p0.y_) < r1);

  Point intersection;
  {
    bool up;
    double up_x{intersect(t, barrier_up).x_};
    double down_x{intersect(t, barrier_down).x_};

    if (up_x > 0 && down_x > 0)
      up = (up_x < down_x);
    else if (up_x > 0. && down_x < 0.)
      up = true;
    else
      up = false;

    if (!up) {
      intersection = intersect(t, barrier_down);

      if (intersection.x_ < l && intersection.x_ > 0) {
        t.p_ = intersection;
        t.m_ = std::tan(2. * std::atan(l / (r2 - r1)) - std::atan(t.m_));
      } else if (intersection.x_ < 0) {
        return Result();
      } else if (intersection.x_ > l) {
        t.exit(l);
        return t.result();
      }
    }
  }

  for (int i{0}; i < 50; ++i) {
    intersection = intersect(t, barrier_up);
    if (intersection.x_ < l && intersection.x_ > 0) {
      t.p_ = intersection;
      t.m_ = std::tan((2. * std::atan(-l / (r2 - r1))) - std::atan(t.m_));
    } else if (intersection.x_ < 0) {
      return Result();
    } else if (intersection.x_ > l) {
      t.exit(l);
      return t.result();
    }

    intersection = intersect(t, barrier_down);
    if (intersection.x_ < l && intersection.x_ > 0) {
      t.p_ = intersection;
      t.m_ = std::tan(2. * std::atan(l / (r2 - r1)) - std::atan(t.m_));
    } else if (intersection.x_ < 0) {
      return Result();
    } else if (intersection.x_ > l) {
      t.exit(l);
      return t.result();
    }
  }

  return t.result();
}

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Point p0,
                                double m0, std::vector<Point> &points)
{
  points.push_back(p0);

  Trajectory t{{0., p0.y_}, m0};
  double l  = barrier_up.max();
  double r1 = barrier_up.pol()(0.);
  double r2 = barrier_up.pol()(l);

  assert(std::abs(p0.y_) < r1);

  Point intersection;
  {
    bool up;
    double up_x{intersect(t, barrier_up).x_};
    double down_x{intersect(t, barrier_down).x_};

    if (up_x > 0 && down_x > 0)
      up = (up_x < down_x);
    else if (up_x > 0. && down_x < 0.)
      up = true;
    else
      up = false;

    if (!up) {
      intersection = intersect(t, barrier_down);
      points.push_back(intersection);
      if (intersection.x_ < l && intersection.x_ > 0) {
        t.p_ = intersection;
        t.m_ = std::tan(2. * std::atan(l / (r2 - r1)) - std::atan(t.m_));
      } else if (intersection.x_ < 0) {
        return Result();
      } else if (intersection.x_ > l) {
        t.exit(l);
        return t.result();
      }
    }
  }

  for (int i{0}; i < 50; ++i) {
    intersection = intersect(t, barrier_up);
    points.push_back(intersection);
    if (intersection.x_ < l && intersection.x_ > 0) {
      t.p_ = intersection;
      t.m_ = std::tan((2. * std::atan(-l / (r2 - r1))) - std::atan(t.m_));
    } else if (intersection.x_ < 0) {
      return Result();
    } else if (intersection.x_ > l) {
      t.exit(l);
      return t.result();
    }

    intersection = intersect(t, barrier_down);
    points.push_back(intersection);
    if (intersection.x_ < l && intersection.x_ > 0) {
      t.p_ = intersection;
      t.m_ = std::tan(2. * std::atan(l / (r2 - r1)) - std::atan(t.m_));
    } else if (intersection.x_ < 0) {
      return Result();
    } else if (intersection.x_ > l) {
      t.exit(l);
      return t.result();
    }
  }

  return t.result();
}