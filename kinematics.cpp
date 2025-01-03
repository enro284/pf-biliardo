#include "kinematics.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

double eps = 1e-10;

Result::Result(Vec2 const& p, double theta)
    : p_{p}
    , theta_{theta}
{}

Result Trajectory::result() const
{
  return Result(p_, std::atan(v_.y_ / v_.x_));
}

void Trajectory::exit(double x)
{
  double time = (x - p_.x_) / v_.x_;

  p_.x_ = x;
  p_.y_ += time * v_.y_;
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

Trajectory::Trajectory(Vec2 const& p, double theta)
    : p_{p}
    , v_{cos(theta), sin(theta)}
{}

Barrier::Barrier(Pol const& pol, double x_max)
    : max_{x_max, pol(x_max)}
    , pol_{pol}
{
  assert(x_max > 0.);
}

Barrier::Barrier(double l, double r1, double r2)
    : max_{l, r2}
    , pol_{std::vector<double>{r1, (r2 - r1) / l}}
{
  assert(l > 0.);
}

double Barrier::max() const
{
  return max_.x_;
}

Pol Barrier::pol() const
{
  return pol_;
}

std::vector<Bounce> intersect(Trajectory const& t, Barrier const* b)
{
  if (t.v_.x_ < eps && t.v_.x_ > -eps) {
    if (t.p_.y_ == b->pol()(t.p_.x_)) {
      return {};
    } else {
      return {{{t.p_.x_, b->pol()(t.p_.x_)}, b}};
    }
  } else {
    double t_m = t.v_.y_ / t.v_.x_;
    Pol t_pol{{t.p_.y_ - t_m * t.p_.x_, t_m}};

    std::vector<double> sol_x = eq_solve(t_pol, b->pol());

    if (t.v_.x_ > 0) {
      std::erase_if(
          sol_x, [&](double x) { return x - t.p_.x_ <= eps || x > b->max(); });
    } else if (t.v_.x_ < 0) {
      std::erase_if(sol_x,
                    [&](double x) { return x - t.p_.x_ >= -eps || x <= 0.; });
    }
    std::vector<Bounce> sol(sol_x.size(), Bounce{});
    std::transform(sol_x.begin(), sol_x.end(), sol.begin(),
                   [=](double x) { return Bounce{{x, t_pol(x)}, b}; });
    return sol;
  }
}

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Trajectory t,
                                std::vector<Vec2>& points)
{
  assert(std::abs(t.p_.y_) < barrier_up.pol()(0.));

  std::vector<Bounce> bounces;
  bounces.reserve(4); // TODO: n barriers * barrier order

  for (int i{0}; i < 50; ++i) {
    points.push_back(t.p_);
    bounces.clear();

    auto up_int   = intersect(t, &barrier_up);
    auto down_int = intersect(t, &barrier_down);
    // TODO: make bounces global or pass them to avoid copying too many times
    bounces.insert(bounces.end(), up_int.begin(), up_int.end());
    bounces.insert(bounces.end(), down_int.begin(), down_int.end());

    if (bounces.size() != 0) { /* choose bounce and update trajectory */

      Bounce bounce =
          *std::min_element(bounces.begin(), bounces.end(),
                            [&](Bounce const& lhs, Bounce const& rhs) {
                              return lhs.p_.dist2(t.p_) < rhs.p_.dist2(t.p_);
                            });

      t.p_ = bounce.p_;

      Vec2 tg = Vec2{1., bounce.b_ptr->pol().der(bounce.p_.x_)};
      tg      = tg * (1. / tg.norm());
      Vec2 n  = tg.ortho();

      auto v_new = n * (-dot(t.v_, n)) + tg * dot(t.v_, tg);
      t.v_       = v_new;

    } else { /* exit right or left */
      if (t.v_.x_ > 0) {
        t.exit(barrier_up.max());
      } else {
        t.exit(0.);
      }
      points.push_back(t.p_);
      return t.result();
    }
  }

  return t.result();
}
