#include "kinematics.hpp"
#include "globals.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

Result::Result(Vec2 const& p, double theta)
    : p_{p}
    , theta_{theta}
{}

Result Trajectory::result() const
{
  return Result(p_, std::atan2(v_.y_, v_.x_));
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

std::vector<Collision> intersect(Trajectory const& t, Barrier const* b)
{
  if (std::abs(t.v_.x_) < Globals::EPS) { /* handle vertical trajectory */
    if (t.p_.y_ == b->pol()(t.p_.x_)) {
      return {};
    } else {
      return {{{t.p_.x_, b->pol()(t.p_.x_)}, b}};
    }
  }

  double t_m = t.v_.y_ / t.v_.x_;
  Pol t_pol{{t.p_.y_ - t_m * t.p_.x_, t_m}};

  std::vector<double> sol_x = eq_solve(t_pol, b->pol());

  // filter solutions based on particle direction
  if (t.v_.x_ > 0) {
    std::erase_if(sol_x, [&](double x) {
      return x - t.p_.x_ <= Globals::EPS || x > b->max();
    });
  } else if (t.v_.x_ < 0) {
    std::erase_if(sol_x, [&](double x) {
      return x - t.p_.x_ >= -Globals::EPS || x <= 0.;
    });
  }
  std::vector<Collision> sol(sol_x.size(), Collision{});
  std::transform(sol_x.begin(), sol_x.end(), sol.begin(),
                 [&](double x) { return Collision{{x, t_pol(x)}, b}; });
  return sol;
}

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Trajectory t,
                                std::vector<Vec2>* bounces)
{
  assert(std::abs(t.p_.y_) < barrier_up.pol()(0.));

  std::vector<Collision> collisions;
  collisions.reserve(4);

  for (int i{0}; i < Globals::MAX_ITERATIONS; ++i) {
    if (bounces)
      bounces->push_back(t.p_);

    collisions.clear();

    auto up_int   = intersect(t, &barrier_up);
    auto down_int = intersect(t, &barrier_down);
    collisions.insert(collisions.end(), up_int.begin(), up_int.end());
    collisions.insert(collisions.end(), down_int.begin(), down_int.end());

    if (collisions.size() != 0) { /* choose bounce and update trajectory */

      Collision bounce =
          *std::min_element(collisions.begin(), collisions.end(),
                            [&](Collision const& lhs, Collision const& rhs) {
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
      if (bounces)
        bounces->push_back(t.p_);

      return t.result();
    }
  }

  return t.result();
}
