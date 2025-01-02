#include "kinematics.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

double eps = 1e-10;

Result::Result(Vec2 const& p, double theta)
    : p_{p}
    , theta_{theta}
{}

Result::Result()
    : p_{-1., 0.}
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
  if (t.v_.x_ == 0) { // TODO: may need eps: y_ calculated from t
    if (t.p_.y_ == b->pol()(t.p_.x_)) {
      std::cout << "same sol\n";
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
      // TODO capture type, equal?
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
// t has to be a copy
{
  assert(std::abs(t.p_.y_) < barrier_up.pol()(0.));

  double l = barrier_up.max(); // TODO: necessary?

  std::vector<Bounce> bounces;
  bounces.reserve(10); // TODO n barriers * barrier order

  for (int i{0}; i < 50; ++i) {
    points.push_back(t.p_);
    bounces.clear();

    auto up_int   = intersect(t, &barrier_up);
    auto down_int = intersect(t, &barrier_down);
    bounces.insert(bounces.end(), up_int.begin(), up_int.end());
    bounces.insert(bounces.end(), down_int.begin(), down_int.end());

    if (bounces.size() != 0) { /* choose bounce and update trajectory */
      std::sort(bounces.begin(), bounces.end(),
                [&](Bounce const& lhs, Bounce const& rhs) {
                  return lhs.p_.dist2(t.p_) < rhs.p_.dist2(t.p_);
                });
      auto bounce = bounces[0];
      std::cout << i << "\tbounce (x,y) " << bounce.p_.x_ << ", "
                << bounce.p_.y_ << '\n';
      std::cout << "\tup size: " << up_int.size()
                << ", down size: " << down_int.size() << '\n';

      for (auto pippo : bounces) {
        std::cout << pippo.p_.x_ << " , " << pippo.p_.y_ << '\n';
      }
      t.p_ = bounce.p_;

      Vec2 tg = Vec2{1., bounce.b_ptr->pol().der(bounce.p_.x_)};
      tg      = tg * (1. / tg.norm());
      std::cout << "  tg: " << tg.x_ << ", " << tg.y_ << '\n';
      Vec2 n = tg.ortho();
      std::cout << "  n: " << n.x_ << ", " << n.y_ << '\n';

      auto v_new = n * (-dot(t.v_, n)) + tg * dot(t.v_, tg);
      t.v_       = v_new;

      std::cout << "  new v: " << t.v_.x_ << ", " << t.v_.y_ << '\n';
    } else { /* exit right or left */
      if (t.v_.x_ > 0) {
        t.exit(l);
        points.push_back(t.p_);
        return t.result(); // TODO: change with break
      } else if (t.v_.x_ < 0) {
        t.exit(0.);
        points.push_back(t.p_);
        return Result();
      } // TODO: =0
    }
  }

  return t.result();
}

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Trajectory t)
// t has to be a copy
{
  assert(std::abs(t.p_.y_) < barrier_up.pol()(0.));

  double l = barrier_up.max(); // TODO: necessary?

  std::vector<Bounce> bounces;
  bounces.reserve(10); // TODO n barriers * barrier order

  for (int i{0}; i < 50; ++i) {
    bounces.clear();

    auto up_int   = intersect(t, &barrier_up);
    auto down_int = intersect(t, &barrier_down);
    bounces.insert(bounces.end(), up_int.begin(), up_int.end());
    bounces.insert(bounces.end(), down_int.begin(), down_int.end());

    if (bounces.size() != 0) { /* choose bounce and update trajectory */
      std::sort(bounces.begin(), bounces.end(),
                [&](Bounce const& lhs, Bounce const& rhs) {
                  return lhs.p_.dist2(t.p_) < rhs.p_.dist2(t.p_);
                });
      auto bounce = bounces[0];
      std::cout << i << "\tbounce (x,y) " << bounce.p_.x_ << ", "
                << bounce.p_.y_ << '\n';
      std::cout << "\tup size: " << up_int.size()
                << ", down size: " << down_int.size() << '\n';

      for (auto pippo : bounces) {
        std::cout << pippo.p_.x_ << " , " << pippo.p_.y_ << '\n';
      }
      t.p_ = bounce.p_;

      Vec2 tg = Vec2{1., bounce.b_ptr->pol().der(bounce.p_.x_)};
      tg      = tg * (1. / tg.norm());
      std::cout << "  tg: " << tg.x_ << ", " << tg.y_ << '\n';
      Vec2 n = tg.ortho();
      std::cout << "  n: " << n.x_ << ", " << n.y_ << '\n';

      auto v_new = n * (-dot(t.v_, n)) + tg * dot(t.v_, tg);
      t.v_       = v_new;

      std::cout << "  new v: " << t.v_.x_ << ", " << t.v_.y_ << '\n';
    } else { /* exit right or left */
      if (t.v_.x_ > 0) {
        t.exit(l);
        return t.result(); // TODO: change with break
      } else if (t.v_.x_ < 0) {
        t.exit(0.);
        return Result();
      } // TODO: =0
    }
  }

  return t.result();
}