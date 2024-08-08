#include "kinematics.hpp"
#include <random>
int main()
{
  int N{100};

  double l{4};
  double r1{1.5};
  double r2{0.7};

  double u_y{0.};
  double sigma_y{1.};
  double u_theta{0.};
  double sigma_theta{3.};

  Pol pol{std::vector<double>{r1, (r2 - r1) / l}};
  Barrier barrier_up{pol, l};
  Barrier barrier_down{-pol, l};

  std::default_random_engine eng;
  std::normal_distribution y_dist{u_y, sigma_y};
  std::normal_distribution theta_dist{u_theta, sigma_theta};

  for (int i{0}; i != N; ++i) {
    double y0{y_dist(eng)};
    double m0{std::tan(theta_dist(eng))};
    simulate_single_particle(barrier_up, barrier_down, Point{0, y0}, m0);
  }
}