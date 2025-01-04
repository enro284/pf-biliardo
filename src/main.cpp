#include "graphics.hpp"
#include "kinematics.hpp"
#include "statistics.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

template<typename T>
void set_from_user_input(T& var, const std::string& var_name)
{
  std::cout << "Enter " << var_name << ": ";
  std::cin >> var;
  if (!std::cin.good()) {
    throw(std::runtime_error("invalid input"));
  }
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int main()
{
  Barrier barrier_up;
  Barrier barrier_down;

  int n_sim{0};
  set_from_user_input(n_sim, "N of particles to simulate (1 with graphics, "
                             ">1 with just statistics)");

  int deg{0};
  set_from_user_input(deg, "barrier equation degree [1,2]");

  double l{0.};
  set_from_user_input(l, "length of the barrier (l)");

  switch (deg) {
  case 1: {
    double r1{0.};
    double r2{0.};
    set_from_user_input(r1, "height at beginning of the barrier (r1)");
    set_from_user_input(r2, "height at end of the barrier (r2)");

    barrier_up   = Barrier{l, r1, r2};
    barrier_down = Barrier{l, -r1, -r2};
    break;
  }
  case 2: {
    std::cout << "Upper barrier equation is: a * x^2 + b * x + c, in the range "
                 "[0,l]\n";
    double a{0.};
    double b{0.};
    double c{0.};
    set_from_user_input(a, "a");
    set_from_user_input(b, "b");
    set_from_user_input(c, "c");

    Pol p{{c, b, a}};
    barrier_up   = Barrier{p, l};
    barrier_down = Barrier{-p, l};
    break;
  }
  }

  if (n_sim == 1) { /* single particle simulation */
    double y0{0.};
    set_from_user_input(y0, "initial height (y0)");
    if (y0 >= barrier_up.pol()(0.) || y0 <= barrier_down.pol()(0.)) {
      throw(std::runtime_error("y0 out of bounds, cannot simulate trajectory"));
    }

    double theta0{0};
    set_from_user_input(theta0, "initial angle (theta0) [rad]");

    Trajectory traj{{0., y0}, theta0};

    std::vector<Vec2> bounces;
    Result res =
        simulate_single_particle(barrier_up, barrier_down, traj, &bounces);

    std::cout << "Simulation result (x, y, theta[rad]): " << res << '\n';

    Plot plot(barrier_up, barrier_down);
    plot.add(bounces);
    std::cout << "Opening the plot in a new window, press \"q\" to quit\n";
    plot.show();

  } else { /* N particle simulation */

    std::cout << "Enter the parameters of the two gaussian distributions:\n";
    double mu_y{0.};
    double sigma_y{1.};
    set_from_user_input(mu_y, "mu_y");
    set_from_user_input(sigma_y, "sigma_y");
    double mu_theta{0.};
    double sigma_theta{3.};
    set_from_user_input(mu_theta, "mu_theta");
    set_from_user_input(sigma_theta, "sigma_theta");

    std::random_device rd;
    std::default_random_engine eng{rd()};
    std::normal_distribution y_dist{mu_y, sigma_y};
    std::normal_distribution theta_dist{mu_theta, sigma_theta};

    Sample statistics_y;
    Sample statistics_theta;

    double r1{barrier_up.pol()(0.)};

    for (int i{0}; i != n_sim; ++i) {
      double y0{y_dist(eng)};
      while (std::abs(y0) > r1) {
        y0 = y_dist(eng);
      }

      double theta0{theta_dist(eng)};

      Trajectory traj{{0., y0}, theta0};
      Result res = simulate_single_particle(barrier_up, barrier_down, traj);

      if (res.get_x() == l) {
        double yf     = res.get_y();
        double thetaf = res.get_theta();
        statistics_y.add(yf);
        statistics_theta.add(thetaf);
      }
    }

    const auto stats_y     = statistics_y.statistics();
    const auto stats_theta = statistics_theta.statistics();

    std::cout
        << "The number of generated particles is " << n_sim
        << ", the number of particles exiting from the right side is "
        << statistics_theta.size()
        << '\n';

    std::cout << "The exit y values have a mean of " << stats_y.mean
              << ", a standard deviation of " << stats_y.std_dev
              << ", a skewnes coefficient of " << stats_y.skewness
              << " and a kurtosis of " << stats_y.kurtosis << '\n';

    std::cout << "The exit angle values have a mean of " << stats_theta.mean
              << ", a standard deviation of " << stats_theta.std_dev
              << ", a skewnes coefficient of " << stats_theta.skewness
              << " and a kurtosis of " << stats_theta.kurtosis << '\n';
  }
  return EXIT_SUCCESS;
}