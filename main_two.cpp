#include "kinematics.hpp"
#include "statistics.hpp"
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

  Sample statistics_y;
  Sample statistics_theta;

  for (int i{0}; i != N; ++i) {
    double y0{y_dist(eng)};
    while (std::abs(y0) > r1) {
      y0 = y_dist(eng);
    }
    double m0{std::tan(theta_dist(eng))};
    Result res =
        simulate_single_particle(barrier_up, barrier_down, Point{0, y0}, m0);
    double yf     = res.get_y();
    double thetaf = res.get_theta();
    statistics_y.add(yf);
    statistics_theta.add(thetaf);
  }

  const auto result_y     = statistics_y.statistics();
  const auto result_theta = statistics_theta.statistics();

  std::cout << "The exit y values have a mean of " << result_y.mean
            << " , a standard deviation of " << result_y.std_dev
            << ", a skewnes coefficient of " << result_y.skewness
            << " and a kurtosis of " << result_y.kurtosis << '\n';

  std::cout << "The exit angle values have a mean of " << result_theta.mean
            << " , a standard deviation of " << result_theta.std_dev
            << ", a skewnes coefficient of " << result_theta.skewness
            << " and a kurtosis of " << result_theta.kurtosis << '\n';
}