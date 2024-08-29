#include "kinematics.hpp"
#include "statistics.hpp"
#include <random>

template<typename T>
void set_from_user_input(T& var, std::string var_name)
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
  std::cout << "Simulation of N particles. r1 is set as follows:\n";
  double r1{1.5};
  std::cout << "barriers: r1 = " << r1 << '\n';
  double l{0.};
  set_from_user_input(l, "length of barrier");
  double r2{0.};
  set_from_user_input(r2, "height at end of barrier");


  int N{0};
  set_from_user_input(N, "number of particles to simulate");

  std::cout << "Enter the parameters of the two gaussian distributions:\n";
  double mu_y{0.};
  set_from_user_input(mu_y, "mu_y");
  double sigma_y{1.};
  set_from_user_input(sigma_y, "sigma_y");
  double mu_theta{0.};
  set_from_user_input(mu_theta, "mu_theta");
  double sigma_theta{3.};
  set_from_user_input(sigma_theta, "sigma_theta");

  Barrier barrier_up{l, r1, r2};
  Barrier barrier_down{l, -r1, -r2};

  std::random_device rd;
  std::default_random_engine eng{rd()};
  std::normal_distribution y_dist{mu_y, sigma_y};
  std::normal_distribution theta_dist{mu_theta, sigma_theta};

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
    if (res.get_x() != -1) {
      double yf     = res.get_y();
      double thetaf = res.get_theta();
      statistics_y.add(yf);
      statistics_theta.add(thetaf);
    }
  }

  const auto stats_y     = statistics_y.statistics();
  const auto stats_theta = statistics_theta.statistics();

  std::cout << "The number of generated particles is " << N
            << ", the number of particles exiting from the right side is "
            << statistics_theta.size()
            << '\n'; // statistics_theta.size() is equal to statistics_y.size()

  std::cout << "The exit y values have a mean of " << stats_y.mean
            << ", a standard deviation of " << stats_y.std_dev
            << ", a skewnes coefficient of " << stats_y.skewness
            << " and a kurtosis of " << stats_y.kurtosis << '\n';

  std::cout << "The exit angle values have a mean of " << stats_theta.mean
            << ", a standard deviation of " << stats_theta.std_dev
            << ", a skewnes coefficient of " << stats_theta.skewness
            << " and a kurtosis of " << stats_theta.kurtosis << '\n';
}