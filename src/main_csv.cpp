#include "kinematics.hpp"
#include <fstream>
#include <random>

std::string filename{"out.csv"};
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
  double r1{0.};
  double r2{0.};
  double l{0.};
  set_from_user_input(r1, "height at beginning of the barrier (r1)");
  set_from_user_input(r2, "height at end of the barrier (r2)");
  set_from_user_input(l, "length of the barrier (l)");

  int N{0};
  set_from_user_input(N, "number of particles to simulate");

  std::cout << "Enter the parameters of the two gaussian distributions:\n";
  double mu_y{0.};
  double sigma_y{1.};
  set_from_user_input(mu_y, "mu_y");
  set_from_user_input(sigma_y, "sigma_y");
  double mu_theta{0.};
  double sigma_theta{3.};
  set_from_user_input(mu_theta, "mu_theta");
  set_from_user_input(sigma_theta, "sigma_theta");

  Barrier barrier_up{l, r1, r2};
  Barrier barrier_down{l, -r1, -r2};

  std::random_device rd;
  std::default_random_engine eng{rd()};
  std::normal_distribution y_dist{mu_y, sigma_y};
  std::normal_distribution theta_dist{mu_theta, sigma_theta};

  std::ofstream csv_file;
  csv_file.open(filename);

  for (int i{0}; i != N; ++i) {
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
      csv_file << yf << ", " << thetaf << '\n';
    }
  }
  csv_file.close();
  std::cout << "Output written to \"" << filename << "\"\n";
  std::cout << "Columns are: yf, thetaf\n";
}