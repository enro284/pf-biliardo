#include "kinematics.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

int main()
{
  double l{4};
  double r1{1.5};
  double r2{0.7};

  Pol pol{std::vector<double>{r1, (r2 - r1) / l}};
  Barrier barrier_up{pol, l};
  Barrier barrier_down{-pol, l};

  while (true) {
    double y0{0};
    std::cout << "enter y0: ";
    std::cin >> y0;
    if (!std::cin.good()) {
      throw(std::runtime_error("invalid input"));
    }
    if (std::abs(y0) > r1) {
      throw(std::runtime_error("y0 > r1, cannot simulate trajectory"));
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    double theta0{0};
    std::cout << "enter theta0[rad]: ";
    std::cin >> theta0;
    if (!std::cin.good()) {
      throw(std::runtime_error("invalid input"));
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    double m0{std::tan(theta0)};
    Result res =
        simulate_single_particle(barrier_up, barrier_down, Point{0, y0}, m0);

    std::cout << "result (x, y, theta[rad]): " << res << '\n';
    std::cout << "Enter 'q' to quit or any other key to continue: ";
    std::string input;
    std::cin >> input;
    if (input == "q")
      return 1;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
}
