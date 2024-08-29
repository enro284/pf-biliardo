#include "kinematics.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

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
  double r1{1.5};
  double l{0.};
  set_from_user_input(l, "length of barrier");
  double r2{0.};
  set_from_user_input(r2, "height at end of barrier");

  std::cout << "barriers: r1 = " << r1 << '\n';

  Barrier barrier_up{l, r1, r2};
  Barrier barrier_down{l, -r1, -r2};

  while (true) {
    double y0{0};
    set_from_user_input(y0, "y0");
    if (std::abs(y0) > r1) {
      throw(std::runtime_error("y0 > r1, cannot simulate trajectory"));
    }

    double theta0{0};
    set_from_user_input(theta0, "theta0");
    double m0{std::tan(theta0)};

    Result res =
        simulate_single_particle(barrier_up, barrier_down, Point{0, y0}, m0);

    std::cout << "result (x, y, theta[rad]): " << res << '\n';

    std::cout << "Type 'q' to quit or any other key to continue: ";
    std::string input;
    std::cin >> input;
    if (input == "q")
      return EXIT_SUCCESS;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
}
