#include "graphics.hpp"
#include "kinematics.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

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
  double r1{0.};
  set_from_user_input(r1, "height at beginning of the barrier (r1)");
  double r2{0.};
  set_from_user_input(r2, "height at end of the barrier (r2)");
  double l{0.};
  set_from_user_input(l, "length of the barrier (l)");

  double y0{0.};
  set_from_user_input(y0, "initial height (y0)");
  if (std::abs(y0) >= r1) {
    throw(std::runtime_error("y0 >= r1, cannot simulate trajectory"));
  }

  double theta0{0};
  set_from_user_input(theta0, "initial angle (theta0) [rad]");
  double m0{std::tan(theta0)};

  Barrier barrier_up{l, r1, r2};
  Barrier barrier_down{l, -r1, -r2};

  std::vector<Point> bounces;
  Result res = simulate_single_particle(barrier_up, barrier_down, Point{0, y0},
                                        m0, bounces);

  std::cout << "Simulation result (x, y, theta[rad]): " << res << '\n';

  Plot plot(barrier_up, barrier_down);
  plot.add(bounces);
  plot.show();

  return EXIT_SUCCESS;
}
