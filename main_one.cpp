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
  Barrier barrier_up{1, 0, 0}; // TODO: default
  Barrier barrier_down{1, 0, 0};

  int deg{0};
  set_from_user_input(deg, "barrier equation degree [1,2]");

  switch (deg) {
  case 1: {
    double r1{0.};
    double r2{0.};
    double l{0.};
    set_from_user_input(r1, "height at beginning of the barrier (r1)");
    set_from_user_input(r2, "height at end of the barrier (r2)");
    set_from_user_input(l, "length of the barrier (l)");

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
    double l{0.};
    set_from_user_input(a, "a");
    set_from_user_input(b, "b");
    set_from_user_input(c, "c");
    set_from_user_input(l, "length of the barrier (l)");

    Pol p{{c, b, a}};
    barrier_up = Barrier{p, l};
    barrier_down = Barrier{-p, l};
    break;
  }
  }

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
      simulate_single_particle(barrier_up, barrier_down, traj, bounces);

  std::cout << "Simulation result (x, y, theta[rad]): " << res << '\n';

  Plot plot(barrier_up, barrier_down);
  plot.add(bounces);
  plot.show();

  return EXIT_SUCCESS;
}