#include <iostream>

#include "kinematics.hpp"

int main()
{
  Pol pol{std::vector<double>{1.5, -0.2}};
  double l{4};
  Barrier barrier_up{pol, l};
  Barrier barrier_down{-pol, l};

  while (true) {
    double y0{0};
    double m{0};
    std::cout << "enter y0: ";
    std::cin >> y0;
    std::cout << "enter m: ";
    std::cin >> m;

    Result res =
        simulate_single_particle(barrier_up, barrier_down, Point{0, y0}, m);

    std::cout << "result (x, y, theta): " << res << '\n';
  }
}
