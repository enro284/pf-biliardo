#include "kinematics.hpp"
#include <cassert>
#include <iostream>

int main()
{
  double l{4};
  double r1{1.5};
  double r2{0.7};

  Pol pol{std::vector<double>{r1, (r2-r1)/l}};
  Barrier barrier_up{pol, l};
  Barrier barrier_down{-pol, l};

  while (true) {
    double y0{0};
    double m0{0};
    std::cout << "enter y0: ";
    std::cin >> y0;
    assert(y0 <= r1 && y0 >= -r1);
    std::cout << "enter m0: ";
    std::cin >> m0;

    Result res =
        simulate_single_particle(barrier_up, barrier_down, Point{0, y0}, m0);

    std::cout << "result (x, y, theta[rad]): " << res << '\n';
  }
}
