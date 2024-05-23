#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "kinematics.hpp"
#include "doctest.h"

TEST_CASE("testing single particle simulation")
{
  Pol barrier_pol{std::vector<double>{1.5, -0.2}};
  double l{4};

  Barrier barrier_up{barrier_pol, l};
  Barrier barrier_down{-barrier_pol, l};
/*
  SUBCASE("0 bouces")
  {
    CHECK(simulate_single_particle(barrier_up, barrier_down, Point{0., 0.}, 0.)
          == Result{Point{l, 0.}, 0.});
  }
*/
  SUBCASE("1 bounces")
  {
    CHECK(simulate_single_particle(barrier_up, barrier_down, Point{0., 0.}, 0.3)
          == Result{Point{l, 0.0809523809524}, -0.686247914178});
  }
  SUBCASE("2 bounces")
  {
    CHECK(simulate_single_particle(barrier_up, barrier_down, Point{0., 0.}, 0.5)
          == Result{Point{l, 0.0932203389752}, 1.2532298484});
  }
  SUBCASE("many bounces, exit left")
  {
    CHECK(simulate_single_particle(barrier_up, barrier_down, Point{0., 0.}, -1.)
          == Result{Point{-1., 0.}, 0.});
  }
}