#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "kinematics.hpp"
#include "doctest.h"

TEST_CASE("testing single particle simulation")
{
  Pol barrier_pol{std::vector<double>{1.5, -0.2}};
  double l{4};

  Barrier barrier_up{barrier_pol, l};
  Barrier barrier_down{-barrier_pol, l};

  SUBCASE("0 bouces")
  {
    Result res =
        simulate_single_particle(barrier_up, barrier_down, {{0., 0.}, 0.});
    CHECK(res.get_x() == doctest::Approx(4.));
    CHECK(res.get_y() == doctest::Approx(0.));
    CHECK(res.get_theta() == doctest::Approx(-0.));
  }

  SUBCASE("1 bounces")
  {
    Result res = simulate_single_particle(barrier_up, barrier_down,
                                          {{0., 0.}, 0.291456794});
    CHECK(res.get_x() == doctest::Approx(4.));
    CHECK(res.get_y() == doctest::Approx(0.0809523809524));
    CHECK(res.get_theta() == doctest::Approx(-0.686247914178));
  }
  SUBCASE("2 bounces")
  {
    Result res = simulate_single_particle(barrier_up, barrier_down,
                                          {{0., 0.}, 0.463647609});
    CHECK(res.get_x() == doctest::Approx(4.));
    CHECK(res.get_y() == doctest::Approx(0.0932203389752));
    CHECK(res.get_theta() == doctest::Approx(1.2532298484));
  }
  SUBCASE("many bounces, exit left")
  {
    CHECK(simulate_single_particle(barrier_up, barrier_down,
                                   {{0., 0.}, -0.785398163})
              .get_x()
          == 0);
  }
}

TEST_CASE("testing single particle simulation with parallel barriers")
{
  Pol barrier_pol{std::vector<double>{1.5, 0.}};
  double l{4};

  Barrier barrier_up{barrier_pol, l};
  Barrier barrier_down{-barrier_pol, l};

  SUBCASE("0 bounces")
  {
    Result res = simulate_single_particle(barrier_up, barrier_down,
                                          {{0., 0.}, 0.291456794});
    CHECK(res.get_x() == doctest::Approx(4.));
    CHECK(res.get_y() == doctest::Approx(1.2));
    CHECK(res.get_theta() == doctest::Approx(0.291456794478));
  }

  SUBCASE("1 bounce")
  {
    Result res = simulate_single_particle(barrier_up, barrier_down,
                                          {{0., 0.}, 0.463647609});
    CHECK(res.get_x() == doctest::Approx(4.));
    CHECK(res.get_y() == doctest::Approx(1.));
    CHECK(res.get_theta() == doctest::Approx(-0.463647609001));
  }

  SUBCASE("particle parallel to both barriers")
  {
    Result res =
        simulate_single_particle(barrier_up, barrier_down, {{0., 0.}, 0.});
    CHECK(res.get_x() == doctest::Approx(4.));
    CHECK(res.get_y() == doctest::Approx(0.));
    CHECK(res.get_theta() == doctest::Approx(-0.));
  }

  SUBCASE("testing too many bounces")
  {
    Result res = simulate_single_particle(barrier_up, barrier_down,
                                          {{0., 0.}, 1.55829698});
    CHECK(res.get_x() < 4);
  }
}