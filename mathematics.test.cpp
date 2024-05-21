#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "mathematics.hpp"
#include "doctest.h"

TEST_CASE("testint pol")
{
  std::vector<double> coeff0{0.};
  SUBCASE("calling pol with no parameters (null) in x")
  {
    Pol pol{coeff0};
    CHECK(pol(0.) == doctest::Approx(0.));
    CHECK(pol(1.) == doctest::Approx(0.));
  }

  std::vector<double> coeff{2.1};
  SUBCASE("calling pol with 1 parameters (null) in x")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(2.1));
    CHECK(pol(2.1) == doctest::Approx(2.1));
  }

  coeff.push_back(1.5);
  SUBCASE("calling pol with 2 parameters (linear) in x")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(3.6));
    CHECK(pol(2.1) == doctest::Approx(5.25));
  }

  std::vector<double> coeff3{2.1, 0.7, 3.4};
  SUBCASE("calling pol with 3 parameters (quadratic) in x")
  {
    Pol pol{coeff3};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(6.2));
    CHECK(pol(2.1) == doctest::Approx(18.564));
  }
}