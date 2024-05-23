#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "mathematics.hpp"
#include "doctest.h"

TEST_CASE("testint pol")
{
  std::vector<double> coeff0{0.};
  SUBCASE("calling pol with 1 parameter (null) in x")
  {
    Pol pol{coeff0};
    CHECK(pol(0.) == doctest::Approx(0.));
    CHECK(pol(1.) == doctest::Approx(0.));
  }

  std::vector<double> coeff{2.1};
  SUBCASE("calling pol with 1 parameter in x")
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

TEST_CASE("testing eq_solve")
{
  std::vector<double> coeff_bar{
      3., -1.}; // first coeff is ==r1, second coeff is slope of barrier

  Pol bar{coeff_bar};

  SUBCASE("calling eq_solve with: linear barrier, flat trajectory")
  {
    std::vector<double> coeff_tra{0.};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(3.));
  }

  SUBCASE("calling eq_solve with: linear barrier, defined trajectory")
  {
    std::vector<double> coeff_tra{1., 2.3};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(0.606060));
  }

  SUBCASE("calling eq_solve with: linear barrier, defined trajectory which "
          "eq_solves negatively")
  {
    std::vector<double> coeff_tra{1., -2.3};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(-1.53846154));
  }

  SUBCASE("calling eq_solve with: linear barrier, defined trajectory which "
          "eq_solves out of bounds")
  {
    std::vector<double> coeff_tra{0., -0.2};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(3.75));
  }

  SUBCASE("calling eq_solve with: linear barrier, parallel trajectory")
  {
    std::vector<double> coeff_tra{2.99, -1.};
    Pol tra{coeff_tra};
    CHECK(std::isinf(eq_solve(tra, bar)));
  }

  SUBCASE("calling eq_solve with: linear barrier, equal trajectory")
  {
    std::vector<double> coeff_tra{3., -1.};
    Pol tra{coeff_tra};
    CHECK(std::isnan(eq_solve(tra, bar)));
  }
}

TEST_CASE("testing Pol::operator-")
{
  std::vector<double> coeff{2.99, -1.};

  Pol pol{coeff};

  Pol inverse_pol{-pol};

  CHECK(inverse_pol.coeff()[0] == -2.99);
  CHECK(inverse_pol.coeff()[1] == 1.);
}
