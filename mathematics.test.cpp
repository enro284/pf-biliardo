#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "mathematics.hpp"
#include "doctest.h"

TEST_CASE("testint pol")
{
  SUBCASE("empty polynomial")
  {
    CHECK_THROWS(Pol(std::vector<double>()));
  }

  SUBCASE("0 degree polynomial")
  {
    std::vector<double> coeff{0.};
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(0.));
    CHECK(pol(1.) == doctest::Approx(0.));

    CHECK(pol.der(0.) == 0.);
    CHECK(pol.der(1.) == 0.);
  }

  std::vector<double> coeff{2.1};
  SUBCASE("1st degree polynomial")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(2.1));
    CHECK(pol(2.1) == doctest::Approx(2.1));

    CHECK(pol.der(0.) == 0.);
    CHECK(pol.der(1.) == 0.);
  }

  coeff.push_back(1.5);
  SUBCASE("2nd degree polynomial")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(3.6));
    CHECK(pol(2.1) == doctest::Approx(5.25));

    CHECK(pol.der(0.) == 1.5);
    CHECK(pol.der(1.) == 1.5);
    CHECK(pol.der(2.1) == 1.5);
  }

  coeff.push_back(3.4);
  SUBCASE("3rd degree polynomial")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(7.));
    CHECK(pol(2.1) == doctest::Approx(20.244));

    CHECK(pol.der(0.) == doctest::Approx(1.5));
    CHECK(pol.der(1.) == doctest::Approx(8.3));
    CHECK(pol.der(2.1) == doctest::Approx(15.78));
  }

  coeff.push_back(5.7);
  SUBCASE("4th degree polynomial")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(12.7));
    CHECK(pol(1.8) == doctest::Approx(49.0584));

    CHECK(pol.der(0.) == doctest::Approx(1.5));
    CHECK(pol.der(1.) == doctest::Approx(25.4));
    CHECK(pol.der(1.8) == doctest::Approx(69.144));
  }
}

TEST_CASE("testing eq_solve")
{
  std::vector<double> coeff_bar{3., -1.};

  Pol bar{coeff_bar};

  SUBCASE("eq_solve with: linear barrier, flat trajectory")
  {
    std::vector<double> coeff_tra{0.};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(3.));
  }

  SUBCASE("eq_solve with: linear barrier, defined trajectory")
  {
    std::vector<double> coeff_tra{1., 2.3};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(0.606060));
  }

  SUBCASE(" eq_solve with: linear barrier, defined trajectory which "
          "eq_solves negatively")
  {
    std::vector<double> coeff_tra{1., -2.3};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(-1.53846154));
  }

  SUBCASE(" eq_solve with: linear barrier, defined trajectory which "
          "eq_solves out of bounds")
  {
    std::vector<double> coeff_tra{0., -0.2};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar) == doctest::Approx(3.75));
  }

  SUBCASE(" eq_solve with: linear barrier, parallel trajectory")
  {
    std::vector<double> coeff_tra{2.99, -1.};
    Pol tra{coeff_tra};
    CHECK(std::isinf(eq_solve(tra, bar)));
  }

  SUBCASE(" eq_solve with: linear barrier, equal trajectory")
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
