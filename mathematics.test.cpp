#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "mathematics.hpp"
#include "doctest.h"

TEST_CASE("testing pol")
{
  std::vector<double> coeff0{0.};
  SUBCASE("calculate pol with 1 parameter (null) in x")
  {
    Pol pol{coeff0};
    CHECK(pol(0.) == doctest::Approx(0.));
    CHECK(pol(1.) == doctest::Approx(0.));
  }

  std::vector<double> coeff{2.1};
  SUBCASE("calculate pol with 1 parameter in x")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(2.1));
    CHECK(pol(2.1) == doctest::Approx(2.1));
  }

  coeff.push_back(1.5);
  SUBCASE("calculate pol with 2 parameters (linear) in x")
  {
    Pol pol{coeff};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(3.6));
    CHECK(pol(2.1) == doctest::Approx(5.25));
  }

  std::vector<double> coeff3{2.1, 0.7, 3.4};
  SUBCASE("calculate pol with 3 parameters (quadratic) in x")
  {
    Pol pol{coeff3};
    CHECK(pol(0.) == doctest::Approx(2.1));
    CHECK(pol(1.) == doctest::Approx(6.2));
    CHECK(pol(2.1) == doctest::Approx(18.564));
  }
}

TEST_CASE("testing eq_solve")
{
  Pol bar{{3., -1.}};

  SUBCASE("calling eq_solve with: linear barrier, flat trajectory")
  {
    Pol tra{{0.}};
    CHECK(eq_solve(tra, bar)[0] == doctest::Approx(3.));
  }

  SUBCASE("calling eq_solve with: linear barrier, defined trajectory")
  {
    std::vector<double> coeff_tra{1., 2.3};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar)[0] == doctest::Approx(0.606060));
  }

  SUBCASE("calling eq_solve with: linear barrier, defined trajectory which "
          "eq_solves negatively")
  {
    std::vector<double> coeff_tra{1., -2.3};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar)[0] == doctest::Approx(-1.53846154));
  }

  SUBCASE("calling eq_solve with: linear barrier, defined trajectory which "
          "eq_solves out of bounds")
  {
    std::vector<double> coeff_tra{0., -0.2};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar)[0] == doctest::Approx(3.75));
  }

  SUBCASE("calling eq_solve with: linear barrier, parallel trajectory")
  {
    std::vector<double> coeff_tra{2.99, -1.};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar).size() == 0);
  }

  SUBCASE("calling eq_solve with: linear barrier, equal trajectory")
  {
    std::vector<double> coeff_tra{3., -1.};
    Pol tra{coeff_tra};
    CHECK(eq_solve(tra, bar).size() == 0);
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

TEST_CASE("Testing Vec2") {
    Vec2 v1{3.0, 4.0};
    Vec2 v2{1.0, 2.0};

    SUBCASE("Equality operator") {
        Vec2 v3{3.0, 4.0};
        CHECK(v1 == v3);
        CHECK(!(v1 == v2));
    }

    SUBCASE("Scalar multiplication") {
        Vec2 scaled = v1 * 2.0;
        CHECK(scaled == Vec2{6.0, 8.0});
    }

    SUBCASE("Norm computation") {
        CHECK(v1.norm() == doctest::Approx(5.0));
    }

    SUBCASE("Orthogonal vector") {
        Vec2 ortho = v1.ortho();
        CHECK(dot(v1, ortho) == doctest::Approx(0.0));
    }

    SUBCASE("Squared distance computation") {
        CHECK(v1.dist2(v2) == doctest::Approx(8.0));
    }

    SUBCASE("Vector addition and subtraction") {
        Vec2 sum = v1 + v2;
        Vec2 diff = v1 - v2;
        CHECK(sum == Vec2{4.0, 6.0});
        CHECK(diff == Vec2{2.0, 2.0});
    }

    SUBCASE("Dot product") {
        CHECK(dot(v1, v2) == doctest::Approx(11.0));
    }
}