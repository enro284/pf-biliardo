#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "statistics.hpp"
#include "doctest.h"

// to add tests https://www.omnicalculator.com/statistics/skewness

TEST_CASE("Testing the class handling a floating point data sample")
{
  Sample sample;
  SUBCASE("Calling statistics() with no points throws")
  {
    CHECK_THROWS(sample.statistics());
  }

  SUBCASE("Calling statistics() with less than 4 point throws")
  {
    sample.add(1.0);
    CHECK_THROWS(sample.statistics());
    sample.add(2.0);
    CHECK_THROWS(sample.statistics());
    sample.add(4.0);
    CHECK_THROWS(sample.statistics());
  }

  SUBCASE("Calling statistics() with 4 points")
  {
    sample.add(1.0);
    sample.add(2.0);
    sample.add(3.0);
    sample.add(4.0);
    REQUIRE(sample.size() == 4);

    const auto result = sample.statistics();
    CHECK(result.mean == doctest::Approx(2.5));
    CHECK(result.std_dev == doctest::Approx(1.291));
    CHECK(result.skewness == doctest::Approx(0.0));
    CHECK(result.kurtosis == doctest::Approx(-1.2));
  }
  SUBCASE("Calling statistics() with 4 points")
  {
    sample.add(14.5);
    sample.add(24.76);
    sample.add(345.);
    sample.add(1.);
    sample.add(50.4);
    sample.add(67.);
    sample.add(88.);

    const auto result = sample.statistics();
    CHECK(result.mean == doctest::Approx(84.38));
    CHECK(result.std_dev == doctest::Approx(118.871));
    CHECK(result.skewness == doctest::Approx(2.2955));
    CHECK(result.kurtosis == doctest::Approx(5.5842));
  }
}