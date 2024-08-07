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

  SUBCASE("Calling statistics() with 3 point throws")
  {
    sample.add(1.0);
    sample.add(2.0);
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
}