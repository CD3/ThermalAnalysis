#include "catch.hpp"

#include <cmath>
#include <iostream>
#include <libField/Field.hpp>

class LinearCombination
{
 public:
  template <typename F, typename T, typename S>
  void add(const F& field, T t, S scale)
  {
  }

  template <typename F>
  void build(F& field)
  {
    field.set_f([](auto t) {
      return (t[0] >= 2) ? sqrt(t[0] / M_PI) - sqrt((t[0] - 2) / M_PI)
                       : sqrt(t[0] / M_PI);
    });
  }
};

TEST_CASE("TempBuilder Tests")
{
  SECTION("Analytic Temperature Profile")
  {
    auto Tinf = [](double t) -> double {
      return sqrt(t / M_PI);
    };  // see LinearCombinations write up
    SECTION("Building finite-pulse profile")
    {
      // setup thermal profile we will build off of.
      Field<double, 1> Tinf_vs_t(100);
      Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

      // setup thermal profile that will be written to
      Field<double, 1> T1_vs_t(100);
      T1_vs_t.setCoordinateSystem(Uniform<double>(0, 5));

      // magic...
      // something like
      LinearCombination tempBuilder;

      tempBuilder.add(Tinf_vs_t, 0, 1);
      tempBuilder.add(Tinf_vs_t, 2, -1);

      tempBuilder.build(T1_vs_t);

      for (int i = 0; i < 100; ++i) {
        double dt = 5. / (100 - 1);
        double t = i * dt;
        CHECK(T1_vs_t.getCoord(i) == Approx(t));
        if (t < 2) {
          CHECK(T1_vs_t(i) == Approx(Tinf(t)));
        } else {
          CHECK(T1_vs_t(i) == Approx(Tinf(t) - Tinf(t - 2)));
        }
      }
    }
  }
}
