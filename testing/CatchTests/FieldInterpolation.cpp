#include <iostream>
#include "catch.hpp"

#include <libField/Field.hpp>
#include <libInterpolate/Interpolate.hpp>

template <class Obj>
struct get_data_type {
  using type = void;
};
template <template <typename> class Obj, typename T>
struct get_data_type<Obj<T>> {
  using type = T;
};
template<typename Obj>
using get_data_type_t = typename get_data_type<Obj>::type;

template <class InterpType>
auto MakeInterpolator(const Field<get_data_type_t<InterpType>,1>& F)
{
  InterpType interp;
  interp.setData(F.size(), F.getCoordinateSystem().getAxis(0).data(), F.data());
  return interp;
}

TEST_CASE("Field Interpolation Tests")
{
  Field<double, 1> F(10);
  F.setCoordinateSystem(Uniform(0, 10));
  F.set_f([](auto x) { return 0.2 * x[0]; });

  SECTION("Manual")
  {
    _1D::LinearInterpolator<double> interp;
    interp.setData(F.size(), F.getCoordinateSystem().getAxis(0).data(),
                   F.data());

    CHECK(interp(0.0) == Approx(0.0));
    CHECK(interp(1.5) == Approx(0.3));
    CHECK(interp(9.5) == Approx(1.9));
    CHECK(interp(10) == Approx(2));
    CHECK(interp(-1) == Approx(0.0));
    CHECK(interp(11) == Approx(0.0));
  }

  SECTION("Helper Function")
  {
    auto interp = MakeInterpolator<_1D::LinearInterpolator<double>>(F);

    CHECK(interp(0.0) == Approx(0.0));
    CHECK(interp(1.5) == Approx(0.3));
    CHECK(interp(9.5) == Approx(1.9));
    CHECK(interp(10) == Approx(2));
    CHECK(interp(-1) == Approx(0.0));
    CHECK(interp(11) == Approx(0.0));
  }
}
