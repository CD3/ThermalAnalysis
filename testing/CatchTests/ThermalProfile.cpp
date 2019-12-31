#include "catch.hpp"

#include <libField/Field.hpp>
#include <BoostUnitDefinitions/Units.hpp>

TEST_CASE("ThermalProfile Tests")
{
  using namespace boost;
  using namespace boost::units;
  SECTION("Configuration")
  {
    //ThermalProfile<double> prof;

    SECTION("Manual Configuration")
    {
      //prof.setOffsetTemperature( 310 * i::kelvin );
      //prof.setScaleFactor( 1.1 * i::dimless );

      Field<double,1> Tvst(3);




    }

  }
}

