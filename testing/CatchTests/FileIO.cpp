#include <iostream>
#include <fstream>
#include "catch.hpp"

#include <libField/Field.hpp>
#include <gputils/io.hpp>

TEST_CASE("Read field from gnuplot text file", "[field][io]")
{

    std::ofstream out("test-2d.txt");
    out << 0.1 << " " << 10 << "\n";
    out << 0.2 << " " << 20 << "\n";
    out << 0.3 << " " << 30 << "\n";
    out << 0.4 << " " << 40 << "\n";
    out.close();


    Field<double, 1> F;
    GP2DData data;
    ReadGPASCII2DDataFile("test-2d.txt", data);

    F.reset(data.x.size());
    for(int i = 0; i < F.size(); ++i)
    {
      F.getAxis(0)[i] = data.x[i];
      F(i) = data.f[i];
    }
    /* std::cout << F << "\n"; */

    CHECK(F.size() == 4);
    CHECK(F.getCoord(0) == Approx(0.1));
    CHECK(F.getCoord(1) == Approx(0.2));
    CHECK(F.getCoord(2) == Approx(0.3));
    CHECK(F.getCoord(3) == Approx(0.4));
    CHECK(F(0) == Approx(10));
    CHECK(F(1) == Approx(20));
    CHECK(F(2) == Approx(30));
    CHECK(F(3) == Approx(40));

}
