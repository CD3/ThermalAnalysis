#include "catch.hpp"
#include "fstream"

#include <cmath>
#include <iostream>
#include <libField/Field.hpp>
#include <libInterpolate/Interpolate.hpp>
using namespace std;

// Writing an implementation, in plain english this implementation needs to:
// have a function add to add (function? field?) values to some container in
// mem, to be stored. There also needs a build method that will assign (not
// return) to an argument (so it must be passed by reference) a Profile
// (field?). It cannot directly add, so the interpolation library (what and
// where is it, must be \#included. Think numpy vector operations but in a real
// physical sense, keeping in mind that these array values show part of a
// larger continuous curve. What should the updated resolution of the field be?
// Be general, not just geometric time spacing may be used
class LinearCombination
{
 public:
  //3 vectors with corresponding indexing. the profiles[i] should be multiplied
  //by a scalar multiple alphas[i], and its times offsett by offsets[i]. These
  //are member variables in order to be captured 'all in one' through capturing
  //of 'this' in the lambda expression in build
  vector<Field<double, 1>> profiles;
  vector<double> alphas;
  vector<double> offsets;
  _1D::MonotonicInterpolator<double> interp;
  
  template <typename F, typename T, typename S>
  void add(const F& field, T t, S scale)
  {
    profiles.push_back(field);
    offsets.push_back(t);
    alphas.push_back(scale);
  }

  template <typename F>
  void build(F& field)
  {
    //Decription: use a lambda expression, capturing this and our argument,
    //field. In order to save some typing, initialize variables to hold the
    //various vectors, accessed through self reference. then iterate through
    //each point in t (should be 1 if 1 dim), then adjust the interpolation x
    //values by shifting them by the correct offset. set the data of the
    //interpolator object and then add it to running "Temp" returned after
    //iteration. TOMORROW: clean up, outer for loop probably unnecessary, look
    //into better way to shift values (slicing maybe?). Increase accuracy by
    //being better at interpolating
    field.set_f([this, &field](auto t) mutable{
    int tem = 0;
    auto local_prof = (this->profiles);
    auto local_alph = (this->alphas);
    auto local_off = (this->offsets);
    auto local_interp = (this->interp);
    for(int i = 0; i < t.size(); i++){
        double Temp = 0;
        for(int j = 0; j < local_prof.size(); j++){
          auto t_data = local_prof[j].getAxis(0);
          for(int k = 0; k < t_data.size(); k++){
            t_data[k] += local_off[j];
          }
          double a = local_alph[j];
//          Field<double, 1> l_p = local_prof[j];
//          local_prof[j].set_f([a, l_p](auto t) {
//            return l_p[0](t) * a;
//          });
//          local_interp.setData(t_data, local_prof[j].getData() * a);
//          Temp += local_interp(t[i]);
          local_interp.setData(t_data, local_prof[j]);
          Temp += local_interp(t[i]) * local_alph[j];
          //^^^^^^^ that is the old way, currently trying multiplying by alpha first
          /*string str_loc = "Currently in:\nprofile #" + std::to_string(j);
          string str_ao  = "\nalpha = " + std::to_string(local_alph[j]) + "\noffset = " + std::to_string(local_off[j]);
          string str_interp = "\nraw interp = " + std::to_string(local_interp(t[i]));
          string str_adj = "\nadj interp = " + std::to_string(local_interp(t[i]) * local_alph[j]);
          string sout = str_loc + str_ao + str_interp + str_adj + "\n\n------------------------------------\n ";
        std::cout << sout;*/
        }
      return Temp;
    }
    });
    return;
  }
};

TEST_CASE("TempBuilder Tests")
{
  SECTION("Analytic Temperature Profile")
  {
    //Tinf seems to be a mapping that is useful?
    auto Tinf = [](double t) -> double {
      return sqrt(t / M_PI);
    };  // see LinearCombinations write up
    SECTION("Building finite-pulse profile")
    {
      // setup thermal profile we will build off of.
      //Tinf_vs_t is a FIELD
      //Tempbuilder is a LINEARCOMBINATINO (METHOD IM DEFINING)
      Field<double, 1> Tinf_vs_t(100);
      Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      //what is the best way to include this^^ in method? not include it at
      //all, make it figure it out
      //just know it sets heat values
      Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

      // setup thermal profile that will be written to
      Field<double, 1> T1_vs_t(100);
      T1_vs_t.setCoordinateSystem(Uniform<double>(0, 5));

      // magic...
      // something like
      LinearCombination tempBuilder;

      tempBuilder.add(Tinf_vs_t, 0, 1);
      tempBuilder.add(Tinf_vs_t, 2, -1);

      tempBuilder.build(T1_vs_t);//by ref
      
      {
      std::ofstream out("T1_vs_t");
      out << T1_vs_t;
      out.close();
      }
      {
      std::ofstream out("Tinf_vs_t");
      out << Tinf_vs_t;
      out.close();
      }
           
      bool f1 = false;
      bool f2 = false;
      bool f3 = false;

      for (int i = 0; i < 100; ++i) {
        //unfinished
        if((f1 && f2) && f3){
          std::cout << "\n\nprevious failure\n\n";
        }
        double dt = 5. / (100 - 1);
        double t = i * dt;
        std::cout << std::to_string(i) + "\n";
        CHECK(T1_vs_t.getCoord(i) == Approx(t).epsilon(0.02));
        f1 = true;
        if (t < 2) {
          CHECK(T1_vs_t(i) == Approx(Tinf(t)).epsilon(0.02));
          f2 = true;
        } else {
          CHECK(T1_vs_t(i) == Approx(Tinf(t) - Tinf(t - 2)).epsilon(0.02));
          f3 = true;
        }
      f1 = false;
      f2 = false;
      f3 = false;
      }
    }
  }
}
