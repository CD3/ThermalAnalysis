#include "catch.hpp"
#include "fstream"

#include <cmath>
#include <iostream>
#include <libField/Field.hpp>
#include <libInterpolate/Interpolate.hpp>
using namespace std;

//get type for Interpolator setup
//
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
  using currentType = _1D::LinearInterpolator<double>;
  vector<double> alphas;
  vector<double> offsets;
  vector<_1D::LinearInterpolator<double>> interpolators;

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

  //simple method to encapsulate making a interpolator and setting data to a
  //field's data. Utitilizes get_data_type_t, etc
  template <class InterpType>
  auto MakeInterpolator(const Field<get_data_type_t<InterpType>,1>& F)
  {
    InterpType interp;
    interp.setData(F.size(), F.getCoordinateSystem().getAxis(0).data(), F.data());
    return interp;
  }
  //3 vectors with corresponding indexing. the profiles[i] should be multiplied
  //by a scalar multiple alphas[i], and its times offsett by offsets[i]. These
  //are member variables in order to be captured 'all in one' through capturing
  //of 'this' in the lambda expression in build

  //look into switching _1d:: (typed twice) to VVV
  //template <typename I>
  //_1D::MonotonicInterpolator<double> interp;

  
  template <typename F, typename T, typename S>
  void add(const F& field, T t, S scale)
  {
    //unimplemented
    //exchange _1D for any interpolator type
    auto interp = MakeInterpolator<_1D::LinearInterpolator<double>>(field);
    interpolators.push_back(interp);
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
    auto local_alph = (this->alphas);
    auto local_off = (this->offsets);
    auto local_interps = (this->interpolators);
    //currently replacing adding field with adding interpolator
    for(int i = 0; i < t.size(); i++){
        // i is always 0, this is just a way to unpack the list while not breaking under 
        double Temp = 0;
        for(int j = 0; j < local_interps.size(); j++){
          //offsets does not work, gnuplot starts at 0
          Temp += (local_interps[j](t[i] - offsets[j]) * local_alph[j]);
          /*string str_loc = "Currently in:\nprofile #" + std::to_string(j);
          string str_ao  = "\nalpha = " + std::to_string(local_alph[j]) + "\noffset = " + std::to_string(local_off[j]);
          string str_interp = "\nraw interp = " + std::to_string(local_interp(t[i]));
          string str_adj = "\nadj interp = " + std::to_string(local_interp(t[i]) * local_alph[j]);
          string sout = str_loc + str_ao + str_interp + str_adj + "\n\n------------------------------------\n ";
        std::cout << sout;*/
        } return Temp; }
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
