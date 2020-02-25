#include <vector>
#include <cmath>
#include <iostream>
#include <libField/Field.hpp>
#include <libInterpolate/Interpolate.hpp>
#include <chrono>

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

template<class currentType>
class LinearCombination
{

 public:
//using currentType = _1D::LinearInterpolator<double>;   
//using currentType = _1D::CubicSplineInterpolator<double>;
//using currentType = _1D::MonotonicInterpolator<double>; 

  vector<double> alphas;
  vector<double> offsets;
  vector<currentType*> interpolators;

// template <class currentType>
  LinearCombination(){

  }

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
  auto& MakeInterpolator(const Field<get_data_type_t<InterpType>,1>& F)
  {
    InterpType& interp;
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
    auto interp& = MakeInterpolator<currentType>(field);
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
    auto t1 = Clock::now();
    field.set_f([this, &field](auto t) mutable{
    int tem = 0;
    auto &local_alph = (this->alphas);
    auto &local_off = (this->offsets);
    auto &local_interps = (this->interpolators);
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
    auto t2 = Clock::now();
    std::cout << "\nTime to complete build process: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " nanoseconds\n";
    return;
  }

  template <typename F>
  double err(F A, F B){
    return (A - B) / B;
  }

};
