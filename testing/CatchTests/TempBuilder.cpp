#include <chrono>
#include <cmath>
#include "fstream"
#include <iostream>
#include "catch.hpp"
#include <libField/Field.hpp>
#include <libInterpolate/Interpolate.hpp>
#include <ThermalAnalysis/LinearCombination.hpp>
#include <libArrhenius/Integration/ArrheniusIntegral.hpp>

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

string print_results(Field<double, 1>& built_field, const std::function<double(double)>& lambda){
      string locStr = "iteration: ";
      string calcValStr = "calculated value: ";
      string actValStr = "expected value: ";
      string timeStr = "time: ";
      string printString = "";
      string divider = "";
      for(int i = 0; i < 20; i++){
        divider += "-";
      }
      divider += "\n";
      for (int i = 0; i < 100; i++){
        locStr = "iteration: ";
        timeStr = "time: ";
        calcValStr = "calculated value: ";
        actValStr = "expected value: ";
        string flagStr;
        double dt = 5. / (100 - 1);
        double t = i * dt;
        /*if(t < 2){
          flagStr = std::to_string(Tinf(t));
        } else {
          flagStr = std::to_string(Tinf(t) - Tinf(t - 2));
        }*/
        flagStr = std::to_string(lambda(t));

        locStr += std::to_string(i);
        locStr += "\n";
        timeStr += std::to_string(t);
        timeStr += "\n";
        calcValStr += std::to_string(built_field(i));
        calcValStr += "\n";
        actValStr += flagStr; 
        actValStr += "\n";
        printString += locStr + timeStr + calcValStr + actValStr + divider;
      }
      return printString;
}

// method for simplicity scaling existing constructed field
double ScaledOmega(Field<double, 1> built_field, double alpha, double T0){
    built_field *= alpha;
    built_field += T0;
    libArrhenius::ArrheniusIntegral<double> Arr(3.1e99, 6.28e5);
    return Arr(built_field.getCoordinateSystem().size(0), built_field.getCoordinateSystem().getAxis(0).data(), built_field.data());
}

TEST_CASE("Pulse Construction Tests")
{
  SECTION("Analytic Temperature Profile")
  {
    //Tinf is the basic temperature of a laser kept on "infinitely", this is a lambda function to  create that
    auto Tinf = [](double t) -> double {
      return sqrt(t / M_PI);
    };  
    SECTION("Building finite-pulse profile")
    {
      //setup thermal profile we will build off of.
      Field<double, 1> Tinf_vs_t(100);
      //Time spacing is Geometric (not linear
      Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      //Mapping 'y' values of field using lamdba Tinf
      Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

      // setup thermal profile that will be written to
      Field<double, 1> T1_vs_t(100);
      T1_vs_t.setCoordinateSystem(Uniform<double>(0, 5));

      // magic...
      // something like
      LinearCombination<_1D::LinearInterpolator<double>>tempBuilder;
//    LinearCombination<_1D::CubicSplineInterpolator<double>>tempBuilder;
//    LinearCombination<_1D::MontonicInterpolator<double>>tempBuilder;

      //add arg1 with a delayed offset of arg2 with a scalar multiplier of arg3
      tempBuilder.add(Tinf_vs_t, 0, 1);
      tempBuilder.add(Tinf_vs_t, 2, -1);

      //with all inperpolators (based on field data) on the stack, create a new profile and write to arg1
      tempBuilder.build(T1_vs_t);//by ref
      
      //perform checks of the agreement with theory of our build method
      for (int i = 0; i < 100; i++) {
        double dt = 5. / (100 - 1);
        double t = i * dt;
        //checking that the time values are correct in built profile
        CHECK(T1_vs_t.getCoord(i) == Approx(t).epsilon(0.02));
        //t = 2 corresponds to the offset where the second inf pulse is "added"
        //with alpha = -1
        if (t < 2) {
          //comparing temperature value with theory
          CHECK(T1_vs_t(i) == Approx(Tinf(t)).epsilon(0.02));
        } else {
          //comparing temperature value with theory
          CHECK(T1_vs_t(i) == Approx(Tinf(t) - Tinf(t - 2)).epsilon(0.02));
        }
      }
    }
  }
}

TEST_CASE("Multiple Pulse Tests"){
//100 pulses, 3 seconds apart
  SECTION("Build single profile"){
      //as in the other test case, use the same initialization in order to
      //create a shark-fin pulse that we will use to combine 100 times
      //Look into catch 2 for Templated test
      //You would want to make the test templated to any interpolator with required methods, and return a numerical value for error
      double pulse_spacing = 3; //seconds
      double pulse_duration = 2;//seconds
      double num_pulses = 5;  //pulses

      auto Tinf = [](double t) -> double {
        return sqrt(t / M_PI);
      };
      // set up a large amount of points for high volume testing
      Field<double, 1> Tinf_vs_t(2000);
      Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

      Field<double, 1> T1_vs_t(2000);
      T1_vs_t.setCoordinateSystem(Uniform<double>(0, pulse_spacing * num_pulses* 1.1)); //340 because 3 second seperation * 100 pulses

      LinearCombination<_1D::LinearInterpolator<double>> tempBuilder;
      tempBuilder.add(Tinf_vs_t, 0, 1);
      tempBuilder.add(Tinf_vs_t, pulse_duration, -1);
      tempBuilder.build(T1_vs_t);

      SECTION("100 pulses"){ 
        LinearCombination<_1D::LinearInterpolator<double>>longBuilder;
        for(int i = 0; i < num_pulses; i++){
          longBuilder.add(T1_vs_t, i * pulse_spacing, 1);
        }
        longBuilder.build(T1_vs_t);
        auto times = T1_vs_t.getCoordinateSystem()[0];
        for(int i = 0; i < times.size(); i++){
          double dt = 16.5 / (2000 - 1);
          double t = i * dt;
          std::cout << i << ":  " << times[i] << "\n";
          std::cout << T1_vs_t << " ->  " << times[i] << "\n";
          if(i == 0){
            continue;
          }
          std::cout << "dt:  " << times[i] - times[i - 1] << "\n";
        }
      }
  }
}

TEST_CASE("Interpolator Comparison"){
  SECTION("Initialize Interpolators"){
    //initialize a number of LinearCombinators for their accuracy to be compared
    //catch 2 template tests
    LinearCombination<_1D::LinearInterpolator<double>> Linear;
    LinearCombination<_1D::CubicSplineInterpolator<double>> Cubic;
    LinearCombination<_1D::MonotonicInterpolator<double>> Monotonic;
    //initialize pulse profiles for the various Combinators to build from (Tinf_vs_t) and write to (T1_vs_t)
    auto Tinf = [](double t) -> double {return sqrt(t / M_PI);};
    Field<double, 1> Tinf_vs_t(100);
    Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

    Field<double, 1> T1_vs_tLin(100);
    Field<double, 1> T1_vs_tCub(100);
    Field<double, 1> T1_vs_tMon(100);
    T1_vs_tLin.setCoordinateSystem(Uniform<double>(0, 5));
    T1_vs_tCub.setCoordinateSystem(Uniform<double>(0, 5));
    T1_vs_tMon.setCoordinateSystem(Uniform<double>(0, 5));

    SECTION("Check Interpolation"){
      Linear.add(Tinf_vs_t, 0, 1);
      Linear.add(Tinf_vs_t, 2, -1);
      Linear.build(T1_vs_tLin);
       
      Cubic.add(Tinf_vs_t, 0, 1);
      Cubic.add(Tinf_vs_t, 2, -1);
      Cubic.build(T1_vs_tCub);

      Monotonic.add(Tinf_vs_t, 0, 1);
      Monotonic.add(Tinf_vs_t, 2, -1);
      Monotonic.build(T1_vs_tMon);
      SECTION("Compare Results"){
        //lambda for what the correct time should be (uses ternary operator for brevity)
        auto function = [Tinf](double time){return time < 2 ? Tinf(time) : (Tinf(time) - Tinf(time - 2));};
        //references method declared at beginning
        string printLinear    = print_results(T1_vs_tLin, function);
        string printCubic     = print_results(T1_vs_tCub, function);
        string printMonotonic = print_results(T1_vs_tMon, function);
        vector<double> LinearDifferences;
        vector<double> CubicDifferences;
        vector<double> MonoDifferences;
        for (int i = 0; i < 100; i++) {
          double dt = 5. / (100 - 1);
          double t = i * dt;
          //t_value is assumed to be correct
          //difference between theory and build for each interpolator
          double LinDiff;
          double CubicDiff;
          double MonoDiff;
          //theoretical value
          double theoVal;
          if (t < 2) {
            theoVal = Tinf(t);
            LinDiff = T1_vs_tLin(i)   - theoVal;
            CubicDiff = T1_vs_tCub(i) - theoVal;
            MonoDiff = T1_vs_tMon(i)  - theoVal;
          } else {
            theoVal = (Tinf(t) - Tinf(t - 2));
            LinDiff = T1_vs_tLin(i)   - theoVal;
            CubicDiff = T1_vs_tCub(i) - theoVal; 
            MonoDiff = T1_vs_tMon(i)  - theoVal; 
          }
          //percent differences
          LinDiff /= theoVal;
          CubicDiff /= theoVal;
          MonoDiff /= theoVal;
          LinearDifferences.push_back(LinDiff);
          CubicDifferences.push_back(CubicDiff);
          MonoDifferences.push_back(MonoDiff);
        }
        double LinMeanErr;
        double CubicMeanErr;
        double MonMeanErr;
        //i = 1 to skipping 1/0 NaN that occurred in prev step, and calculating averages
        for(int i = 1; i < 100; i++){
          LinMeanErr += LinearDifferences[i];
          CubicMeanErr += CubicDifferences[i];
          MonMeanErr += MonoDifferences[i];
        }
        //99 because we skipped 0 in the previous step
        LinMeanErr /= 99.0;
        CubicMeanErr /= 99.0;
        MonMeanErr /= 99.0;
        string reportString = "Mean Perc. Err Linear: ";
        reportString += std::to_string(LinMeanErr) + "\n";
        reportString += "Mean Perc. Err Cubic: ";
        reportString += std::to_string(CubicMeanErr) + "\n";
        reportString += "Mean Perc. Err Monotonic: ";
        reportString += std::to_string(MonMeanErr) + "\n";
        std::cout << reportString;
      }
    }
  }
}

TEST_CASE("Arrhenius Integral Test"){
  SECTION("Set Up Temperature Profiles"){
    double bodytemp = 310;
    double t_const  = 11.95800;
    libArrhenius::ArrheniusIntegral<double> Arr(3.1e99, 6.28e5);
    SECTION("Constant Temp Profile Evaulation"){
      auto Tconst = [t_const](double t) -> double {
        return t_const;
      };  
      //setup constant temp thermal profile we will evaluate
      Field<double, 1> Tconst_vs_t(100);
      Tconst_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      Tconst_vs_t.set_f([&Tconst](auto t) -> double { return Tconst(t[0]); });
      //set up arrhenius integral

      Tconst_vs_t += bodytemp;
      double Omega = Arr(Tconst_vs_t.getCoordinateSystem().size(0), Tconst_vs_t.getCoordinateSystem().getAxis(0).data(), Tconst_vs_t.data());
      double Omega_a = (3.1e99) * exp(-6.28e5 / (8.314 * (bodytemp + t_const))) * 248.479;
      CHECK(Omega == Approx(Omega_a).epsilon(0.02));
    }

    // setting up single pulse profile
    auto Tinf = [](double t) -> double {
      return sqrt(t / M_PI);
    };
    Field<double, 1> Tinf_vs_t(800);
    Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });
    Field<double, 1> T1_vs_t(100);
    T1_vs_t.setCoordinateSystem(Uniform<double>(0, 5));
    LinearCombination<_1D::LinearInterpolator<double>>tempBuilder;

    tempBuilder.add(Tinf_vs_t, 0, 1);
    tempBuilder.add(Tinf_vs_t, 2, -1);
    tempBuilder.build(T1_vs_t);
    
    //Damage threshold calculation should be agnostic of profile, it only happens to be passed a single pulse with the given parameters
    SECTION("Calculating Damage Threshold"){
      double U, M, L; //upper bound
      double O_U = 0;
      double O_M = 0;
      double O_L = 0;
      U = 1;
      L = 0; //arr(0) -> 0
      O_U = ScaledOmega(T1_vs_t, U, 310);
      while(O_U < 1){
        U *= 2;
        O_U = ScaledOmega(T1_vs_t, U, 310);
      }
      M = U;
      double dOmega_dAlpha;
      while(1 != Approx(O_M).epsilon(0.01)){
        O_U = ScaledOmega(T1_vs_t, M * 1.01, 310);
        O_L = ScaledOmega(T1_vs_t, M * 0.99, 310);
        dOmega_dAlpha = (O_U - O_L) / (M * 0.02);
        O_M = ScaledOmega(T1_vs_t, M, 310);
        M = M - (((ScaledOmega(T1_vs_t, M, 310) - 1) / dOmega_dAlpha));
        CHECK(M > 0);
        std::cout << Arr(T1_vs_t.getCoordinateSystem().size(0), T1_vs_t.getCoordinateSystem().getAxis(0).data(), T1_vs_t.data());
        std::cout << O_M << "\n";
      }
      std::cout << "final Omega: " << O_M << " at alpha = " << M << "\n";
    }
  }
}


