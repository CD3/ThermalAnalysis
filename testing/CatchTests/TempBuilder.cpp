#include "catch.hpp"
#include "fstream"
#include <ThermalAnalysis/LinearCombinator.hpp>
#include <cmath>
#include <iostream>
#include <libField/Field.hpp>
#include <libInterpolate/Interpolate.hpp>
#include <chrono>

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
      //attempted to catch only when a value was past epsilon, abandoned, should be removed
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
      /*string interpolatorName = "GENERATED";
      string graphName = "./graphs/T1_vs_t-" + interpolatorName;
      //printing results and saving graphs
      {
      std::ofstream out("./interpResults/" + interpolatorName + ".dat");
      out << printString;
      out.flush();
      out.close();
      }

      {
      std::ofstream out(graphName);
      out << T1_vs_t;
      out.close();
      }

      {
      std::ofstream out("./graphs/Tinf_vs_t");
      out << Tinf_vs_t;
      out.close();
      }*/
     
}

TEST_CASE("Multiple Pulse Tests"){
//100 pulses, 3 seconds apart
  SECTION("Setting up profile"){
    //as in the other test case, use the same initialization in order to create a shark-fin pulse that we will use to combine 100 times
     
       
      double pulse_spacing = 3; //seconds
      double pulse_duration = 2; //seconds
      double num_pulses = 100; //pulses

      auto Tinf = [](double t) -> double {
        return sqrt(t / M_PI);
      };
      Field<double, 1> Tinf_vs_t(800);
      Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

      Field<double, 1> T1_vs_t(4000);
      T1_vs_t.setCoordinateSystem(Uniform<double>(0, pulse_spacing * num_pulses* 1.1)); //340 because 3 second seperation * 100 pulses

      //need to create orginal profile first
      LinearCombination<_1D::LinearInterpolator<double>>tempBuilder;

      tempBuilder.add(Tinf_vs_t, 0, 1);
      tempBuilder.add(Tinf_vs_t, pulse_duration, -1);


      //auto t1 = Clock::now();
      tempBuilder.build(T1_vs_t);//by ref
      //auto t2 = Clock::now();
      //std::cout << "\nTime to complete build process: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " nanoseconds\n";

      SECTION("100 pulses"){
        LinearCombination<_1D::LinearInterpolator<double>>longBuilder;
        for(int i = 0; i < num_pulses; i++){
          longBuilder.add(T1_vs_t, i * pulse_spacing, 1);
        }
        longBuilder.build(T1_vs_t);

        string interpolatorName = "100Pulses";
        string graphName = "./graphs/T1_vs_t-" + interpolatorName;
        //printing results and saving graphs
/*      {
        std::ofstream out("./interpResults/" + interpolatorName + ".dat");
        out << printString;
        out.flush();
        out.close();
        }
*/
        {
        std::ofstream out(graphName); out << T1_vs_t; out.close();
        }

        //plot the temperature profile created
        string plotString = "gnuplot -p -e \"plot '"+ graphName +"'\"";
        system((plotString).c_str());
      }
     
  }
}

TEST_CASE("Pulse Construction Tests")
{
  SECTION("Analytic Temperature Profile")
  {
    //Tinf seems to be a mapping that is useful?
    auto Tinf = [](double t) -> double {
      return sqrt(t / M_PI);
    };  // see LinearCombinations write up
    SECTION("Building finite-pulse profile")
    {
      //setup thermal profile we will build off of.
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
      LinearCombination<_1D::LinearInterpolator<double>>tempBuilder;
//    LinearCombination<_1D::CubicSplineInterpolator<double>>tempBuilder;
//    LinearCombination<_1D::MontonicInterpolator<double>>tempBuilder;

      tempBuilder.add(Tinf_vs_t, 0, 1);
      tempBuilder.add(Tinf_vs_t, 2, -1);

      tempBuilder.build(T1_vs_t);//by ref
      for (int i = 0; i < 100; i++) {
        double dt = 5. / (100 - 1);
        double t = i * dt;
        CHECK(T1_vs_t.getCoord(i) == Approx(t).epsilon(0.02));
        if (t < 2) {
          CHECK(T1_vs_t(i) == Approx(Tinf(t)).epsilon(0.02));
        } else {
          CHECK(T1_vs_t(i) == Approx(Tinf(t) - Tinf(t - 2)).epsilon(0.02));
        }
      }
     
    }
  }
}

TEST_CASE("Interpolator Comparison"){
  SECTION("Initialize Interpolators"){
    //initialize a number of LinearCombinators for their accuracy to be compared
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

    SECTION("Run Interpolation"){
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
        auto function = [Tinf](double time){return time < 2 ? Tinf(time) : (Tinf(time) - Tinf(time - 2));};
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
          double LinDiff;
          double CubicDiff;
          double MonoDiff;
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
        //skipping 1/0
        for(int i = 1; i < 100; i++){
          LinMeanErr += LinearDifferences[i];
          CubicMeanErr += CubicDifferences[i];
          MonMeanErr += MonoDifferences[i];
        }
        LinMeanErr /= 100.0;
        CubicMeanErr /= 100.0;
        MonMeanErr /= 100.0;
        string reportString = "Mean Perc. Err Linear: ";
        reportString += std::to_string(LinMeanErr) + "\n";
        reportString += "Mean Perc. Err Cubic: ";
        reportString += std::to_string(CubicMeanErr) + "\n";
        reportString += "Mean Perc. Err Monotonic: ";
        reportString += std::to_string(MonMeanErr) + "\n";
        std::cout << reportString;
        //Setting up strings to keep track of iteration number, etc
        //plot the temperature profile created
        //string plotString = "gnuplot -p -e \"plot '"+ graphName +"'\"";
        //system((plotString).c_str());
      }
    }
  }
}
