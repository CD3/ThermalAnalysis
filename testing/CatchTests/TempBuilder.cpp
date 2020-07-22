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
// Limitations of the current tests:
//  only dealing with 1 dimensional pulses
//  physical parameters reflect retinal values
//  no way to "cut off" non ideal fields for finding damage threshold (cooling profiles/0 temp profiles)

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
    libArrhenius::ArrheniusIntegral<double> Arr(3.1e99, 6.28e5); // A and E_a should be passed in as an argument if materials change
    return Arr(built_field.getCoordinateSystem().size(0), built_field.getCoordinateSystem().getAxis(0).data(), built_field.data());
}

double damageThreshold(Field<double, 1> built_field, double T0, double precision, int rounding = 0){
    // Newton's method to minimize the Arrhenius integral, so f' is the
    // integrand and f'' is approximated through a secant line centered around
    // x
    double U, M, L; //Upper, Middle, Lower bounds
    double O_U, O_M, O_L; //"Old" UML
    double dOmega_dAlpha; // secant approximation of slope at f(M)
    double dM = 0.01; // range for secant calculation
    double goal = 1;
    switch(rounding){
      case -1:
        goal -= dM; 
      break;
      case 1:
        goal += dM; 
      break;
      default:
        //do nothing
      break;   
    }
    O_U = 0;
    O_M = 0;
    O_L = 0;
    U = 1;
    L = 0; //arr(0) -> 0
    O_U = ScaledOmega(built_field, U, T0);
    while(O_U < 1){
      U *= 2;
      O_U = ScaledOmega(built_field, U, T0);
    }
    M = U;
    while(1 != Approx(O_M).epsilon(precision)){
      O_U = ScaledOmega(built_field, M * (1. + dM), T0);
      O_L = ScaledOmega(built_field, M * (1. - dM), T0);
      dOmega_dAlpha = (O_U - O_L) / (M * (dM * 2.));
      O_M = ScaledOmega(built_field, M, T0);
      M = M - (((ScaledOmega(built_field, M, T0) - 1.) / dOmega_dAlpha));
      CHECK(M > 0);
    }
    return M;
}

template<typename T>
class GeometricPeriodImp
{
 public:
  GeometricPeriodImp(T min, T dx, T stretch, T period)
  {
    this->min     = min;
    this->dx      = dx;
    this->stretch = stretch;
    this->period = period;
    assert(stretch > 1);
    assert(period > min);
  }

  auto operator()(size_t i, size_t N) const
  {
    // x[0] = xmin
    // x[1] = xmin + dx
    // x[2] = xmin + dx + s*dx
    // x[3] = xmin + dx + s*dx + s*s*dx
    // x[i] = xmin + dx*\sigma s^(i-1)
    // which is the geometric series.
    // determine the maximum i value of the period of the function
    // which makes it periodic
    
    if(i >= N){
      i = N-1;
    }

    int i_p = static_cast<int>(log(1 - ((period - min) * (1 - stretch) / dx)) / log(stretch));
    std::cout << "maximum:" << i_p << "\n";
    std::cout << "current:" << i << "\n";
    std::cout << "modulized:" << i % i_p << "\n";
    return (period * (i / i_p)) + min + 1. * dx * (1 - pow(1. * stretch, i % i_p)) / (1 - stretch);
  }

 protected:
  T min, dx, period;
  double stretch;
};

// using a factory function here so that argument types can be deduced.
template<typename T, typename S>
GeometricPeriodImp<T> GeometricPeriod(T min, T dx, S stretch, T period)
{
  return GeometricPeriodImp<T>(min, dx, stretch, period);
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
      Field<double, 1> Tinf_vs_t(1000);
      Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
      Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

      Field<double, 1> T1_vs_t(1000);
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
          double dt = 16.5 / (1000 - 1);
          double t = i * dt;
//        std::cout << i << ":  " << times[i] << "\n";
//        std::cout << T1_vs_t << " ->  " << times[i] << "\n";
          
//        std::cout << "dt:  " << times[i] - times[i - 1] << "\n";
        }
      }
  }
}

/*TEMPLATE_TEST_CASE("Accurate interpolation", "", (_1D::LinearInterpolator<double>), (_1D::CubicSplineInterpolator<double>), (_1D::MonotonicInterpolator<double>)){
  //catch 2 template tests
  SECTION("Initialize Interpolators"){
    LinearCombination<TestType> Combinator;
    //initialize pulse profiles for the various Combinators to build from (Tinf_vs_t) and write to (T1_vs_t)
    auto Tinf = [](double t) -> double {return sqrt(t / M_PI);};
    Field<double, 1> Tinf_vs_t(100);
    Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    Tinf_vs_t.set_f([&Tinf](auto t) -> double { return Tinf(t[0]); });

    Field<double, 1> T1_vs_t(100);
    T1_vs_t.setCoordinateSystem(Uniform<double>(0, 5));

    SECTION("Check Interpolation"){
      Combinator.add(Tinf_vs_t, 0, 1);
      Combinator.add(Tinf_vs_t, 2, -1);
      Combinator.build(T1_vs_t);

      SECTION("Compare Results"){
        auto function = [Tinf](double time){return time < 2 ? Tinf(time) : (Tinf(time) - Tinf(time - 2));};
        //references method declared at beginning
        string printRes = print_results(T1_vs_t, function);
        vector<double> Differences;
        for (int i = 0; i < 100; i++) {
          double dt = 5. / (100 - 1);
          double t = i * dt;
          //t_value is assumed to be correct
          //difference between theory and build for each interpolator
          double Diff;
          //theoretical value
          double theoVal;
          if (t < 2) {
            theoVal = Tinf(t);
            Diff = T1_vs_t(i)   - theoVal;
          } else {
            theoVal = (Tinf(t) - Tinf(t - 2));
            Diff = T1_vs_t(i)   - theoVal;
          }
          CHECK(theoVal == Approx(T1_vs_t(i)).epsilon(0.02));
          //percent differences
          Diff /= theoVal;
        }
        double MeanErr;
        //i = 1 to skipping 1/0 NaN that occurred in prev step, and calculating averages
        for(int i = 1; i < 100; i++){
          MeanErr += Differences[i];
        }
        //99 because we skipped 0 in the previous step
        MeanErr /= 99.0;
        string reportString = "Mean Perc. Err Linear: ";
        reportString += std::to_string(MeanErr) + "\n";
        //std::cout << reportString;
      }
    }
  }
}
*/

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

TEST_CASE("Arrhenius"){
  SECTION("Set Up Temperature Profiles"){
    //set up physical constants
    double bodytemp = 310.0; // K
    double t_const  = 11.95800; // K
    double A    = 3.1e99;
    double E_a  = 6.28e5;// J / mol
    double R    = 8.314; // J / mol K
    double Omega; // Value for ArrheniusIntegral (Class)
    double Omega_a; //Value of ArrheniusIntegral (Analytical)
    double m;
    // enstantiate objects (Fields, ArrheniusIntegral, TempBuilder)
    Field<double, 1> Tconst_vs_t(100);
    Field<double, 1> T_const_pulse_vs_t(800);
    libArrhenius::ArrheniusIntegral<double> Arrhenius(A, E_a);
    LinearCombination<_1D::LinearInterpolator<double>> tempBuilder;
    // lambda function to return constant value
    auto Tconst = [t_const](double t) -> double {
      return t_const;
    };  
    //setup constant temp thermal profile
    Tconst_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    Tconst_vs_t.set_f([&Tconst](auto t) -> double { return Tconst(t[0]); });
    // setup field to build to
    T_const_pulse_vs_t.setCoordinateSystem(Uniform<double>(0, 5));
    //build pulse
    tempBuilder.add(Tconst_vs_t, 0, 1);
    tempBuilder.add(Tconst_vs_t, 2, -1);
    tempBuilder.build(T_const_pulse_vs_t);
    T_const_pulse_vs_t += bodytemp; //??? value is close but field should represent DeltaT???
    //Verify calculation matches analytical solution
    Omega = Arrhenius(T_const_pulse_vs_t.getCoordinateSystem().size(0), T_const_pulse_vs_t.getCoordinateSystem().getAxis(0).data(), T_const_pulse_vs_t.data());
    Omega_a = (A) * exp(-E_a / (R * (bodytemp + t_const))) * 2;
    CHECK(Omega == Approx(Omega_a).epsilon(0.02));
    T_const_pulse_vs_t -= bodytemp; //??? value is close but field should represent DeltaT???

    ofstream output;
    output.open("graphs/ConstantTest.txt");
    output << T_const_pulse_vs_t;
    output.close();

    //just to set last mark
    //Damage threshold calculation should be agnostic of profile, it only happens to be passed a single pulse with the given parameters
    
    SECTION("Calculating Damage Threshold"){
      //keeping the nearly identical SECTION instead of replacing with method
      //in case a different method is achieved, needs to be tested, more
      //verbose output is required, not because I am lazy :)
      double U, M, L; //Upper, Middle, Lower bounds
      double O_U, O_M, O_L; //"Old" UML
      double dOmega_dAlpha; // secant approximation of slope at f(M)
      double dM = 0.01; // range for secant calculation
      string minimization_buffer;
      O_U = 0;
      O_M = 0;
      O_L = 0;
      U = 1;
      L = 0; //arr(0) -> 0
      O_U = ScaledOmega(T_const_pulse_vs_t, U, bodytemp);
      while(O_U < 1){
        U *= 2;
        O_U = ScaledOmega(T_const_pulse_vs_t, U, bodytemp);
      }
      M = U;
      while(1 != Approx(O_M).epsilon(0.01)){
        O_U = ScaledOmega(T_const_pulse_vs_t, M * (1. + dM), bodytemp);
        O_L = ScaledOmega(T_const_pulse_vs_t, M * (1. - dM), bodytemp);
        dOmega_dAlpha = (O_U - O_L) / (M * (dM * 2.));
        O_M = ScaledOmega(T_const_pulse_vs_t, M, bodytemp);
        M = M - (((ScaledOmega(T_const_pulse_vs_t, M, bodytemp) - 1.) / dOmega_dAlpha));
        CHECK(M > 0);
      }
      m = M;
      std::cout << "final Omega: " << O_M << " at alpha = " << M << "\n";
    }
    double m_method = damageThreshold(T_const_pulse_vs_t, bodytemp, 0.01);

    double tau = 2;
    double T_max = (1.0 / m) * ((E_a / (R * log(A * tau))) - bodytemp);
    double log_test = log(10);
    CHECK(T_max == Approx((t_const)).epsilon(0.01));
  }
  SECTION("Multiple Pulse Arrhenius Tests"){
    double N = 64; // max number of pulses tested
    double t_o = 3; // rising edge to rising edge seperation (s)
    double tau = 2;
    double lambda = 10; // linear point density (1 / s)
    double bodytemp = 310.0; // K
    double t_const  = 11.95800; // K
    double A    = 3.1e99;
    double E_a  = 6.28e5;// J / mol
    double R    = 8.314; // J / mol K
    double m_i;
    string minimization_buffer;

    libArrhenius::ArrheniusIntegral<double> Arrhenius(A, E_a);
    LinearCombination<_1D::LinearInterpolator<double>> tempBuilder;

    Field<double, 1> Tconst_vs_t(100);
    Field<double, 1> Tinf_vs_t(100);
    Field<double, 1> T1_vs_t(200);
    Field<double, 1> T2_vs_t(200);
    Field<double, 1> N_vs_alpha(8);

    Tconst_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    Tconst_vs_t.set_f([t_const](auto t) -> double { return t_const; });
    T1_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    N_vs_alpha.setCoordinateSystem(Geometric<int>(0, 2, 2));

    tempBuilder.add(Tconst_vs_t, 0, 1);
    tempBuilder.add(Tconst_vs_t, tau, -1);
    tempBuilder.build(T1_vs_t);

    //have to start at 1 so threshold calculation hits base case
    //gonna have to define some point density and reset coordinate system each time
    //probably gonna have to create a new Field each time, unless constructor is redefinable
    
    for(int i = 2; i < N + 1; i *= 2){
      tempBuilder.clear();

      Field<double, 1> Ti_vs_t((int)(lambda * (i + 1) * t_o));
      Ti_vs_t.setCoordinateSystem(Uniform<double>(0, (i + 1) * t_o));

      for(int j = 1; j < i + 1; j++){
        tempBuilder.add(T1_vs_t, j * t_o, 1);
      }
      tempBuilder.build(Ti_vs_t);
      m_i = damageThreshold(Ti_vs_t, bodytemp, 0.01);
      
      N_vs_alpha(i) = m_i;
      minimization_buffer += std::to_string(i) + " " + std::to_string(m_i) + "\n";
      
      if(i * 2 > N){
      }
      std::cout << "\ni: " << i << "\n alpha: " << m_i << "\n";
    }
    ofstream output;
    output.open("graphs/N_vs_alpha_square.txt");
    output << minimization_buffer;
//  output << N_vs_alpha;
    output.close();

    system("echo \"PROCESS COMPLETE\"");
    
    //Create a constant field.
    //Create a single constant pulse
    //Create string and 2d vector to store pulses
    //
    //for 1 - N
    //  create field with i pulses
    //  run minimization in order to find m (alpha)
    //  add i and m to external buffer

  }
  SECTION("MPE Limit Tests"){
    std::cout << "Entered Test\n";
    double T = 0.1; //s
    double tau = 5.0e-6; // s
    int N = 1682; // max number of pulses
//  int N = 160; // max number of pulses
    double pulse_period = 5.95e-5; //s
    double PRF = 1.0 / pulse_period; // 1 / s
    double t_off = (T - N * tau) / N - 1; // falling edge to rising edge seperation (s)
    double lambda = 10; // linear point density (1 / s)
    double bodytemp = 310.0; // K
    double t_const  = 11.95800; // K
    double A    = 3.1e99;
    double E_a  = 6.28e5;// J / mol
    double R    = 8.314; // J / mol K
    double m_i;
    double alpha1, alpha2;

    Field<double, 1> Tconst_vs_t((int)((pulse_period / tau) * lambda)); // 1 / s
    Field<double, 1> Tinf_vs_t((int)(lambda * lambda)); // 1 / s^2
//  Field<double, 1> Tinf_vs_t((int)((pulse_period / tau) * lambda * lambda));
    Field<double, 1> T1_vs_t((int)((pulse_period / tau) * lambda)); // 1 / s
    Field<double, 1> T2_vs_t((int)((lambda * lambda))); // 1 / s^2
//  Field<double, 1> T2_vs_t((int)((pulse_period / tau) * lambda * lambda));
    Field<double, 1> MPE_limit_square((int)((pulse_period / tau) * lambda * N)); // 1 / s
//  Field<double, 1> MPE_limit_sharkfin((int)((pulse_period / tau) * lambda * N * 0.75)); // 1 / s
    Field<double, 1> MPE_limit_sharkfin((int)(lambda * lambda * N)); // 1 / s

    libArrhenius::ArrheniusIntegral<double> Arrhenius(A, E_a);
    LinearCombination<_1D::LinearInterpolator<double>> tempBuilder;

    Tconst_vs_t.setCoordinateSystem(Uniform<double>(0, 2 * pulse_period));
    Tconst_vs_t.set_f([t_const](auto t) -> double { return t_const; });

    Tinf_vs_t.setCoordinateSystem(Geometric<double>(0, tau * 0.5, 1.15));
//  Tinf_vs_t.setCoordinateSystem(Uniform<double>(0, N * pulse_period));
    Tinf_vs_t.set_f([t_const](auto t) -> double { return t_const * sqrt(t[0] / M_PI); });

    T1_vs_t.setCoordinateSystem(Uniform<double>(0, pulse_period));

//  T2_vs_t.setCoordinateSystem(Uniform<double>(0, N * pulse_period * 0.5));
    T2_vs_t.setCoordinateSystem(Geometric<double>(0, tau * 0.5, 1.15));

    MPE_limit_square.setCoordinateSystem(Uniform<double>(0, N * pulse_period));
//  MPE_limit_sharkfin.setCoordinateSystem(Uniform<double>(0, N * pulse_period));
    MPE_limit_sharkfin.setCoordinateSystem(GeometricPeriod<double>(0, tau * 0.5, 1.15, pulse_period));

    tempBuilder.add(Tconst_vs_t, 0, 1);
    tempBuilder.add(Tconst_vs_t, tau, -1);
    tempBuilder.build(T1_vs_t);
    tempBuilder.clear();
    std::cout << "T1_vs_t constructed";

    tempBuilder.add(Tinf_vs_t, 0, 1);
    tempBuilder.add(Tinf_vs_t, tau, -1);
    tempBuilder.build(T2_vs_t);
    tempBuilder.clear();
    std::cout << "T2_vs_t constructed";

    alpha1 = damageThreshold(T1_vs_t, bodytemp, 0.01, -1);
    alpha2 = damageThreshold(T1_vs_t, bodytemp, 0.01, -1);
    T1_vs_t *= alpha1;
    T2_vs_t *= alpha2;

    tempBuilder.add(T1_vs_t, 0, 1);
    for(int i = 1; i < N; i++){
      tempBuilder.add(i * pulse_period, 1);
    }
    tempBuilder.build(MPE_limit_square);
    tempBuilder.clear();
    std::cout << "MPE_limit_square constructed";

    tempBuilder.add(T2_vs_t, 0, 1);
    for(int i = 1; i < N; i++){
      tempBuilder.add(i * pulse_period, 1);
    }
    tempBuilder.build(MPE_limit_sharkfin);
    tempBuilder.clear();
    std::cout << "MPE_limit_sharkfin constructed";

    {
      ofstream output;
      output.open("graphs/TEST.txt");
      output << Tinf_vs_t;
      output.close();
    }
  
    {
      ofstream output;
      output.open("graphs/PotentialHazardSharkfin.txt");
      output << MPE_limit_sharkfin;
      output.close();
    }
    {
      ofstream output;
      output.open("graphs/PotentialHazardSquare.txt");
      output << MPE_limit_square;
      output.close();
    }

    system("gnuplot -p -e \'plot \"graphs/TEST.txt\" with linespoints \'");
    
    double alpha_ = damageThreshold(MPE_limit_square, bodytemp, 0.01);
    double alpha2_ = damageThreshold(MPE_limit_sharkfin, bodytemp, 0.01);
    std::cout << "\nalphas (square, shark): ( " << alpha_ << ", " << alpha2_ << ")\n";
    // double alpha_square = damageThreshold(MPE_limit_square, bodytemp, 0.01);
    //double alpha_sharkfin = damageThreshold(MPE_limit_sharkfin, bodytemp, 0.01);
  }
}

