TEST_CASE("Arrhenius Integral Test"){
  SECTION("Set Up Temperature Profiles"){
    double bodytemp = 310.0; // K
    double t_const  = 11.95800; // K
    double A    = 3.1e99;
    double E_a  = 6.28e5;// J / mol
    double R    = 8.314; // J / mol K
    libArrhenius::ArrheniusIntegral<double> Arrhenius(A, E_a);

    auto Tconst = [t_const](double t) -> double {
      return t_const;
    };  
    //setup constant temp thermal profile we will evaluate
    Field<double, 1> Tconst_vs_t(100);
    Tconst_vs_t.setCoordinateSystem(Geometric<double>(0, 0.1, 1.05));
    Tconst_vs_t.set_f([&Tconst](auto t) -> double { return Tconst(t[0]); });
    //set up arrhenius integral
    Field<double, 1> T_const_pulse_vs_t(800);
    T_const_pulse_vs_t.setCoordinateSystem(Uniform<double>(0, 5));
    LinearCombination<_1D::LinearInterpolator<double>> tempBuilder;

    tempBuilder.add(Tconst_vs_t, 0, 1);
    tempBuilder.add(Tconst_vs_t, 2, -1);
    tempBuilder.build(T_const_pulse_vs_t);

    T_const_pulse_vs_t += bodytemp;

    SECTION("Constant Temp Profile Evaulation"){

      double Omega = Arrhenius(Tconst_vs_t.getCoordinateSystem().size(0), Tconst_vs_t.getCoordinateSystem().getAxis(0).data(), Tconst_vs_t.data());
      // sub a (_a) probably means analytical but what the hell is 248.479???
      // has to be some unitless scale factor, maybe from just being off but
      // where did this come in??? was this numerical or is this a unit
      // conversion or something???
      double Omega_a = (A) * exp(-E_a / (R * (bodytemp + t_const))) * 248.479;
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
/*
    tempBuilder.add(Tinf_vs_t, 0, 1);
    tempBuilder.add(Tinf_vs_t, 2, -1);
    tempBuilder.build(T1_vs_t);
*/  

    tempBuilder.add(Tconst_vs_t, 0, 1);
    tempBuilder.add(Tconst_vs_t, 2, -1);
    tempBuilder.build(T1_vs_t);

    ofstream output;
    output.open("graphs/ConstantTest.txt");
    output << T1_vs_t;
    output.close();

    //just to set last mark
    //Damage threshold calculation should be agnostic of profile, it only happens to be passed a single pulse with the given parameters
    SECTION("Calculating Damage Threshold"){
      double U, M, L; //upper bound
      double O_U = 0;
      double O_M = 0;
      double O_L = 0;
      U = 1;
      L = 0; //arr(0) -> 0
      O_U = ScaledOmega(T1_vs_t, U, bodytemp);
      while(O_U < 1){
        U *= 2;
        O_U = ScaledOmega(T1_vs_t, U, bodytemp);
      }
      M = U;
      double dOmega_dAlpha;
      //Save something for dope #$% figurs
      string minimization_buffer;
      //redundant calculations here vvv
      while(1 != Approx(O_M).epsilon(0.01)){
        O_U = ScaledOmega(T1_vs_t, M * 1.01, bodytemp);
        O_L = ScaledOmega(T1_vs_t, M * 0.99, bodytemp);
        dOmega_dAlpha = (O_U - O_L) / (M * 0.02);
        O_M = ScaledOmega(T1_vs_t, M, bodytemp);
       M = M - (((ScaledOmega(T1_vs_t, M, bodytemp) - 1) / dOmega_dAlpha));
        CHECK(M > 0);
        std::cout << Arrhenius(T1_vs_t.getCoordinateSystem().size(0), T1_vs_t.getCoordinateSystem().getAxis(0).data(), T1_vs_t.data());
        std::cout << O_M << "\n";
      }
      std::cout << "final Omega: " << O_M << " at alpha = " << M << "\n";
      // T_o is bodytemp, M is alpha, tau is undefinied
      double tau = 2;
      double T_max = (1.0 / M) * (((E_a / R) * log(A * tau)) - bodytemp);
      std::cout << "T_max: " << std::to_string(T_max) << ", Temp rise " << std::to_string(bodytemp + t_const);
      CHECK(T_max == Approx((bodytemp + t_const)).epsilon(0.01));
    }
  }
}
