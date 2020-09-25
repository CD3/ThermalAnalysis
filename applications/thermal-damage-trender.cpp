#include <ThermalAnalysis/LinearCombination.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <libField/HDF5.hpp>
#include <gputils/io.hpp>
#include <libArrhenius/Arrhenius.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

namespace po = boost::program_options;

Field<double,1> read_profile(const std::string& filename)
{
  Field<double, 1> F;
  GP2DData data;
  ReadGPASCII2DDataFile(filename, data);

  F.reset(data.x.size());
  for(int i = 0; i < F.size(); ++i)
  {
    F.getAxis(0)[i] = data.x[i];
    F(i) = data.f[i];
  }

  return F;
}

int main(int argc, const char** argv)
{
  // define command line options
  po::options_description po_opts("Options");
  po_opts.add_options()("help,h", "print help message.")(
      "verbose,v", po::value<int>()->implicit_value(0), "verbose level.")
      ("Ea", po::value<double>()->default_value(6.28e5), "Activation energy.")
      ("A",  po::value<double>()->default_value(3.1e99), "Frequency factor.")
      ;

  // now define our arguments.
  po::options_description po_args("Arguments");
  po_args.add_options()
    ("thermal-profiles"  , po::value<std::vector<std::string>>()->composing(), "Text files containing thermal profiles.")
    ;

  // combine the options and arguments into one option list.
  // this is what we will use to parse the command line, but
  // when we output a description of the options, we will just use
  // po_opts
  po::options_description all_options("Options and Arguments");
  all_options.add(po_opts).add(po_args);

  // tell boost how to translate positional options to named options
  po::positional_options_description args;
  args.add("thermal-profiles", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
      .options(all_options)
      .positional(args)
      .run(),
      vm);
  po::notify(vm);

  if (argc == 1 || vm.count("help")) {
    std::cout << R"EOL(Coming soon.
)EOL";
    return 0;
  }


  for( auto file : vm["thermal-profiles"].as<std::vector<std::string>>() )
  {
    if( !boost::filesystem::exists(file) )
    {
      std::cerr << file << " does not appear to exists. Skipping" << std::endl;
      continue;
    }

    Field<double,1> Tvst = read_profile(file);
    if(Tvst.size() < 10)
    {
      std::cerr << "Thermal profile has less than 10 data points. Cannot make any reasonable predictions about threshold trends... Skipping." << std::endl;
      continue;
    }


    auto T0 = Tvst[0];
    auto tmin = Tvst.getCoord(0);
    auto tmax = Tvst.getCoord(Tvst.size()-1);

    Tvst -= T0;

    // calculate how many 
    std::vector<double> taus;
    taus.push_back(tmax-tmin);
    while( taus[taus.size()-1] > 100e-6 )
    {
      taus.push_back(taus[taus.size()-1]/2);
    }

    double dt = 10e-6;

    libArrhenius::ThresholdCalculator< libArrhenius::ArrheniusIntegral<double> > calc(vm["A"].as<double>(),vm["Ea"].as<double>());
    for( auto tau : taus)
    {
      double tmax2 = std::min(tmin + 2*tau,tmax);
      int N = tmax2/dt;
      Field<double,1> Tvst2(N);
      Tvst2.setCoordinateSystem( Uniform(tmin,tmax2) );
      LinearCombination<_1D::MonotonicInterpolator<double>> tempBuilder;

      tempBuilder.add(Tvst, 0, 1.);
      tempBuilder.add(Tvst, tau, -1.);
      tempBuilder.build(Tvst2);

      Tvst2 += T0;

      double threshold = calc(Tvst2.size(), Tvst2.getAxis(0).data(),Tvst2.data() );
      std::cout << tau << " " << threshold << std::endl;

    }







    












  }


  return 0;
}
