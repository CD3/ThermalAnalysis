#include <ThermalAnalysis/LinearCombination.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>
#define UNITCONVERT_NO_BACKWARD_COMPATIBLE_NAMESPACE
#include <UnitConvert.hpp>
#include <UnitConvert/GlobalUnitRegistry.hpp>
#include <BoostUnitDefinitions/Units.hpp>

namespace po = boost::program_options;
using namespace boost::units;

int main(int argc, const char** argv)
{

  // define command line options
  po::options_description po_opts("Options");
  po_opts.add_options()("help,h", "print help message.")(
      "verbose,v", po::value<int>()->implicit_value(0), "verbose level.");

  // now define our arguments.
  po::options_description po_args("Arguments");
  po_args.add_options()
    ("config-files"  , po::value<std::vector<std::string>>()->composing(), "Configuration files to process.")
    ;

  // combine the options and arguments into one option list.
  // this is what we will use to parse the command line, but
  // when we output a description of the options, we will just use
  // po_opts
  po::options_description all_options("Options and Arguments");
  all_options.add(po_opts).add(po_args);

  // tell boost how to translate positional options to named options
  po::positional_options_description args;
  args.add("config-files", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
      .options(all_options)
      .positional(args)
      .run(),
      vm);
  po::notify(vm);

  if (argc == 1 || vm.count("help")) {
    std::cout << R"EOL(A command line application for adding thermal profiles (Time-temperature history) to build complex profiles.

tempBuilder can be used to construct thermal profiles from complex temporal laser exposures from profiles generated from simple exposures. For
example, the thermal profile for a multiple-pulse exposure can be generated from the profile for a single-pulse exposure.
)EOL";
    return 0;
  }

  auto &ureg = UnitConvert::getGlobalUnitRegistry();

  for( auto file : vm["config-files"].as<std::vector<std::string>>() )
  {
    if( boost::filesystem::exists(file) )
    {

      YAML::Node config = YAML::LoadFile(file);


      if(!config["combinations"])
      {
        std::cout << "WARNING: No \"combinations\" section found in configuration file. Nothing will be done." << std::endl;
      }

      for( auto combination : config["combinations"] )
      {
        std::string input = combination["input"].as<std::string>();
        std::string output = combination["output"].as<std::string>();

        if(!combination["terms"])
        {
          std::cout << "WARNING: No \"terms \" section found in combination. Nothing will be done." << std::endl;
        }
        for( auto term : combination["terms"] )
        {
          // use UnitConvert library to allow users to specify quantities in any unit they want
          quantity<t::s> start_time = ureg.makeQuantity<double>( term["time"].as<std::string>() ).to<t::s>();
          auto scale = term["scale"].as<double>();


          /* read thermal profile in from input */
          /* build linear combination */
          /* write to output */
          /* vvv delete these lines vvv */
          std::cout << "start_time: " << start_time << std::endl;
          std::cout << "scale: " << scale << std::endl;
        }


      }


    }
    else{
      std::cerr << file << " does not appear to exists. Skipping" << std::endl;
    }

  }


  return 0;
}
