#include <ThermalAnalysis/LinearCombination.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>
#define UNITCONVERT_NO_BACKWARD_COMPATIBLE_NAMESPACE
#include <UnitConvert.hpp>
#include <UnitConvert/GlobalUnitRegistry.hpp>
#include <BoostUnitDefinitions/Units.hpp>
#include <libField/HDF5.hpp>
#include <iostream>
#include <fstream>

namespace po = boost::program_options;
using namespace boost::units;
/*
string operator>>(string &input, string name){
  string lineText;
  ifstream readFile(input);
  while(getline (readFile, lineText)){
    std::cout << lineText;
  }
  lineText << input;
  int b = 0;
  return lineText; 
}
*/

int main(int argc, const char** argv)
{
  // quick shortcut to create an input file
  /*
  Field<double, 1> delete_later(100);
  delete_later.setCoordinateSystem(Geometric<double>(0, 0.001, 1.05));
  delete_later.set_f([](auto t){return sqrt(t[0]);});
  {
    ofstream output_stream;
    output_stream.open("Tvst.txt");
    output_stream << delete_later;
    output_stream.close();
  }
  hdf5write("Tvst.h5", delete_later);
  */

  // define command line options
  po::options_description po_opts("Options");
  po_opts.add_options()("help,h", "print help message.")(
      "verbose,v", po::value<int>()->implicit_value(0), "verbose level.")(
      "write-example-config,w", po::value<std::string>(), "Write an example configuration file and exit.");

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


  if( vm.count("write-example-config") )
  {
    std::string text = R"EOL(
combinations:
  - input: 
      file_name: Tvst.h5
      format: h5
    output: 
      file_name: Tvst-1.txt
      resolution: 
        dt: 1 ms
        t_min: 0 s
        t_max: 2 s
    terms:
      - time: 0.1 s
        scale: 1
      - time: 0.2 s
        scale: 3
      - time: 0.5 s
        scale: 0.5
)EOL";

    if( vm["write-example-config"].as<std::string>() == "-" )
    {
      std::cout << text;
    }
    else{
      std::ofstream out(vm["write-example-config"].as<std::string>().c_str());
      out << text;
    }


    return 0;

  }

  auto &ureg = UnitConvert::getGlobalUnitRegistry();

  for( auto file : vm["config-files"].as<std::vector<std::string>>() )
  {

//    ReadFile >> a;

    if( boost::filesystem::exists(file) )
    {

      YAML::Node config = YAML::LoadFile(file);


      if(!config["combinations"])
      {
        std::cout << "WARNING: No \"combinations\" section found in configuration file. Nothing will be done." << std::endl;
      }

      for( auto combination : config["combinations"] )
      {
        // yaml shenanigans to obtain relevant values from config files
        auto output_field = combination["output"];
        auto input_field = combination["input"];
        auto disc_field = output_field["discretization"];
        std::string input = input_field["file_name"].as<std::string>();
        std::string input_format = input_field["format"].as<std::string>();
        std::string output = output_field["file_name"].as<std::string>();
        std::string output_format = output_field["format"].as<std::string>();
        std::string disc_type = disc_field["type"].as<std::string>();
        // discretization initializations
        //needs organization
        //as of rn all discretizers have a dt and t_min
        int N;
        Field<double, 1> T_output;
        quantity<t::s> dt = ureg.makeQuantity<double>( disc_field["dt"].as<std::string>() ).to<t::s>();
        quantity<t::s> t_min = ureg.makeQuantity<double>( disc_field["t_min"].as<std::string>() ).to<t::s>();
        quantity<t::s> t_max = ureg.makeQuantity<double>( disc_field["t_max"].as<std::string>() ).to<t::s>();
        double val_dt = dt.value();
        double val_t_min = t_min.value();
        double val_t_max = t_max.value();
        if (disc_type == "uniform"){
          N = (int)((t_max.value() - t_min.value()) / dt.value());
          T_output = Field<double, 1>(N);
          T_output.setCoordinateSystem( Uniform<double>(t_min.value(), t_max.value()) );
        } else if (disc_type == "geometric") {
          quantity<t::s> stretch = ureg.makeQuantity<double>( disc_field["stretch"].as<std::string>() ).to<t::s>();
          double s = stretch.value();
          N = (int) (log(1.0 - (( (val_t_max - val_t_min) * (1.0 - s)) / val_dt)) / log(s));
        } else if (disc_type == "geometric periodic") {
          quantity<t::s> stretch = ureg.makeQuantity<double>( disc_field["stretch"].as<std::string>() ).to<t::s>();
          quantity<t::s> period = ureg.makeQuantity<double>( disc_field["period"].as<std::string>() ).to<t::s>();
          double s = stretch.value();
          double T = period.value();
          N = (int) (((val_t_max - val_t_min) / T) * log(1.0 - (( (T - val_t_min) * (1.0 - s)) / val_dt)) / log(s));
        } else {
          std::cout << "WARNING: No recognized type found in discretization. Nothing will be done." << std::endl;
          continue;
        }
        
        LinearCombination<_1D::CubicSplineInterpolator<double>> tempBuilder;

        Field<double, 1> T_i;


        if(!combination["terms"])
        {
          std::cout << "WARNING: No \"terms \" section found in combination. Nothing will be done." << std::endl;
        }
        //need to set coordinate system before building
        
        // reading Field from file
        std::cout << "Reading from " << input;
        if("txt" == input_format ){
          //no >> operator
          //T_i << std::cin(input);
          // insert way to read from file
          std::cout << "no segfault!" << "\n";
          //T_i << std::cin(input);
          std::cout << "no segfault!" << "\n";

        } else if(input_format == std::string("h5")){
          hdf5read(input, T_i);
        } else {
          std::cout << "WARNING: No known input file format found in output. Nothing will be done." << std::endl;
          continue;
        }
        
        for( auto term : combination["terms"] )
        {
          // use UnitConvert library to allow users to specify quantities in any unit they want
          quantity<t::s> start_time = ureg.makeQuantity<double>( term["time"].as<std::string>() ).to<t::s>();
          auto scale = term["scale"].as<double>();
          tempBuilder.add(T_i, start_time.value(), scale);
        }

        tempBuilder.build(T_output);

        // Writing generated Field data to file
        if(output_format == std::string("txt")){
          ofstream output_stream;
          output_stream.open(output);
          output_stream << T_output;
          output_stream.close();
        } else if(output_format == std::string("h5")){
          hdf5write(output, T_i);
        } else {
          std::cout << "WARNING: No known output file format found in output. Nothing will be done." << std::endl;
          continue;
        }

      }


    }
    else{
      std::cerr << file << " does not appear to exists. Skipping" << std::endl;
    }

  }


  return 0;
}
