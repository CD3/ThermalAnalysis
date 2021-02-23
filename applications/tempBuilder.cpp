#include <ThermalAnalysis/LinearCombination.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <yaml-cpp/yaml.h>
#define UNITCONVERT_NO_BACKWARD_COMPATIBLE_NAMESPACE
#include <UnitConvert.hpp>
#include <UnitConvert/GlobalUnitRegistry.hpp>
#include <BoostUnitDefinitions/Units.hpp>
#include <libField/HDF5.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <gputils/io.hpp>

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

struct ItemNotFound
{
  std::string item;
};

template<typename R, typename N, typename I>
auto get_(const N& node, I b, I e, boost::optional<R> def = boost::none)
{
  auto elem = *b;
  b++;
  if(node[elem])
  {
    if( b == e )
      return node[elem].template as<R>();
    
    return get_(node[elem],b,e,def);
  }

  if(def)
  {
    return def.value();
  }

  throw ItemNotFound{elem};

  return def.value();

}
template<typename R, typename N>
auto get(const N& node, std::string key, boost::optional<R> def = boost::none)
{
  std::vector<std::string> key_elements;
  boost::split(key_elements, key, boost::is_any_of("/"));
  try{
  return get_<R>(node, key_elements.begin(), key_elements.end(), def);
  }
  catch( ItemNotFound e)
  {
    std::stringstream out;
    out << "Could not find element '"<< e.item <<"' when looking for key '" << key <<"', which is required. Exiting now.";
    std::cout << out.str() << std::endl;
    exit(1);
  }
  return R();
}

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
    std::cout << po_opts << std::endl;
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
    discretization: 
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
        auto input = get<std::string>( combination, "input/file_name" );
        auto input_format = get<std::string>(combination, "input/format" );
        auto output = get<std::string>(combination, "output/file_name");
        auto output_format = get<std::string>(combination, "output/format");
        auto disc_type = get<std::string>(combination, "output/discretization/type");

        // discretization initializations
        //needs organization
        //as of rn all discretizers have a dt and t_min
        quantity<t::s> dt = ureg.makeQuantity<double>( get<std::string>(combination, "output/discretization/dt") ).to<t::s>();
        quantity<t::s> t_min = ureg.makeQuantity<double>( get<std::string>(combination, "output/discretization/t_min") ).to<t::s>();
        quantity<t::s> t_max = ureg.makeQuantity<double>( get<std::string>(combination, "output/discretization/t_max") ).to<t::s>();

        Field<double, 1> T_output;

        if (disc_type == "uniform"){
          int N = (int)((t_max.value() - t_min.value()) / dt.value());
          T_output = Field<double, 1>(N);
          T_output.setCoordinateSystem( Uniform<double>(t_min.value(), t_max.value()) );
        } else if (disc_type == "geometric") {
          quantity<t::s> stretch = ureg.makeQuantity<double>( combination["output"]["discretization"]["stretch"].as<std::string>() ).to<t::s>();
          double s = stretch.value();
          int N = (int) (log(1.0 - (( (t_max.value() - t_min.value()) * (1.0 - s)) / dt.value())) / log(s));
        } else if (disc_type == "geometric periodic") {
          quantity<t::s> stretch = ureg.makeQuantity<double>( combination["output"]["discretization"]["stretch"].as<std::string>() ).to<t::s>();
          quantity<t::s> period = ureg.makeQuantity<double>( combination["output"]["discretization"]["period"].as<std::string>() ).to<t::s>();
          double s = stretch.value();
          double T = period.value();
          int N = (int) (((t_max.value() - t_min.value()) / T) * log(1.0 - (( (T - t_min.value()) * (1.0 - s)) / dt.value())) / log(s));
        } else {
          std::cout << "WARNING: No recognized type found in discretization. Nothing will be done." << std::endl;
          continue;
        }
        
        LinearCombination<_1D::CubicSplineInterpolator<double>> tempBuilder;

        Field<double, 1> T_i;

        
        // reading Field from file
        std::cout << "Reading from " << input << std::endl;
        if("txt" == input_format ){
          GP2DData data;
          ReadGPASCII2DDataFile(input, data);

          T_i.reset(data.x.size());
          for(int i = 0; i < T_i.size(); ++i)
          {
            T_i.getAxis(0)[i] = data.x[i];
            T_i(i) = data.f[i];
          }

        } else if(input_format == std::string("h5")){
          hdf5read(input, T_i);
        } else {
          std::cout << "WARNING: No known input file format found in output. Nothing will be done." << std::endl;
          continue;
        }



        auto T0 = T_i[0];
        T_i -= T0;

        
        if(!combination["terms"])
        {
          std::cout << "WARNING: No \"terms \" section found in combination. Nothing will be done." << std::endl;
        }
        for( auto term : combination["terms"] )
        {
          // use UnitConvert library to allow users to specify quantities in any unit they want
          quantity<t::s> start_time = ureg.makeQuantity<double>( term["time"].as<std::string>() ).to<t::s>();
          auto scale = term["scale"].as<double>();
          tempBuilder.add(T_i, start_time.value(), scale);
        }

        tempBuilder.build(T_output);
        T_output += T0;

        // Writing generated Field data to file
        if(output_format == std::string("txt")){
          ofstream output_stream;
          output_stream.open(output);
          output_stream << std::setprecision(10);
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
