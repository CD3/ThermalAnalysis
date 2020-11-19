from conans import ConanFile, CMake, tools
import os, io, re

class ConanBuild(ConanFile):
    generators = "cmake"
    requires = "libField/0.8@cd3/devel", "libArrhenius/0.2.1@cd3/devel", "libInterpolate/2.3.1@cd3/devel", "libIntegrate/0.4.2@cd3/devel", "BoostUnitDefinitions/0.1.1@cd3/devel", "UnitConvert/0.7.1@cd3/devel", "boost/1.72.0", "yaml-cpp/0.6.3", "gp-utils/0.1@cd3/devel"

    default_options = {
        'boost:without_atomic': True,
        'boost:without_chrono': True,
        'boost:without_container': True,
        'boost:without_context': True,
        'boost:without_contract': True,
        'boost:without_coroutine': True,
        'boost:without_date_time': True,
        'boost:without_exception': True,
        'boost:without_fiber': True,
        'boost:without_filesystem': False,
        'boost:without_graph': True,
        'boost:without_graph_parallel': True,
        'boost:without_iostreams': True,
        'boost:without_locale': True,
        'boost:without_log': False,
        'boost:without_math': True,
        'boost:without_mpi': True,
        'boost:without_program_options': False,
        'boost:without_python': True,
        'boost:without_random': True,
        'boost:without_regex': True,
        'boost:without_serialization': True,
        'boost:without_stacktrace': True,
        'boost:without_system': False,
        'boost:without_test': True,
        'boost:without_thread': True,
        'boost:without_timer': True,
        'boost:without_type_erasure': True,
        'boost:without_wave': True,
        }
    def build(self):
      cmake = CMake(self)
      cmake.configure()
      cmake.build()
