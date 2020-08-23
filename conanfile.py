from conans import ConanFile, CMake, tools
import os, io, re

class ConanBuild(ConanFile):
    generators = "cmake"
    requires = "libField/0.7@cd3/devel", "libArrhenius/0.1@cd3/devel", "libInterpolate/2.3.1@cd3/devel", "libIntegrate/0.4.2@cd3/devel", "BoostUnitDefinitions/0.1.1@cd3/devel", "UnitConvert/0.7.1@cd3/devel", "boost/1.69.0@conan/stable", "yaml-cpp/0.6.3"

    def build_requirements(self):
      if tools.which("cmake") is None:
        self.build_requires("cmake_installer/3.16.0@conan/stable")
      else:
        output = io.StringIO()
        self.run("cmake --version",output=output)
        version = re.search( "cmake version (?P<version>\S*)", output.getvalue())
        if tools.Version( version.group("version") ) < "3.12.0":
          self.build_requires("cmake_installer/3.16.0@conan/stable")


    def build(self):
      cmake = CMake(self)
      cmake.configure()
      cmake.build()
