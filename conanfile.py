from conans import ConanFile, CMake, tools
import os

class ConanBuild(ConanFile):
    generators = "cmake_paths", "virtualenv"
    requires = "libField/0.7@cd3/devel", "libArrhenius/0.1@cd3/devel", "libInterpolate/2.3@cd3/devel", "BoostUnitDefinitions/0.1.1@cd3/devel"

    def build(self):
      cmake = CMake(self)
      cmake.configure()
      cmake.build()
