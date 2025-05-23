import os
from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.files import copy


class Followme(ConanFile):
    name = "followme"
    version = "1.0"
    license = ""
    settings = "os", "compiler", "build_type", "arch"

    build_requires = "cmake/[*]"

    options = {
        "cmake_build_type": ["Debug", "Release", "RelWithDebInfo", "MinSizeRel"]
    }

    default_options = {
        "cmake_build_type": "Release",
        "range-v3/*:header_only": True,
        # "fmt/*:header_only": True,
        # "spdlog/*:header_only": True,
        # "tomlplusplus/*:header_only": True,
        # "nlohmann_json/*:header_only": True,
        # "geographiclib/*:shared": True
    }

    def requirements(self):
        print("Executing requirements()")
        self.requires("range-v3/0.12.0")
        self.requires("fmt/11.2.0", override=True)
        self.requires("spdlog/1.15.1")
        self.requires("eigen/3.4.0")
        # self.requires("tomlplusplus/3.4.0")
        # self.requires("nlohmann_json/3.12.0")
        # self.requires("geographiclib/2.4")

    def valid(self):
        print("Executing valid()")
        check_min_cppstd(self, 20)

    def layout(self):
        print("Executing layout()")
        cmake_layout(self)

    def generate(self):
        print("Executing generate()")
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.generate()

        dest_bin_folder = os.path.join(self.build_folder, "bin")  # 适用于单配置生成器
        for dep in self.dependencies.values():
            for bindir in dep.cpp_info.bindirs:
                abs_bindir = os.path.join(dep.package_folder, bindir)
                print(f"绝对路径: {abs_bindir}")
                if self.settings.os == "Windows":
                    copy(self, "*.dll", src=abs_bindir, dst=dest_bin_folder, keep_path=False)
                elif self.settings.os == "Macos":
                    copy(self, "*.dylib", src=abs_bindir, dst=dest_bin_folder, keep_path=False)
                elif self.settings.os == "Linux":
                    copy(self, "*.so*", src=abs_bindir, dst=dest_bin_folder, keep_path=False)  # so* 包含版本号

    def build(self):
        print("Executing build()")
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
