from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Command, setup
from setuptools.command.build import build as _build
from setuptools.command.install import install as _install


ROOT = Path(__file__).resolve().parent
CMAKE_LISTS = ROOT / "CMakeLists.txt"


def resolve_cmake_command() -> list[str]:
    cmake_exe = shutil.which("cmake")
    if cmake_exe:
        return [cmake_exe]

    scripts_dir = Path(sys.executable).resolve().parent
    for candidate in ("cmake.exe", "cmake"):
        cmake_path = scripts_dir / candidate
        if cmake_path.exists():
            return [str(cmake_path)]

    try:
        import cmake  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Unable to find a CMake executable. Install CMake or the Python 'cmake' package "
            "in the active environment."
        ) from exc

    cmake_bin_dir = Path(cmake.CMAKE_BIN_DIR)  # type: ignore[attr-defined]
    for candidate in ("cmake.exe", "cmake"):
        cmake_path = cmake_bin_dir / candidate
        if cmake_path.exists():
            return [str(cmake_path)]

    raise RuntimeError("The Python 'cmake' package is installed, but no CMake executable was found.")


def cmake_command(*args: str) -> list[str]:
    return [*resolve_cmake_command(), *args]


def read_version() -> str:
    content = CMAKE_LISTS.read_text(encoding="utf-8")
    match = re.search(r"project\(CoMISo VERSION ([0-9.]+)", content)
    if not match:
        raise RuntimeError("Unable to determine version from CMakeLists.txt")
    return match.group(1)


def split_cmake_args(raw: str) -> list[str]:
    return shlex.split(raw, posix=os.name != "nt") if raw else []


def run_command(args: list[str], cwd: Path) -> None:
    print(f"+ {' '.join(args)}")
    subprocess.run(args, cwd=str(cwd), check=True)


class CMakeCommand(Command):
    user_options = [
        ("build-dir=", None, "Out-of-source CMake build directory"),
        ("install-prefix=", None, "CMake install prefix"),
        ("cmake-generator=", None, "CMake generator name"),
        ("cmake-args=", None, "Additional arguments passed to cmake configure"),
        ("build-type=", None, "CMake build type"),
    ]

    def initialize_options(self) -> None:
        self.build_dir = None
        self.install_prefix = None
        self.cmake_generator = None
        self.cmake_args = None
        self.build_type = None

    def finalize_options(self) -> None:
        self.build_dir = self.build_dir or os.environ.get(
            "COMISO_BUILD_DIR", str(ROOT / "build" / "python")
        )
        self.install_prefix = self.install_prefix or os.environ.get(
            "COMISO_INSTALL_PREFIX", sys.prefix
        )
        self.cmake_generator = self.cmake_generator or os.environ.get("CMAKE_GENERATOR")
        self.cmake_args = self.cmake_args or os.environ.get("COMISO_CMAKE_ARGS", "")
        self.build_type = self.build_type or os.environ.get("CMAKE_BUILD_TYPE", "Release")

        self.build_dir = str(Path(self.build_dir).resolve())
        self.install_prefix = str(Path(self.install_prefix).resolve())

    def configure_args(self) -> list[str]:
        args = cmake_command(
            "-S",
            str(ROOT),
            "-B",
            self.build_dir,
            f"-DCMAKE_INSTALL_PREFIX={self.install_prefix}",
            f"-DCMAKE_BUILD_TYPE={self.build_type}",
            "-DCOMISO_FETCH_EIGEN=ON",
            "-DCOMISO_BUILD_EXAMPLES=OFF",
            "-DCOMISO_ENABLE_UNITTESTS=OFF",
        )
        if self.cmake_generator:
            args.extend(["-G", self.cmake_generator])
        args.extend(split_cmake_args(self.cmake_args))
        return args

    def build_args(self) -> list[str]:
        args = cmake_command("--build", self.build_dir, "--config", self.build_type)
        if os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL") is None:
            args.extend(["--parallel", str(os.cpu_count() or 1)])
        return args

    def install_args(self) -> list[str]:
        return cmake_command(
            "--install",
            self.build_dir,
            "--config",
            self.build_type,
        )


class BuildNative(CMakeCommand):
    description = "Configure and build the native CoMISo library with CMake"

    def run(self) -> None:
        run_command(self.configure_args(), ROOT)
        run_command(self.build_args(), ROOT)


class InstallNative(CMakeCommand):
    description = "Configure, build, and install the native CoMISo library with CMake"

    def run(self) -> None:
        run_command(self.configure_args(), ROOT)
        run_command(self.build_args(), ROOT)
        run_command(self.install_args(), ROOT)


class BuildCommand(_build):
    sub_commands = [("build_native", None)] + _build.sub_commands


class InstallCommand(_install):
    def run(self) -> None:
        self.run_command("install_native")
        super().run()


setup(
    name="comiso-native",
    version=read_version(),
    description="Setuptools wrapper for building and installing the CoMISo native library",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="CoMISo contributors",
    url="https://github.com/OpenFlipper-Free/CoMISo",
    license="GPL-3.0-or-later",
    packages=[],
    include_package_data=False,
    cmdclass={
        "build": BuildCommand,
        "build_native": BuildNative,
        "install": InstallCommand,
        "install_native": InstallNative,
    },
)
