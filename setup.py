from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Command, setup
try:
    from setuptools.command.build import build as _build
except ImportError:
    from distutils.command.build import build as _build
try:
    from setuptools.command.install import install as _install
except ImportError:
    from distutils.command.install import install as _install


ROOT = Path(__file__).resolve().parent
CMAKE_LISTS = ROOT / "CMakeLists.txt"


def default_build_dir() -> Path:
    system = os.environ.get("COMISO_BUILD_PLATFORM")
    if not system:
        if sys.platform.startswith("linux"):
            system = "linux"
        elif sys.platform == "darwin":
            system = "macos"
        elif os.name == "nt":
            system = "windows"
        else:
            system = sys.platform.replace(os.sep, "-")
    return ROOT / "build" / f"python-{system}"


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


def select_wsl_distro(explicit_name: str | None) -> str:
    if explicit_name:
        return explicit_name

    result = subprocess.run(
        ["wsl.exe", "-l", "-q"],
        check=True,
        capture_output=True,
        text=True,
    )
    distros = [line.strip().lstrip("\ufeff") for line in result.stdout.splitlines() if line.strip()]
    if not distros:
        raise RuntimeError(
            "No WSL distributions are available. Set COMISO_WSL_DISTRO or pass --wsl-distro."
        )
    return distros[0]


def windows_to_wsl_path(path: str) -> str:
    resolved = Path(path).resolve()
    drive = resolved.drive.rstrip(":").lower()
    parts = [part for part in resolved.parts[1:] if part not in ("\\", "/")]
    if not drive:
        raise RuntimeError(f"Cannot convert non-drive path to WSL format: {resolved}")
    return f"/mnt/{drive}/" + "/".join(parts)


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
            "COMISO_BUILD_DIR", str(default_build_dir())
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


class WSLCommand(Command):
    user_options = [
        ("wsl-distro=", None, "WSL distribution name"),
        ("python-executable=", None, "Python executable inside WSL"),
        ("build-dir=", None, "Out-of-source CMake build directory inside WSL"),
        ("install-prefix=", None, "CMake install prefix inside WSL"),
        ("cmake-generator=", None, "CMake generator name inside WSL"),
        ("cmake-args=", None, "Additional arguments passed to cmake configure inside WSL"),
        ("build-type=", None, "CMake build type"),
    ]

    subcommand = ""
    description = "Run a setup.py native command inside WSL"

    def initialize_options(self) -> None:
        self.wsl_distro = None
        self.python_executable = None
        self.build_dir = None
        self.install_prefix = None
        self.cmake_generator = None
        self.cmake_args = None
        self.build_type = None

    def finalize_options(self) -> None:
        self.wsl_distro = select_wsl_distro(self.wsl_distro or os.environ.get("COMISO_WSL_DISTRO"))
        self.python_executable = self.python_executable or os.environ.get(
            "COMISO_WSL_PYTHON", "python3"
        )
        self.cmake_args = self.cmake_args or os.environ.get("COMISO_CMAKE_ARGS", "")
        self.build_type = self.build_type or os.environ.get("CMAKE_BUILD_TYPE", "Release")

    def relay_args(self) -> list[str]:
        args = [self.subcommand]
        if self.build_dir:
            args.extend(["--build-dir", windows_to_wsl_path(self.build_dir)])
        if self.install_prefix:
            args.extend(["--install-prefix", windows_to_wsl_path(self.install_prefix)])
        if self.cmake_generator:
            args.extend(["--cmake-generator", self.cmake_generator])
        if self.cmake_args:
            args.extend(["--cmake-args", self.cmake_args])
        if self.build_type:
            args.extend(["--build-type", self.build_type])
        return args

    def run(self) -> None:
        if os.name != "nt":
            raise RuntimeError(f"{self.subcommand} is intended to be launched from Windows.")

        repo_root = windows_to_wsl_path(str(ROOT))
        relay = " ".join(shlex.quote(arg) for arg in self.relay_args())
        shell_command = (
            f"cd {shlex.quote(repo_root)} && "
            f"{shlex.quote(self.python_executable)} setup.py {relay}"
        )
        run_command(
            ["wsl.exe", "-d", self.wsl_distro, "sh", "-lc", shell_command],
            ROOT,
        )


class BuildNativeLinux(WSLCommand):
    description = "Run the native CoMISo Linux build inside a WSL distribution"
    subcommand = "build_native"


class InstallNativeLinux(WSLCommand):
    description = "Run the native CoMISo Linux install inside a WSL distribution"
    subcommand = "install_native"


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
        "build_native_linux": BuildNativeLinux,
        "install": InstallCommand,
        "install_native": InstallNative,
        "install_native_linux": InstallNativeLinux,
    },
)
