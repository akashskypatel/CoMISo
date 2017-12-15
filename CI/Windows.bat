git submodule init
git submodule update --remote


:: set up Libraty Paths
set LIBPATH_BASE=E:/libs/%COMPILER%
set CMAKE_WINDOWS_LIBS_DIR=E:/libs


mkdir rel
cd rel

IF "%ARCHITECTURE%" == "x64" (
  set ARCH_VS= Win64
  set STRING_ARCH=64-Bit
) else (
  set ARCH_VS=
  set STRING_ARCH=32-Bit
)


IF "%QT_VERSION%" == "Qt5.3.1" (
 set QT_REV=5.3
 set QT_SUFFIX=_opengl
)

IF "%QT_VERSION%" == "Qt5.5.1" (
 set QT_REV=5.5
 set QT_SUFFIX=
)


IF "%BUILD_PLATFORM%" == "VS2013" (
    set LIBPATH=E:/libs/VS2013
    set GTESTVERSION=gtest-1.6.0
    set GENERATOR=Visual Studio 12%ARCH_VS%
    set VS_PATH="C:\Program Files (x86)\Microsoft Visual Studio 12.0\Common7\IDE\devenv.com"
    IF "%ARCHITECTURE%" == "x64" (
      set QT_BASE_CONFIG=-DQT5_INSTALL_PATH=C:\Qt\%QT_VERSION%-vs2013-%STRING_ARCH%\%QT_REV%\msvc2013_64%QT_SUFFIX% 
    )

    IF "%ARCHITECTURE%" == "x32" (
      set QT_BASE_CONFIG=-DQT5_INSTALL_PATH=C:\Qt\%QT_VERSION%-vs2013-%STRING_ARCH%\%QT_REV%\msvc2013%QT_SUFFIX%
    )

    SET BOOST_ROOT=
) 

IF "%BUILD_PLATFORM%" == "VS2015" (
    set LIBPATH=E:/libs/VS2015
    set GTESTVERSION=gtest-1.7.0
    set GENERATOR=Visual Studio 14%ARCH_VS%
    set VS_PATH="C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\IDE\devenv.com"
    set QT_BASE_CONFIG=TODO
) 

set CMAKE_CONFIGURATION=%QT_BASE_CONFIG% -DGLUT_INCLUDE_DIR="%LIBPATH%\%ARCHITECTURE%\freeglut-2.8.1\include" -DGLUT_glut_LIBRARY="%LIBPATH%\%ARCHITECTURE%\freeglut-2.8.1\lib\freeglut.lib" -DGLEW_INCLUDE_DIR="%LIBPATH%\%ARCHITECTURE%\glew-1.10.0\include" -DGLEW_LIBRARY="%LIBPATH%\%ARCHITECTURE%\glew-1.10.0\lib\glew32.lib" -DBOOST_ROOT="%LIBPATH%/%ARCHITECTURE%/boost_1_59_0" -DBOOST_LIBRARYDIR="%LIBPATH%/%ARCHITECTURE%/boost_1_59_0/lib64-msvc-12.0" -DCGAL_INCLUDE_DIR="%LIBPATH%/%ARCHITECTURE%/CGAL-4.7/include" -DCGAL_LIBRARY_DIR="%LIBPATH%/%ARCHITECTURE%/CGAL-4.7/lib" -DCGAL_BIN_DIR="%LIBPATH%/%ARCHITECTURE%/CGAL-4.7/bin" -DMUMPS_INCLUDE_DIR="%LIBPATH%/%ARCHITECTURE%/Ipopt-3.11.9/Ipopt/MSVisualStudio/v8-ifort/installed/include" -DMUMPS_LIBRARY="%LIBPATH%/%ARCHITECTURE%/Ipopt-3.11.9/Ipopt/MSVisualStudio/v8-ifort/installed/lib/CoinMumpsC.lib" 

"C:\Program Files (x86)\CMake\bin\cmake.exe" -DGTEST_PREFIX="%LIBPATH%\%ARCHITECTURE%\%GTESTVERSION%" -G "%GENERATOR%"  -DCMAKE_BUILD_TYPE=Release -DOPENFLIPPER_BUILD_UNIT_TESTS=TRUE -DCMAKE_WINDOWS_LIBS_DIR=%CMAKE_WINDOWS_LIBS_DIR% %CMAKE_CONFIGURATION% ..

IF %errorlevel% NEQ 0 exit /b %errorlevel%


%VS_PATH% /Build "Release" CoMISo.sln /Project "ALL_BUILD"

IF %errorlevel% NEQ 0 exit /b %errorlevel%

