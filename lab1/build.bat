@echo off
if not exist "build" mkdir build
cd build
cmake ..
if %ERRORLEVEL% neq 0 (
    echo CMake configure failed!
    exit /b %ERRORLEVEL%
)
cmake --build .
if %ERRORLEVEL% neq 0 (
    echo Build failed!
    exit /b %ERRORLEVEL%
)
echo Build successful!
