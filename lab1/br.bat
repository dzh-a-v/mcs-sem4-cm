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
cd ..
echo Build successful!
echo.
if exist "build\Debug\interpolation.exe" (
    build\Debug\interpolation.exe
) else if exist "build\interpolation.exe" (
    build\interpolation.exe
)
