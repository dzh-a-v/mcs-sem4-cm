@echo off
if exist "build\Debug\interpolation.exe" (
    build\Debug\interpolation.exe
) else if exist "build\interpolation.exe" (
    build\interpolation.exe
) else (
    echo Program not found. Run build.bat first!
    exit /b 1
)
