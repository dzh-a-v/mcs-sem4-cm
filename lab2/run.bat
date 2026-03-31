@echo off
if exist "build\Debug\main.exe" (
    build\Debug\main.exe
) else if exist "build\main.exe" (
    build\main.exe
) else (
    echo Program not found. Run build.bat first!
    exit /b 1
)
