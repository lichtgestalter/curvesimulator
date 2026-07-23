@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..\..") do set "PROJECT_ROOT=%%~fI\"

if exist "%PROJECT_ROOT%.venv\Scripts\python.exe" (
    set "PYTHON_EXE=%PROJECT_ROOT%.venv\Scripts\python.exe"
    set "USE_VENV=1"
) else (
    set "PYTHON_EXE=py -3"
    set "USE_VENV=0"
)

if defined PYTHONPATH (
    set "PYTHONPATH=%PROJECT_ROOT%src;%PYTHONPATH%"
) else (
    set "PYTHONPATH=%PROJECT_ROOT%src"
)

call :run_example "1"
if errorlevel 1 exit /b %errorlevel%

call :run_example "2"
if errorlevel 1 exit /b %errorlevel%

echo.
echo All requested examples finished successfully.
exit /b 0

:run_example
echo.
echo === Running %~1\run_example.py ===
pushd "%SCRIPT_DIR%%~1" || (
    echo Failed to change directory to %SCRIPT_DIR%%~1
    exit /b 1
)

if "%USE_VENV%"=="1" (
    "%PYTHON_EXE%" "run_example.py"
) else (
    %PYTHON_EXE% "run_example.py"
)

set "EXIT_CODE=%ERRORLEVEL%"
popd

if not "%EXIT_CODE%"=="0" (
    echo Example %~1 failed with exit code %EXIT_CODE%.
    exit /b %EXIT_CODE%
)

exit /b 0


