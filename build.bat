@echo off
echo Installing dependencies...
pip install -r requirements.txt

echo.
echo Building executable...
pyinstaller --onefile --windowed --name FluorescenceAnalysis ^
    --add-data "fluorescence_analysis.py;." ^
    app.py

echo.
echo Done! Executable is in dist\FluorescenceAnalysis.exe
pause
