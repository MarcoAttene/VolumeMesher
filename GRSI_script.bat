ECHO OFF

WHERE /Q cmake
IF %ERRORLEVEL% NEQ 0 (
	ECHO ERROR! Could not find CMake in your system path. Please verify your CMake installation.
	EXIT /B 1
)

mkdir build
cd build
cmake ..

IF %ERRORLEVEL% NEQ 0 (
	ECHO ERROR! CMake failed.
	EXIT /B 1
)

cmake --build . --config Release

IF %ERRORLEVEL% NEQ 0 (
	ECHO ERROR! Building failed.
	EXIT /B 1
)

cd Release
mesh_generator.exe -v -s -t ..\..\models\wood_fish.off D ..\..\models\hilbert.off
ren skin.off C.off
mesh_generator.exe -v -s -t ..\..\models\wood_fish.off I ..\..\models\hilbert.off
ren skin.off D.off

echo Files C.off and D.off have been created in directory build/Release
set /p DUMMY=Hit ENTER to terminate this script...
