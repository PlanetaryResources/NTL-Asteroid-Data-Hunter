@echo off
ECHO ------------------------------------------------------------------
ECHO  Automatic deployment tool - Local Build
ECHO.
ECHO  TCSASSEMBLER
ECHO ------------------------------------------------------------------

:: The build dist path
set distpath=%~dp0\test-dist
:: The app folder
set appfolder=asteroid
:: Full app dir
set fulldir=%distpath%\%appfolder%


if exist "%distpath%" (
  rmdir /s /q "%distpath%"
)

mkdir "%distpath%"
mkdir "%fulldir%"
mkdir "%fulldir%\detector"
mkdir "%distpath%\data"
mkdir "%distpath%\data\hunting-data"
mkdir "%distpath%\data\db"


ECHO Copy detector to %fulldir%\detector
copy algo\dist\* "%fulldir%\detector\"

set dburl=jdbc:h2:file:%distpath%\data\db\hunter
:: Build Tester
ECHO ------------------------
ECHO  Start Building Tester
ECHO ------------------------
cd tester
call mvn clean install

ECHO ------------------------
ECHO  DONE Building Tester
ECHO ------------------------

cd ..

:: Build database
ECHO ------------------------
ECHO  Strat building database
ECHO ------------------------
cd hunter
call mvn sql:execute -Dh2.jdbc.url="%dburl%"

ECHO ------------------------
ECHO  DONE building database
ECHO ------------------------

:: Run test
ECHO ------------------------
ECHO  Start Running Test
ECHO ------------------------
call mvn clean test -DskipDB=true  -Dh2.jdbc.url="%dburl%"  -Ddetection.exe="%fulldir%\detector\detector.exe" -Ddata.base.directory="%distpath%\data\hunting-data"
ECHO -----------------------
ECHO  DONE Running Test
ECHO -----------------------

cd ..
