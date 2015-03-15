@echo off
ECHO ------------------------------------------------------------------
ECHO  Asteroid Data Hunter Frontend Assembly v1.0
ECHO
ECHO  TCSASSEMBLER
ECHO ------------------------------------------------------------------

:: The build dist path
set distpath=%~dp0\dist
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

ECHO  Copy start.bat to %distpath%
copy start.bat "%distpath%"\start.bat
copy start.sh "%distpath%"\start.sh

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
ECHO  start building database
ECHO ------------------------
cd hunter
call mvn sql:execute -Dh2.jdbc.url="%dburl%"

ECHO ------------------------
ECHO  DONE building database
ECHO ------------------------

:: Build war
ECHO ------------------------
ECHO  start building war
ECHO ------------------------
call mvn clean package -DskipDB=true -Dmaven.test.skip=true -Dwar.output.directory="%fulldir%"
ECHO -----------------------
ECHO  DONE building war
ECHO -----------------------

cd ..

:: Build jetty runner
ECHO ------------------------
ECHO  Start building jetty 
ECHO ------------------------
cd run-jetty
call mvn clean package -Djetty.output.directory="%fulldir%"
cd ..
ECHO ------------------------
ECHO  DONE building jetty
ECHO ------------------------

if "%1" == "installer" (
	ECHO ----------------------
	ECHO Start building Installer
	ECHO ----------------------
	
	if "%JAVA_HOME%" == "" (
		ECHO "JAVA_HOME is not set in your machine"
		exit
	)
	
	cd izpack
	call mvn clean
	mkdir target\staging
	xcopy "%JAVA_HOME%" target\staging\jre /s /e /h /i /q
	copy Unix_shortcutSpec.xml target\staging\Unix_shortcutSpec.xml
	copy Win_shortcutSpec.xml target\staging\Win_shortcutSpec.xml
	copy logo.ico target\staging\logo.ico
	copy logo.png target\staging\logo.png
	call mvn install
	cd ..
	
	if exist "installer-dist" (
		rmdir /s /q "installer-dist"
	)
	mkdir installer-dist
	copy izpack\target\izpack-maven-plugin-example-standard.jar installer-dist\asteroid-installer-windows.jar
	xcopy "%JAVA_HOME%" installer-dist\jre /s /e /h /i /q
	copy izpack\install.bat "installer-dist\Asteroid Data Hunter Installer.bat"
	
	ECHO ----------------------
	ECHO Done building Installer
	ECHO ----------------------
)

