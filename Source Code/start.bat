@echo off
ECHO ------------------------------------------------------------------
ECHO  Automatic deployment tool - Start
ECHO.
ECHO  TCSASSEMBLER
ECHO.
ECHO ------------------------------------------------------------------

:: The application dist data path
set applicationpath=%~dp0
:: The asteroid app folder
set appfolder=asteroid
:: Full app dir
set fulldir=%applicationpath%\%appfolder%
set datadir=%applicationpath%\data
set dbdir=%datadir%\db

set warurl=%fulldir%\hunter.war
set dburl=jdbc:h2:file:%dbdir%\hunterdb;DB_CLOSE_DELAY=10

:: check if there is a built-in JRE with the app
if exist "%applicationpath%\jre\bin\java.exe" (
	set JAVA_HOME=%applicationpath%\jre
)

echo JAVA_HOME=%JAVA_HOME%
echo warurl=%warurl%
echo applicationpath=%applicationpath%
"%JAVA_HOME%\bin\java" -classpath "%fulldir%\run-jetty-1.0\lib\*" -Xms256m -Xmx2048m -XX:PermSize=128m -XX:MaxPermSize=2048m gov.nasa.asteroid.hunter.AsteroidHunterApplication "%warurl%" "%datadir%" "%fulldir%\detector\detector.exe"