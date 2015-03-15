@echo off
:: The application dist data path
set applicationpath=%~dp0
:: check if there is a built-in JRE with the app
if exist "%applicationpath%\jre\bin\java.exe" (
	set JAVA_HOME=%applicationpath%\jre
)

echo JAVA_HOME=%JAVA_HOME%

"%JAVA_HOME%\bin\java" -jar "%applicationpath%\asteroid-installer-windows.jar"