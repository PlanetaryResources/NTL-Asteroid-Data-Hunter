#!/bin/bash

echo ------------------------------------------------------------------
echo  Asteroid Hunter Frontend Assembly v1.0
echo
echo  TCSASSEMBLER
echo ------------------------------------------------------------------

# The build dist path
BASEDIR=$(cd "$(dirname "$0")"; pwd)

distpath=$BASEDIR/dist
# The app folder
appfolder=asteroid
# Full app dir
fulldir=$distpath/$appfolder

rm -rf $distpath

mkdir -p "$distpath"
mkdir -p "$fulldir"
mkdir -p "$fulldir/detector"
mkdir -p "$distpath/data"
mkdir -p "$distpath/data/hunting-data"
mkdir -p "$distpath/data/db"

echo  Copy start.bat to $distpath
cp start.sh "$distpath"/start.sh
cp start.bat "$distpath"/start.bat

echo Copy detector to $fulldir/detector
if [[ "$OSTYPE" == "darwin"* ]];
then
    cp algo/dist/detector.mac "$fulldir/detector/detector"
    chmod +x "$fulldir/detector/detector"
else
    cp algo/dist/detector.ubuntu "$fulldir/detector/detector"
    chmod +x "$fulldir/detector/detector"
fi

dburl=jdbc:h2:file:$distpath/data/db/hunter
# Build Tester
echo ------------------------
echo  Start Building Tester
echo ------------------------
cd tester
 mvn clean install

echo ------------------------
echo  DONE Building Tester
echo ------------------------

cd ..

# Build database
echo ------------------------
echo  start building database
echo ------------------------
cd hunter
 mvn sql:execute -Dh2.jdbc.url="$dburl"

echo ------------------------
echo  DONE building database
echo ------------------------

# Build war
echo ------------------------
echo  start building war
echo ------------------------
 mvn clean package -DskipDB=true -Dmaven.test.skip=true -Dwar.output.directory="$fulldir"
echo -----------------------
echo  DONE building war
echo -----------------------

cd ..

# Build jetty runner
echo ------------------------
echo  Start building jetty 
echo ------------------------
cd run-jetty
 mvn clean package -Djetty.output.directory="$fulldir"
cd ..
echo ------------------------
echo  DONE building jetty
echo ------------------------

if [ "$1" == "installer" ];
then
	echo ----------------------
	echo Start building Installer
	echo ----------------------
	
	if [ "$JAVA_HOME" == "" ];
	then
		echo "JAVA_HOME is not set in your machine"
		exit
	fi
	
	cd izpack
	mvn clean
	mkdir -p target/staging
	cp -r "$JAVA_HOME" target/staging/jre 
	cp logo.* target/staging/
	cp *Spec.xml target/staging/

	mvn install
	cd ..
	
	rm -rf installer-dist
	mkdir installer-dist
	cp izpack/target/izpack-maven-plugin-example-standard.jar installer-dist/asteroid-installer-linux.jar
	cp -r "$JAVA_HOME" installer-dist/jre 
	cp izpack/install.sh "installer-dist/Asteroid_Data_Hunter_Installer.sh"
	
	echo ----------------------
	echo Done building Installer
	echo ----------------------
fi

