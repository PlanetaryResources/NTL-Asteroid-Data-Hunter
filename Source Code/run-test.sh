#!/bin/sh
echo ------------------------------------------------------------------
echo  Automatic deployment tool - Local Build
echo  TCSASSEMBLER
echo ------------------------------------------------------------------

# The build dist path

BASEDIR=$(cd "$(dirname "$0")"; pwd)
distpath=$BASEDIR/dist

# The app folder
appfolder=asteroid
# Full app dir
fulldir=$distpath/$appfolder

rm -rf "$distpatch"

mkdir -p "$distpath"
mkdir -p "$fulldir"
mkdir -p "$fulldir/detector"
mkdir -p "$distpath/data"
mkdir -p "$distpath/data/hunting-data"
mkdir -p "$distpath/data/db"


echo Copy detector to $fulldir/detector
cp algo/dist/* "$fulldir/detector/"

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
echo  Strat building database
echo ------------------------
cd hunter
mvn sql:execute -Dh2.jdbc.url="$dburl"

echo ------------------------
echo  DONE building database
echo ------------------------

# Run test
echo ------------------------
echo  Start Running Test
echo ------------------------
mvn clean test -DskipDB=true  -Dh2.jdbc.url="$dburl"  -Ddetection.exe="$fulldir/detector/detector" -Ddata.base.directory="$distpath/data/hunting-data"
echo -----------------------
echo  DONE Running Test
echo -----------------------

cd ..
