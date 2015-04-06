#!/bin/bash
echo ------------------------------------------------------------------
echo  Automatic deployment tool - Start
echo  TCSASSEMBLER
echo ------------------------------------------------------------------

# The application dist data path
BASEDIR=$(cd "$(dirname "$0")"; pwd)
applicationpath=$BASEDIR
# The asteroid app folder
appfolder=asteroid
# Full app dir
fulldir=$applicationpath/$appfolder
datadir=$applicationpath/data
dbdir=$datadir/db

warurl=$fulldir/hunter.war
dburl=jdbc:h2:file:$dbdir/hunterdb;DB_CLOSE_DELAY=10

# check if there is a built-in JRE with the app
if [ -f "$applicationpath/jre/bin/java" ] && [[ "$OSTYPE" != "darwin"* ]];
then
	export JAVA_HOME="$applicationpath/jre"
elif [[ "$OSTYPE" == "darwin"* ]];
then
    export JAVA_HOME=`/usr/libexec/java_home`
fi

echo JAVA_HOME=$JAVA_HOME
echo warurl=$warurl
echo applicationpath=$applicationpath

#stop the existing apps
INSTANCES=`ps aux | grep "java.*ASTEROID_APP"  | grep -v "grep" | awk '{print $2}'`
for INSTANCE in $INSTANCES
do
    kill -9 $INSTANCE
done

"$JAVA_HOME/bin/java" -classpath "$fulldir/run-jetty-1.0/lib/*" -DASTEROID_APP -Xms256m -Xmx2048m -XX:PermSize=128m -XX:MaxPermSize=2048m gov.nasa.asteroid.hunter.AsteroidHunterApplication "$warurl" "$datadir" "$fulldir/detector/detector"
