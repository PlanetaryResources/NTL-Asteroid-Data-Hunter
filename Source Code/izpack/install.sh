#!/bin/bash
# The application dist data path
BASEDIR=$(cd "$(dirname "$0")"; pwd)
applicationpath=$BASEDIR

# check if there is a built-in JRE with the app
if [ -f "$applicationpath/jre/bin/java" ];
then
	export JAVA_HOME=$applicationpath/jre
fi

echo JAVA_HOME=$JAVA_HOME

"$JAVA_HOME/bin/java" -jar "$applicationpath/asteroid-installer-linux.jar"
