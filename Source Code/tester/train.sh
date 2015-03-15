#!/bin/sh

TRAIN_DATA_DIR=../../data/train/
TRAIN_FILE=../../data/traindata-5.txt
DETECTOR_BIN=../algo/dist/detector
OUTPUT_FOLDER=../../data/trained-data/

TRAINER_JAR=./target/tester-1.0-SNAPSHOT.jar

java -cp $TRAINER_JAR gov.nasa.asteroid.tester.AsteroidDetectorTester -folder $TRAIN_DATA_DIR"/" -train $TRAIN_FILE -exec "$DETECTOR_BIN --mode train" -output $OUTPUT_FOLDER"/"

#copy the known Neos files
cp $TRAIN_DATA_DIR/*.txt $OUTPUT_FOLDER"/"
