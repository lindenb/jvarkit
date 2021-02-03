#!/bin/bash
set -e
test -s ${HOME}/src/jsandbox/dist/velocityjson.jar

java -jar ${HOME}/src/jsandbox/dist/velocityjson.jar -Dprimitive=double array.vm > DoubleArray.java
java -jar ${HOME}/src/jsandbox/dist/velocityjson.jar -Dprimitive=float array.vm > FloatArray.java
java -jar ${HOME}/src/jsandbox/dist/velocityjson.jar -Dprimitive=byte array.vm > ByteArray.java
java -jar ${HOME}/src/jsandbox/dist/velocityjson.jar -Dprimitive=short array.vm > ShortArray.java
java -jar ${HOME}/src/jsandbox/dist/velocityjson.jar -Dprimitive=int array.vm > IntArray.java
java -jar ${HOME}/src/jsandbox/dist/velocityjson.jar -Dprimitive=long array.vm > LongArray.java
