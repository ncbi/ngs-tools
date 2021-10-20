#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
NGSJAR="$HOME/ncbi-outdir/sra-tools/ngs-java/jar/ngs-java.jar"
DEPENDS="$SPARKJAR:$NGSJAR"
GOVCLASS="$HOME/ncbi-outdir/ngs-java/cls/gov/*"
NGSCLASS="$HOME/ncbi-outdir/ngs-java/cls/ngs/*"

MAINJAR="./ngs_sql.jar"

compile()
{
	#compile the java-sources into classes
	rm *.class

	echo "START COMPILING"	
	javac -Xlint:unchecked -cp $DEPENDS -d . src/*.java
	echo "END COMPILING"		
}

package()
{
	#package all compiled classes into a jar ( for execution we still need classpaths to SPARK )
	
	rm -rf ./gov ./ngs
	
	echo "START PACKING"
	
	mkdir gov ngs
	cp -r $GOVCLASS ./gov
	cp -r $NGSCLASS ./ngs
	
	#make the jar we package up contain all it's dependencies...
	jar cf $MAINJAR ./gov ./ngs *.class
	
	echo "END PACKING"
	
	rm *.class
	rm -rf ./gov ./ngs
}

compile
package
