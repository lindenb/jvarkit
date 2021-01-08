/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.FileFilter;
import java.io.PrintWriter;
import java.io.Reader;
import java.sql.Clob;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Types;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

/**

BEGIN_DOC




### Input VCF

The tool is optimized for storing very similar VCF files into an apache derby database , for example a big VCF file which would have been splitted into one VCF per transcript.



### Schema

At the time of writing this document, the current schema is:


```

CREATE TABLE ROWCONTENT(ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,MD5SUM CHAR(32) UNIQUE,CONTENT CLOB,CONTIG VARCHAR(20),FILTERED SMALLINT NOT NULL,START INT,STOP INT,ALLELE_REF VARCHAR(50));
CREATE TABLE VCF(ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,NAME VARCHAR(255));
CREATE TABLE VCFROW(ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,VCF_ID INTEGER CONSTRAINT row2vcf REFERENCES VCF,ROW_ID INTEGER CONSTRAINT row2content REFERENCES ROWCONTENT);

```


The database is created the first time the database is created. It can be a slow process.
Whole VCF lines are stored in a CBLOB.
The embedded database is local and can be removed by a simple 

```
rm -rf database.db 
```





### Inserting VCFs into the database

Inserting one VCF:


```
$ java -jar dist/vcfderby01.jar -a read -d database.db input.vcf input2.vcf.gz

#ID	NAME
1	input.vcf
2	input2.vcf.gz
```



You can insert a VCF any number of times:


```
$ java -jar dist/vcfderby01.jar -a read -d database.db input.vcf input.vcf input2.vcf.gz
#ID	NAME
3	input.vcf
4	input2.vcf.gz
```



The program also accepts concatenated VCF files:


```
$ gunzip -c input2.vcf.gz input2.vcf.gz input2.vcf.gz input2.vcf.gz | java -jar dist/vcfderby01.jar -a read -d database.db 2> /dev/null 
#ID	NAME
5	vcf1461860749175
6	vcf1461860749176
7	vcf1461860749177
8	vcf1461860749178
```





### Listing the available VCFs



```
$ java -jar dist/vcfderby01.jar -d database.db -a list
#ID	NAME	COUNT_VARIANTS
1	input.vcf	35
2	input2.vcf.gz	35
3	input.vcf	35
4	input2.vcf.gz	35
5	vcf1461860749175	35
6	vcf1461860749176	35
7	vcf1461860749177	35
8	vcf1461860749178	35

```





### Export one or more VCFs by ID

if more that one ID is given, the output is a stream of concatenated VCF.


```
$ java -jar dist/vcfderby01.jar -d database.db -a dump 5 8 3 2> /dev/null | grep "CHROM"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
```



if the output name ends with '*.zip', each VCF is saved in the zip as a new entry.


```
$ java -jar dist/vcfderby01.jar -d database.db -a dump -o out.zip 5 8 3
$ unzip -t out.zip

Archive:  out.zip
    testing: input.vcf.vcf            OK
    testing: vcf1461860749175.vcf     OK
    testing: vcf1461860749178.vcf     OK
No errors detected in compressed data of out.zip.

```




### Dumping all VCFs

In a zip:


```
$ java -jar dist/vcfderby01.jar -d database.db -a dumpall -o out.zip
$ unzip -t out.zip
Archive:  out.zip
    testing: input.vcf.vcf            OK
    testing: input2.vcf.gz.vcf        OK
    testing: 00001.ID3.vcf            OK
    testing: 00001.ID4.vcf            OK
    testing: vcf1461860749175.vcf     OK
    testing: vcf1461860749176.vcf     OK
    testing: vcf1461860749177.vcf     OK
    testing: vcf1461860749178.vcf     OK
No errors detected in compressed data of out.zip.

```



as a concatenated stream of VCFs:


```
$ java -jar dist/vcfderby01.jar -d database.db -a dumpall|\
  grep CHROM
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3

```





### Delete some VCFs by ID



```
$ java -jar dist/vcfderby01.jar -d database.db -a delete 1 8 10 
```





### Dumping all VCFs




### Accessing the database using ij

ij is an interactive SQL scripting tool that comes with Derby.
The libraries for ij is available in http://mvnrepository.com/artifact/org.apache.derby/derbytools or you can use something like:


```
sudo apt-get install derby-tools
```





```

 echo " connect 'jdbc:derby:database.db'; select count(*) from VCF; exit" |
 java -cp ./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij> ij> 1          
-----------
8          

1 row selected

```





```
$ echo " connect 'jdbc:derby:database.db'; select ID,NAME from VCF; exit" | java -cp ./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij> ij> ID         |NAME                                                                                                                            
--------------------------------------------------------------------------------------------------------------------------------------------
1          |input.vcf                                                                                                                       
2          |input2.vcf.gz                                                                                                                   
3          |input.vcf                                                                                                                       
4          |input2.vcf.gz                                                                                                                   
5          |vcf1461860749175                                                                                                                
6          |vcf1461860749176                                                                                                                
7          |vcf1461860749177                                                                                                                
8          |vcf1461860749178                                                                                                                

8 rows selected
ij>>
```



Adding a column in the VCF:


```

$ echo " connect 'jdbc:derby:database.db'; ALTER TABLE VCF ADD COLStatValue DOUBLE DEFAULT -1.0; select * from VCF; exit" | java -cp ./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij> ij> 0 rows inserted/updated/deleted
ij> ID         |NAME                                                                                                                            |MYSTATVALUE             
---------------------------------------------------------------------------------------------------------------------------------------------------------------------
1          |input.vcf                                                                                                                       |-1.0                    
2          |input2.vcf.gz                                                                                                                   |-1.0                    
3          |input.vcf                                                                                                                       |-1.0                    
4          |input2.vcf.gz                                                                                                                   |-1.0                    
5          |vcf1461860749175                                                                                                                |-1.0                    
6          |vcf1461860749176                                                                                                                |-1.0                    
7          |vcf1461860749177                                                                                                                |-1.0                    
8          |vcf1461860749178                                                                                                                |-1.0                    

8 rows selected

```



Updating the new column MYSTATVALUE:


```
echo " connect 'jdbc:derby:database.db'; UPDATE VCF SET MyStatValue =999 WHERE ID=1 OR ID=5 OR ID=8; select ID,MyStatValue from VCF; exit" | java -cp ./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij>
ij> 3 rows inserted/updated/deleted
ij> ID         |MYSTATVALUE             
------------------------------------
1          |999.0                   
2          |-1.0                    
3          |-1.0                    
4          |-1.0                    
5          |999.0                   
6          |-1.0                    
7          |-1.0                    
8          |999.0                   

8 rows selected

```



When a new column is added in a VCF, it is handled by vcderby01 and the action: 'list'



```
$ java -jar dist/vcfderby01.jar -d database.db -a list
#ID	MYSTATVALUE	NAME	COUNT_VARIANTS
1	999.0	input.vcf	35
2	-1.0	input2.vcf.gz	35
3	-1.0	input.vcf	35
4	-1.0	input2.vcf.gz	35
5	999.0	vcf1461860749175	35
6	-1.0	vcf1461860749176	35
7	-1.0	vcf1461860749177	35
8	999.0	vcf1461860749178	35
```





END_DOC
*/

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**

BEGIN_DOC


### Input VCF
The tool is optimized for storing very similar VCF files into an apache derby database , for example a big VCF file which would have been splitted int
o one VCF per transcript.


### Schema
At the time of writing this document, the current schema is:

```

CREATE TABLE ROWCONTENT(ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,MD5SUM CHAR(32) UNIQUE,CONTENT CLO
B,CONTIG VARCHAR(20),FILTERED SMALLINT NOT NULL,START INT,STOP INT,ALLELE_REF VARCHAR(50));
CREATE TABLE VCF(ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,NAME VARCHAR(255));
CREATE TABLE VCFROW(ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,VCF_ID INTEGER CONSTRAINT row2vcf REFE
RENCES VCF,ROW_ID INTEGER CONSTRAINT row2content REFERENCES ROWCONTENT);

```

The database is created the first time the database is created. It can be a slow process.
Whole VCF lines are stored in a CBLOB.
The embedded database is local and can be removed by a simple 
```
rm -rf database.db 
```



### Inserting VCFs into the database
Inserting one VCF:

```
$ java -jar dist/vcfderby01.jar -a read -d database.db input.vcf input2.vcf.gz

#ID	NAME
1	input.vcf
2	input2.vcf.gz
```


You can insert a VCF any number of times:

```
$ java -jar dist/vcfderby01.jar -a read -d database.db input.vcf input.vcf input2.vcf.gz
#ID	NAME
3	input.vcf
4	input2.vcf.gz
```


The program also accepts **concatenated** VCF files:

```
$ gunzip -c input2.vcf.gz input2.vcf.gz input2.vcf.gz input2.vcf.gz | java -jar dist/vcfderby01.jar -a read -d database.db 2> /dev/null 
#ID	NAME
5	vcf1461860749175
6	vcf1461860749176
7	vcf1461860749177
8	vcf1461860749178
```



### Listing the available VCFs

```
$ java -jar dist/vcfderby01.jar -d database.db -a list
#ID	NAME	COUNT_VARIANTS
1	input.vcf	35
2	input2.vcf.gz	35
3	input.vcf	35
4	input2.vcf.gz	35
5	vcf1461860749175	35
6	vcf1461860749176	35
7	vcf1461860749177	35
8	vcf1461860749178	35

```



### Export one or more VCFs by ID
if more that one ID is given, the output is a stream of concatenated VCF.

```
$ java -jar dist/vcfderby01.jar -d database.db -a dump 5 8 3 2> /dev/null | grep "CHROM"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
```


if the output name ends with '*.zip', each VCF is saved in the zip as a new entry.

```
$ java -jar dist/vcfderby01.jar -d database.db -a dump -o out.zip 5 8 3
$ unzip -t out.zip

Archive:  out.zip
    testing: input.vcf.vcf            OK
    testing: vcf1461860749175.vcf     OK
    testing: vcf1461860749178.vcf     OK
No errors detected in compressed data of out.zip.

```


### Dumping all VCFs
In a zip:

```
$ java -jar dist/vcfderby01.jar -d database.db -a dumpall -o out.zip
$ unzip -t out.zip
Archive:  out.zip
    testing: input.vcf.vcf            OK
    testing: input2.vcf.gz.vcf        OK
    testing: 00001.ID3.vcf            OK
    testing: 00001.ID4.vcf            OK
    testing: vcf1461860749175.vcf     OK
    testing: vcf1461860749176.vcf     OK
    testing: vcf1461860749177.vcf     OK
    testing: vcf1461860749178.vcf     OK
No errors detected in compressed data of out.zip.

```


as a concatenated stream of VCFs:

```
$ java -jar dist/vcfderby01.jar -d database.db -a dumpall|\
  grep CHROM
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1.variant	S2.variant4	S3.variant2	S4.variant3

```



### Delete some VCFs by ID

```
$ java -jar dist/vcfderby01.jar -d database.db -a delete 1 8 10 
```



### Dumping all VCFs


### Accessing the database using ij
**ij** is an interactive SQL scripting tool that comes with Derby.
The libraries for **ij** is available in [http://mvnrepository.com/artifact/org.apache.derby/derbytools](http://mvnrepository.com/artifact/org.apache.
derby/derbytools) or you can use something like:

```
sudo apt-get install derby-tools
```



```

 echo " connect 'jdbc:derby:database.db'; select count(*) from VCF; exit" |
 java -cp ./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/ap
ache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij> ij> 1          
-----------
8          

1 row selected

```



```
$ echo " connect 'jdbc:derby:database.db'; select ID,NAME from VCF; exit" | java -cp ./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.
jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derb
y.tools.ij
ij version 10.12
ij> ij> ID         |NAME                                                                                                                            
--------------------------------------------------------------------------------------------------------------------------------------------
1          |input.vcf                                                                                                                       
2          |input2.vcf.gz                                                                                                                   
3          |input.vcf                                                                                                                       
4          |input2.vcf.gz                                                                                                                   
5          |vcf1461860749175                                                                                                                
6          |vcf1461860749176                                                                                                                
7          |vcf1461860749177                                                                                                                
8          |vcf1461860749178                                                                                                                

8 rows selected
ij>>
```


Adding a column in the VCF:

```

$ echo " connect 'jdbc:derby:database.db'; ALTER TABLE VCF ADD COLStatValue DOUBLE DEFAULT -1.0; select * from VCF; exit" | java -cp ./lib/org/apache/
derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby/derbyclient/10.1
2.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij> ij> 0 rows inserted/updated/deleted
ij> ID         |NAME                                                                                                                            |MYSTA
TVALUE             
------------------------------------------------------------------------------------------------------------------------------------------------------
---------------
1          |input.vcf                                                                                                                       |-1.0     
               
2          |input2.vcf.gz                                                                                                                   |-1.0     
               
3          |input.vcf                                                                                                                       |-1.0     
               
4          |input2.vcf.gz                                                                                                                   |-1.0     
               
5          |vcf1461860749175                                                                                                                |-1.0     
               
6          |vcf1461860749176                                                                                                                |-1.0     
               
7          |vcf1461860749177                                                                                                                |-1.0     
               
8          |vcf1461860749178                                                                                                                |-1.0     
               

8 rows selected

```


Updating the new column MYSTATVALUE:

```
echo " connect 'jdbc:derby:database.db'; UPDATE VCF SET MyStatValue =999 WHERE ID=1 OR ID=5 OR ID=8; select ID,MyStatValue from VCF; exit" | java -cp 
./lib/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar:./lib/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar:./lib/org/apache/derby
/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar  org.apache.derby.tools.ij
ij version 10.12
ij>
ij> 3 rows inserted/updated/deleted
ij> ID         |MYSTATVALUE             
------------------------------------
1          |999.0                   
2          |-1.0                    
3          |-1.0                    
4          |-1.0                    
5          |999.0                   
6          |-1.0                    
7          |-1.0                    
8          |999.0                   

8 rows selected

```


When a new column is added in a VCF, it is handled by vcderby01 and the action: 'list'


```
$ java -jar dist/vcfderby01.jar -d database.db -a list
#ID	MYSTATVALUE	NAME	COUNT_VARIANTS
1	999.0	input.vcf	35
2	-1.0	input2.vcf.gz	35
3	-1.0	input.vcf	35
4	-1.0	input2.vcf.gz	35
5	999.0	vcf1461860749175	35
6	-1.0	vcf1461860749176	35
7	-1.0	vcf1461860749177	35
8	999.0	vcf1461860749178	35
```




END_DOC

 */
@Program(name="vcfderby01",
	description="Insert similar VCFs into an Apache Derby Database",
	keywords={"vcf","sql","derby","burden"})
public class VcfDerby01
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfDerby01.class).make();
	private StringToMd5 toMd5 = new StringToMd5();
	
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	
	
	@Parameter(names={"-d","--derby"},description="REQUIRED. path to Derby database storage directory.")
	private String derbyFilePath = "";
	
	@Parameter(names={"-a","--action"},description="REQUIRED. action to perform. 'read': read a zip or a concatenated stream of vcf files and insert it into a derby database. 'list': list the available vcf. 'dump' dump one or more VCF. 'dumpall' dump all VCFs. 'dumpuniq' dum all as a one and only uniq vcf. 'delete' : delete one or more VCF by ID")
	private String actionStr = "";
	
	@Parameter(names={"-t","--title"},description="Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the name of the inserted VCF file")
	private String titleHeaderStr = "";

	
	private static int MAX_REF_BASE_LENGTH=50;
	private long ID_GENERATOR = System.currentTimeMillis();
	private Connection conn=null;
	private static final String VCF_HEADER_FILE_ID="##VcfDerby01VcfId=";
	private static final String VCF_HEADER_FILE_NAME="##VcfDerby01VcfName=";
	public VcfDerby01()
		{
		}
	 
	
	private File getDerbyDirectory() {
		return new File(this.derbyFilePath);
	}
	
	
	private void openDerby() {
		try {
			boolean create;
			final Properties props = new Properties();
			final File derbyDir = getDerbyDirectory();
			LOG.info("open derby :" + getDerbyDirectory());
			if(derbyDir.exists()) {
				if(!derbyDir.isDirectory()) {
					throw new RuntimeIOException("derby database is not a directory : "+derbyDir);
					}
				
				if( derbyDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						if(pathname.isFile()){
							if( pathname.getName().equals("service.properties")) return true;
							if( pathname.getName().equals("README_DO_NOT_TOUCH_FILES.txt")) return true;
							}
						if(pathname.isDirectory()){
							if( pathname.getName().startsWith("log")) return true;
						}
						return false;
					}
				}).length==0) {
					throw new RuntimeIOException("derby database exist but doesn't look like a derby directory : "+derbyDir);
				}
				
				create=false;
				}
			else
				{
				create=true;
				}
			props.setProperty("create", String.valueOf(create));
			this.conn = DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
			
			if(create) {
				final String tableId = "ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY";
				final Statement stmt= this.conn.createStatement();
				final String sqls[]={
						"CREATE TABLE ROWCONTENT("+tableId+",MD5SUM CHAR(32) UNIQUE,CONTENT CLOB,CONTIG VARCHAR(20),FILTERED SMALLINT NOT NULL,START INT,STOP INT,ALLELE_REF VARCHAR("+MAX_REF_BASE_LENGTH+"))",
						"CREATE TABLE VCF("+tableId+",NAME VARCHAR(255))",
						"CREATE TABLE VCFROW("+tableId+",VCF_ID INTEGER CONSTRAINT row2vcf REFERENCES VCF,ROW_ID INTEGER CONSTRAINT row2content REFERENCES ROWCONTENT)"
						};
				for(final String sql:sqls) {
					LOG.warn(sql);
					stmt.execute(sql);
				}
				stmt.close();
			}
			this.conn.setAutoCommit(true);
		} catch (Exception e) {
			CloserUtil.close(this.conn);
			this.conn = null;
			throw new RuntimeException(e);
		}
	}
	
	private static long getLastGeneratedId(final PreparedStatement pstmt) throws SQLException {
			ResultSet keys = null;
			long id = -1L;
			try {
			keys = pstmt.getGeneratedKeys();
		    while (keys.next()) {
		       if(id!=-1L) throw new IllegalStateException();
		       id= keys.getLong(1);
		       }	 
			if(id==-1L) throw new IllegalStateException("No SQL ID was found");
			return id;
			}
			finally {
				CloserUtil.close(keys);
			}
		}
	
	private void closeDerby() {
		CloserUtil.close(this.conn);
		this.conn = null;
		try {
			final Properties props = new Properties();
			props.setProperty("shutdown", "true");
			final File derbyDir = getDerbyDirectory();
			DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
		} catch (Exception e) {
			
			}
		}
	
	private void compress() {
		PreparedStatement pstmt =null;
		
		for(final String sqlStmt :new String[]{
				"call SYSCS_UTIL.SYSCS_COMPRESS_TABLE(?,?, 1)",	
				"call SYSCS_UTIL.SYSCS_INPLACE_COMPRESS_TABLE(?,?, 1, 1, 1)"
		}) {
		try {
		/* compress data */
		pstmt = this.conn.prepareStatement(sqlStmt);
		pstmt.setString(1, "APP");
		for(final String table: new String[]{"ROWCONTENT","VCF","VCFROW"}) {
			LOG.info("compressing "+table+" with "+sqlStmt);
			pstmt.setString(2,table);
			pstmt.execute();
		}
		pstmt.close();
		} catch(Exception err) {
			LOG.warn("Cannot compress" + err.getMessage());
		}
		finally 
		{
			CloserUtil.close(pstmt);
		}
		}

		
		
		
	}
	
	private int doCommandDumpAll(final List<String> args){
		Statement stmt = null;
		ResultSet row = null;
		final Set<Long> vcfids = new TreeSet<>();
		if(!args.isEmpty()) {
			LOG.error("Too many arguments");
			return -1;
		}
		try {
			stmt = this.conn.createStatement();
			row = stmt.executeQuery("SELECT ID from VCF");
			while(row.next()) {
				vcfids.add(row.getLong(1));
			}
			row.close();row=null;
			stmt.close();stmt=null;
			return dump(vcfids);
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(row);
			CloserUtil.close(stmt);
		}
	}
	
	private int doCommandDumpUniq(final List<String> args){
		if(!args.isEmpty()) {
			LOG.error("Too many arguments");
			return -1;
		}
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		try {
			boolean chrom_line_seen=false;
			for(int side=0;side<2;++side)
				{
				final String sql=(side==0?
						"SELECT ROWCONTENT.CONTENT FROM ROWCONTENT WHERE ROWCONTENT.CONTIG IS NULL ORDER BY ROWCONTENT.ID " :
						"SELECT ROWCONTENT.CONTENT FROM ROWCONTENT WHERE ROWCONTENT.CONTIG IS NOT NULL ORDER BY ROWCONTENT.CONTIG,ROWCONTENT.START,ROWCONTENT.ALLELE_REF "
						);
				LOG.info(sql);
				pstmt2 = this.conn.prepareStatement(sql);
				pwOut = openFileOrStdoutAsPrintWriter(this.outputFile);
				row =  pstmt2.executeQuery();
				while(row.next()) {
					final Clob clob = row.getClob(1);
					final Reader r= clob.getCharacterStream();
					if(side==0)
						{
						final String s = IOUtils.copyToString(r);
						if(s.startsWith("#CHROM")) {
							if(chrom_line_seen) {
								r.close();
								continue;
							}
							chrom_line_seen=true;
						} else
							{
							if(chrom_line_seen) {
								r.close();
								continue;
								}
							}
						pwOut.print(s);
						}
					else
						{
						IOUtils.copyTo(r,pwOut);
						}	
					
					
					r.close();
					pwOut.println();
					}
				row.close();
					
				pstmt2.close();pstmt2=null;
				
				pwOut.flush();
				}
			pwOut.close();pwOut=null;
			
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt2);
		}
	}

	@SuppressWarnings("unused")
	private int ___dumpAll(List<String> args){
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		int num_vcf_exported=0;
		try {
			pstmt2 = this.conn.prepareStatement(
					"SELECT ROWCONTENT.CONTENT,VCF.ID,VCF.NAME FROM VCF,VCFROW,ROWCONTENT WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND ORDER BY VCF.ID,VCFROW.ID ");
			pwOut = openFileOrStdoutAsPrintWriter(this.outputFile);
			row =  pstmt2.executeQuery();
			final String CHROM_prefix="#CHROM\t";
			final StringBuilder chrom_header_line = new StringBuilder(CHROM_prefix.length());
			while(row.next()) {
				final Clob clob = row.getClob(1);
				final Reader r= clob.getCharacterStream();
				/* read the first bytes to check if it's the #CHROM line
				 * if true, add a VCF header line with VCF ID and NAME
				 *  */
				chrom_header_line.setLength(0);
				int c;
				while(chrom_header_line.length()< CHROM_prefix.length() &&
					((c=r.read())!=1) )
					{
					chrom_header_line.append((char)c);
					}
				final String amorce = chrom_header_line.toString();
				if(amorce.equals(CHROM_prefix))
					{
					num_vcf_exported++;
					final long vcf_id = row.getLong(2);
					final String vcfName= row.getString(3);
					pwOut.println(VCF_HEADER_FILE_ID+vcf_id);
					pwOut.println(VCF_HEADER_FILE_NAME+vcfName);
					}
				pwOut.print(amorce);
				IOUtils.copyTo(r,pwOut);
				r.close();
				pwOut.println();
			}
			row.close();
				
			pstmt2.close();
			
			pwOut.flush();
			pwOut.close();
			
			if(num_vcf_exported==0 ) {
				LOG.warn("NO VCF WAS EXPORTED");
				}
			LOG.info("count(VCF) exported "+num_vcf_exported);
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt2);
		}
	}
	

	private int dump(final Set<Long> vcfIds){
		final double timeStart = System.currentTimeMillis();

		PreparedStatement pstmt = null;
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		int num_vcf_exported=0;
		try {
			pstmt = this.conn.prepareStatement("SELECT NAME from VCF where ID=?");
			pstmt2 = this.conn.prepareStatement("SELECT ROWCONTENT.CONTENT FROM VCF,VCFROW,ROWCONTENT WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND VCF.ID=? ORDER BY VCFROW.ID ");
			
			
			pwOut = openFileOrStdoutAsPrintWriter(this.outputFile);
			
			
			for(final long vcf_id : vcfIds)
				{
				final double remain = ((((System.currentTimeMillis()-timeStart)/(1+num_vcf_exported)))*(vcfIds.size()-(1+num_vcf_exported)))/1000.0;
				
				LOG.info("Getting VCF "+vcf_id+" "+(num_vcf_exported)+"/"+vcfIds.size() +" . Remains "+(long)remain+" seconds.");
				
				pstmt.setLong(1, vcf_id);
				String vcfName=null;
				row = pstmt.executeQuery();
				while(row.next()) {
					vcfName = row.getString(1);
				}
				row.close();
				if(vcfName==null) 
					{
					LOG.error("Cannot find VCF ID "+vcf_id);
					return -1;
					}
				LOG.info("dumping "+vcfName+" ID:"+vcf_id);
				pstmt2.setLong(1, vcf_id);
				row =  pstmt2.executeQuery();
				
				final String CHROM_prefix="#CHROM\t";
				while(row.next()) {
					final Clob clob = row.getClob(1);
					final Reader r= clob.getCharacterStream();
					/* read the first bytes to check if it's the #CHROM line
					 * if true, add a VCF header line with VCF ID and NAME
					 *  */
					final StringBuilder chrom_header_line = new StringBuilder(CHROM_prefix.length());
					int c;
					while(chrom_header_line.length()< CHROM_prefix.length() &&
						((c=r.read())!=1) )
						{
						chrom_header_line.append((char)c);
						}
					final String amorce = chrom_header_line.toString();
					if(amorce.equals(CHROM_prefix))
						{
						pwOut.println(VCF_HEADER_FILE_ID+vcf_id);
						pwOut.println(VCF_HEADER_FILE_NAME+vcfName);
						}
					pwOut.print(amorce);
					IOUtils.copyTo(r,pwOut);
					r.close();
					pwOut.println();
				}
				
				
				row.close();
				num_vcf_exported++;
				
				}
				
			pstmt.close();
			pstmt2.close();
			
			pwOut.flush();
			pwOut.close();
			
			if(num_vcf_exported==0 ) {
				LOG.warn("NO VCF WAS EXPORTED");
			}
			LOG.info("count(VCF) exported "+num_vcf_exported);
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			CloserUtil.close(pstmt2);
		}
	}

	
	
	private int doCommandDump(List<String> args){
		final Pattern comma = Pattern.compile("[,]");
		final Set<Long> vcfIds = new TreeSet<>();
		try {
			for(final String token: args) {
				for(final String idStr: comma.split(token))
					{
					if(idStr.trim().isEmpty()) continue;
					LOG.info("Getting VCF "+idStr);
					long vcf_id = -1L;
					try {
						vcf_id = Long.parseLong(idStr);
					} catch(final NumberFormatException err) {
						LOG.error("Bad VCF ID :"+idStr);
						return -1;
					}
					vcfIds.add(vcf_id);
					}
				}
			return dump(vcfIds);
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
		}
	}

	
	private int doReadConcatenatedVcf(List<String> args){
		int number_of_ref_allele_truncated=0;
		PreparedStatement pstmt = null;
		PreparedStatement pstmt2 = null;
		PreparedStatement pstmt3 = null;
		ResultSet row = null;
		PrintWriter pw = null;
		args = new ArrayList<>(IOUtils.unrollFiles(args));
		LOG.info(args.toString());
		LineIterator lineIter=null;
		final String titleHeaderTag = (
				this.titleHeaderStr==null || this.titleHeaderStr.trim().isEmpty()?
				null:
				"##"+titleHeaderStr+"="
				);
		try {
			int fileidx=0;
			
			pw = openFileOrStdoutAsPrintWriter(this.outputFile);
			pw.println("#ID\tNAME");

			do
			{
				if(fileidx==0 && args.isEmpty()) {
					lineIter = IOUtils.openStreamForLineIterator(stdin());
				} else
				{
					lineIter = IOUtils.openURIForLineIterator(args.get(fileidx));
				}
				int num_vcf_in_this_stream = 0;
				while(lineIter.hasNext()) {
					String filename= "vcf"+(++ID_GENERATOR);
					if(num_vcf_in_this_stream==0 && !args.isEmpty()) {
						filename = args.get(fileidx);
					}
					
					final List<String> headerLines = new ArrayList<>();
					while(lineIter.hasNext() && lineIter.peek().startsWith("#")) {
						final String h= lineIter.next();
						if( h.startsWith(VCF_HEADER_FILE_ID) ||h.startsWith(VCF_HEADER_FILE_NAME)) {
							LOG.info("Ignoring line "+h);
							continue;
						}
						/* find filename in vcf header */
						if( titleHeaderTag!=null &&
							h.startsWith(titleHeaderTag) &&
							h.trim().length()>titleHeaderTag.length()) {
							filename = h.substring(titleHeaderTag.length()).trim();
						}
						
						headerLines.add(h);
					}
					final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(headerLines);
					
					pstmt = this.conn.prepareStatement("INSERT INTO VCF(NAME) VALUES(?)",PreparedStatement.RETURN_GENERATED_KEYS);
					pstmt.setString(1, filename);
					if(pstmt.executeUpdate()!=1) {
						LOG.error("Cannot insert VCF ?");
						return -1;
					}
					final long vcf_id =getLastGeneratedId(pstmt);
					pstmt.close();
					
					pw.print(vcf_id);
					pw.print("\t");
					pw.println(filename);
					pw.flush();
					
					pstmt = this.conn.prepareStatement("SELECT ID FROM ROWCONTENT WHERE MD5SUM=?");
					pstmt2 = this.conn.prepareStatement("INSERT INTO ROWCONTENT(MD5SUM,CONTENT,CONTIG,START,STOP,ALLELE_REF,FILTERED) VALUES (?,?,?,?,?,?,?)",PreparedStatement.RETURN_GENERATED_KEYS);
					pstmt3 = this.conn.prepareStatement("INSERT INTO VCFROW(VCF_ID,ROW_ID) VALUES (?,?)");
					pstmt3.setLong(1, vcf_id);
					
					final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(cah.header);
					/* insert VCF header lines */
					for(final String line:headerLines) {
						final String md5=this.toMd5.apply(line);
						long content_id = -1L;
						pstmt.setString(1, md5);
						row = pstmt.executeQuery();
						while(row.next()) {
							content_id = row.getLong(1);
						}
						row.close();				
						
						/* vcf content was not found, create it */
						if(content_id==-1L) {
							pstmt2.setString(1, md5);
							
							pstmt2.setString(2,line);
							pstmt2.setNull(3,Types.VARCHAR);
							pstmt2.setNull(4,Types.INTEGER);
							pstmt2.setNull(5,Types.INTEGER);
							pstmt2.setNull(6,Types.VARCHAR);
							pstmt2.setShort(7, (short)1);
							if(pstmt2.executeUpdate()!=1) {
								LOG.error("Cannot insert ROWCONTENT ?");
								return -1;
							}
							content_id =getLastGeneratedId(pstmt2);
						}
						
						/* insert new VCF row */
						pstmt3.setLong(2, content_id);
						if(pstmt3.executeUpdate()!=1) {
							LOG.error("Cannot insert VCFROW ?");
							return -1;
						}
					}
					
					LOG.info("Inserted "+filename+" ID="+vcf_id);
					while(lineIter.hasNext() && !lineIter.peek().startsWith("#")) {
						final String line = lineIter.next();
						final String md5 = this.toMd5.apply(line);
						
						long content_id = -1L;
						pstmt.setString(1, md5);
						row = pstmt.executeQuery();
						while(row.next()) {
							content_id = row.getLong(1);
						}
						row.close();
						/* vcf variants content was not found, create it */
						if(content_id==-1L) {
							/* decode to get chrom/start/end/ref */
							final VariantContext ctx = progress.watch(cah.codec.decode(line));
							
							pstmt2.setString(1, md5);
							
							pstmt2.setString(2,line);
							pstmt2.setString(3, ctx.getContig());
							pstmt2.setInt(4, ctx.getStart());
							pstmt2.setInt(5, ctx.getEnd());
							String refBase =ctx.getReference().getBaseString();
							/* sql table for Ref_allele is a varchar(MAX_REF_BASE_LENGTH) */
							if(refBase.length()>50) {
								LOG.warn("Warning: TRUNCATING LARGE REF BASE TO FIT IN DATABASE : VARCHAR("+MAX_REF_BASE_LENGTH+") characters:"+refBase);
								refBase = refBase.substring(0,MAX_REF_BASE_LENGTH);
								++number_of_ref_allele_truncated;
							}
							pstmt2.setString(6,refBase );
							pstmt2.setShort(7, (short)(ctx.isFiltered()?1:0));
							if(pstmt2.executeUpdate()!=1) {
								LOG.error("Cannot insert ROWCONTENT ?");
								return -1;
							}
							content_id =getLastGeneratedId(pstmt2);
						}
						
						/* insert new VCF row */
						pstmt3.setLong(2, content_id);
						if(pstmt3.executeUpdate()!=1) {
							LOG.error("Cannot insert VCFROW ?");
							return -1;
						}
					}
					pstmt2.close();
					pstmt3.close();
					pstmt.close();
					progress.finish();
					num_vcf_in_this_stream++;
					} /* end of while iter has next */
				CloserUtil.close(lineIter);
				lineIter=null;
				fileidx++;
			} while(fileidx < args.size());
			
			pw.flush();
			pw.close();
			
			compress();
			LOG.warn("Number of REF alleles length(REF)> VARCHAR("+MAX_REF_BASE_LENGTH+") truncated:"+number_of_ref_allele_truncated);
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(pw);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			CloserUtil.close(pstmt2);
			CloserUtil.close(pstmt3);
			CloserUtil.close(lineIter);
		}
	}

			
	
	private int doCommandList(List<String> args){
		Statement pstmt = null;
		ResultSet row = null;
		PrintWriter pw = null;
		if(!args.isEmpty()) {
			LOG.error("Too many arguments");
			return -1;
		}
		try {
			/** user can add extra columns in the VCF, we collect the names of the columns */
			final List<String> cols =new ArrayList<>();
			pstmt = this.conn.createStatement();
			row = pstmt.executeQuery("SELECT  sys.SYSCOLUMNS.COLUMNNAME FROM  sys.SYSCOLUMNS, sys.SYSTABLES WHERE  sys.SYSTABLES.TABLEID= sys.SYSCOLUMNS.REFERENCEID AND  sys.SYSTABLES.TABLENAME=\'VCF\'");
			while(row.next()) {
				cols.add("VCF."+row.getString(1));
				}
			row.close();
			pstmt.close();
			final String sql=new StringBuilder("SELECT ").
					append(String.join(",", cols)).
					append(",COUNT(VCFROW.ID)  as \"COUNT_VARIANTS\" " 
						+ "FROM VCF,VCFROW,ROWCONTENT "
						+ "WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND ROWCONTENT.CONTIG IS NOT NULL "
						+ "GROUP BY "
						).
				append(String.join(",", cols)).
				toString()
				;
			
			LOG.info(sql);
			pstmt = this.conn.createStatement();
			row = pstmt.executeQuery(sql);
			final ResultSetMetaData meta = row.getMetaData();
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			for(int i=0;i< meta.getColumnCount();i++)
				{
				pw.print((i==0?"#":"\t"));
				pw.print(meta.getColumnLabel(i+1));
				}
			pw.println();
			while(row.next()) {
				for(int i=0;i< meta.getColumnCount();i++)
				{
				pw.print((i==0?"":"\t"));
				pw.print(row.getString(i+1));
				}
			pw.println();
			}
			
			pw.flush();
			return RETURN_OK;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(pw);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
		}
	}
	
	private int doCommandDelete(List<String> args){
		final Pattern comma = Pattern.compile("[,]");
		PreparedStatement pstmt1=null;
		PreparedStatement pstmt2=null;
		try {
			pstmt1 = this.conn.prepareStatement("DELETE FROM VCFROW WHERE VCF_ID=?");
			pstmt2 = this.conn.prepareStatement("DELETE FROM VCF WHERE ID=?");
			for(final String token: args) {
				for(final String idStr: comma.split(token))
					{
					if(idStr.trim().isEmpty()) continue;
					long vcf_id=-1L;
					try {
						vcf_id = Long.parseLong(idStr);
					} catch (NumberFormatException e) {
						vcf_id=-1L;
					}
					if(vcf_id<1L) {
						LOG.warn("Bad VCF ID :"+idStr);
						continue;
					}
					pstmt1.setLong(1, vcf_id);
					pstmt2.setLong(1, vcf_id);
					LOG.info("deleting vcf id:"+vcf_id);
					pstmt1.executeUpdate();
					pstmt2.executeUpdate();
					}
				}
			compress();
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(pstmt1);
			CloserUtil.close(pstmt2);
		}
	}

	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.derbyFilePath==null ||this.derbyFilePath.isEmpty())
				{
				LOG.error("Undefined Derby DB Directory");
				return -1;
				}
			try {
				Class.forName("org.apache.derby.jdbc.ClientDriver").newInstance();
			} catch(final Exception err ){
				LOG.error("Cannot get derby driver",err);
				return -1;
			}
			
			
			openDerby();
			final String command  = String.valueOf(this.actionStr);
			
			if(command.equals("read")) {
				return doReadConcatenatedVcf(args);
			} else if(command.equals("list")) {
				return doCommandList(args);
			}else if(command.equals("dump")) {
				return doCommandDump(args);
			}
			else if(command.equals("dumpall")) {
				return doCommandDumpAll(args);
			}
			else if(command.equals("dumpuniq")) {
				return doCommandDumpUniq(args);
			}
			else if(command.equals("delete")) {
				return doCommandDelete(args);
			}
			else {
				LOG.error("unknown command : "+command);
				return -1;
			}			
		} 
		catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			closeDerby();
		}
		}
	
	
	public static void main(String[] args)
		{
		new VcfDerby01().instanceMainWithExit(args);
		}
	}
