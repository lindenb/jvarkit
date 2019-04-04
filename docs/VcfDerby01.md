# VcfDerby01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Insert similar VCFs into an Apache Derby Database


## Usage

```
Usage: vcfderby01 [options] Files
  Options:
    -a, --action
      REQUIRED. action to perform. 'read': read a zip or a concatenated stream 
      of vcf files and insert it into a derby database. 'list': list the 
      available vcf. 'dump' dump one or more VCF. 'dumpall' dump all VCFs. 
      'dumpuniq' dum all as a one and only uniq vcf. 'delete' : delete one or 
      more VCF by ID
      Default: <empty string>
    -d, --derby
      REQUIRED. path to Derby database storage directory.
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -t, --title
      Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the 
      name of the inserted VCF file
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * vcf
 * sql
 * derby
 * burden


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfderby01
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfDerby01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfDerby01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfderby01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





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




