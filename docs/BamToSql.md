# BamToSql

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert a SAM/BAM to sqlite statements


## Usage

```
Usage: bam2sql [options] Files
  Options:
    -c, --cigar
      print cigar data
      Default: false
    -f, --flag
      expands details about sam flag
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -r, --region
      An interval as the following syntax : "chrom:start-end" or 
      "chrom:middle+extend"  or "chrom:start-end+extend" or 
      "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * sql
 * sqlite


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2sql
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BamToSql.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BamToSql.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2sql** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Motivation

Inserting a BAM in a SQL is not a good idea of course !

But it might be interesting to get some informations about the bases in a segment of bam.

## Schema
The schema can change if some options (-c , -f) are used.
At the time of writing the schema is :

```sql
CREATE TABLE IF NOT EXISTS SamFile
(
id INTEGER PRIMARY KEY,
filename TEXT
);

CREATE TABLE IF NOT EXISTS Dictionary
(
id INTEGER PRIMARY KEY,
name TEXT NOT NULL,
length INT NOT NULL,
tid INT NOT NULL,
samfile_id INT NOT NULL,
FOREIGN KEY(samfile_id) REFERENCES SamFile(id)
);

CREATE TABLE IF NOT EXISTS ReadGroup
(
id INTEGER PRIMARY KEY,
groupId TEXT NOT NULL,
sample TEXT NOT NULL,
samfile_id INT NOT NULL,
FOREIGN KEY(samfile_id) REFERENCES SamFile(id)
);

CREATE TABLE IF NOT EXISTS Read
(
id INTEGER PRIMARY KEY,
name TEXT NOT NULL,
flag INTEGER NOT NULL,
rname TEXT,
pos INTEGER,
mapq INTEGER NOT NULL,
cigar TEXT,
rnext TEXT,
pnext INTEGER,
tlen INTEGER,
sequence TEXT NOT NULL,
qualities TEXT NOT NULL,
samfile_id INT NOT NULL,
group_id INT,
FOREIGN KEY(samfile_id) REFERENCES SamFile(id),
FOREIGN KEY(group_id) REFERENCES ReadGroup(id)
);

CREATE TABLE IF NOT EXISTS Cigar
(
id INTEGER PRIMARY KEY,
read_pos INT ,
read_base TEXT,
read_qual INT ,
ref_pos INT ,
ref_base TEXT,
operator TEXT NOT NULL,
read_id INT NOT NULL,
FOREIGN KEY(read_id) REFERENCES Read(id)
);
```

## Example

Build a sqlite3 database for a set of BAM files in the region "rotavirus:1-10""

```
$java -jar dist/bam2sql.jar -r 'rotavirus:1-10' -R  ref.fa -c S*.bam |\
sqlite3 database.sqlite
```

Select data from sqlite database where the genomic position is "rotavirus:5" 

```sql
select  SamFile.filename,
		ReadGroup.sample,
		Read.flag,
		Read.rname,
		Cigar.operator,
		Cigar.read_pos,
		Cigar.read_base,
		Cigar.read_qual,
		Cigar.ref_pos,
		Cigar.ref_base
from
		SamFile,Read,Cigar,ReadGroup
where
		SamFile.id = Read.samfile_id AND
		ReadGroup.id = Read.group_id AND 
		Cigar.read_id = Read.id and
		Read.rname = "rotavirus" and 
		Cigar.ref_pos= 5
		;
```

query:

```
$ sqlite3 -header -separator '   ' database.sqlite &lt; query.sql  | column -t 
```

output:

```
filename  sample  flag  rname      operator  read_pos  read_base  read_qual  ref_pos  ref_base
S1.bam    S1      99    rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      99    rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         3         T          10         5        T
S1.bam    S1      163   rotavirus  M         3         T          10         5        T
S1.bam    S1      99    rotavirus  M         2         T          10         5        T
S1.bam    S1      99    rotavirus  M         2         T          10         5        T
S1.bam    S1      163   rotavirus  M         2         T          10         5        T
S1.bam    S1      163   rotavirus  M         2         T          10         5        T
S1.bam    S1      99    rotavirus  M         1         T          10         5        T
S1.bam    S1      99    rotavirus  M         1         T          10         5        T
S1.bam    S1      163   rotavirus  M         1         T          10         5        T
S1.bam    S1      163   rotavirus  M         1         T          10         5        T
S1.bam    S1      163   rotavirus  M         4         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      163   rotavirus  M         0         T          10         5        T
S1.bam    S1      99    rotavirus  S         3         T          10         5        T
S1.bam    S1      163   rotavirus  S         0         T          10         5        T
S1.bam    S1      99    rotavirus  S         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         4         T          10         5        T
S2.bam    S2      163   rotavirus  M         4         T          10         5        T
S2.bam    S2      99    rotavirus  M         4         T          10         5        T
S2.bam    S2      99    rotavirus  M         4         T          10         5        T
S2.bam    S2      163   rotavirus  M         4         T          10         5        T
S2.bam    S2      163   rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         2         A          10         5        T
S2.bam    S2      163   rotavirus  M         2         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      163   rotavirus  M         3         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      99    rotavirus  M         1         T          10         5        T
S2.bam    S2      99    rotavirus  M         0         T          10         5        T
S2.bam    S2      99    rotavirus  M         0         T          10         5        T
S2.bam    S2      163   rotavirus  S         4         T          10         5        T
S2.bam    S2      99    rotavirus  S         2         A          10         5        T
S3.bam    S3      99    rotavirus  M         4         A          10         5        T
S3.bam    S3      163   rotavirus  M         4         T          10         5        T
S3.bam    S3      99    rotavirus  M         4         T          10         5        T
S3.bam    S3      99    rotavirus  M         3         T          10         5        T
S3.bam    S3      99    rotavirus  M         3         T          10         5        T
S3.bam    S3      99    rotavirus  M         3         T          10         5        T
S3.bam    S3      163   rotavirus  M         3         T          10         5        T
S3.bam    S3      163   rotavirus  M         3         T          10         5        T
S3.bam    S3      163   rotavirus  M         2         T          10         5        T
S3.bam    S3      163   rotavirus  M         2         T          10         5        T
S3.bam    S3      99    rotavirus  M         2         T          10         5        T
S3.bam    S3      163   rotavirus  M         2         T          10         5        T
S3.bam    S3      99    rotavirus  M         1         T          10         5        T
S3.bam    S3      163   rotavirus  M         1         A          10         5        T
S3.bam    S3      99    rotavirus  M         1         A          10         5        T
S3.bam    S3      99    rotavirus  M         1         A          10         5        T
S3.bam    S3      99    rotavirus  M         1         T          10         5        T
S3.bam    S3      163   rotavirus  M         1         T          10         5        T
S3.bam    S3      99    rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         A          10         5        T
S3.bam    S3      99    rotavirus  M         0         T          10         5        T
S3.bam    S3      163   rotavirus  M         0         T          10         5        T
S3.bam    S3      99    rotavirus  S         2         A          10         5        T
S4.bam    S4      163   rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         4         T          10         5        T
S4.bam    S4      99    rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         3         T          10         5        T
S4.bam    S4      163   rotavirus  M         3         T          10         5        T
S4.bam    S4      99    rotavirus  M         3         T          10         5        T
S4.bam    S4      163   rotavirus  M         2         T          10         5        T
S4.bam    S4      99    rotavirus  M         1         T          10         5        T
S4.bam    S4      99    rotavirus  M         0         T          10         5        T
S4.bam    S4      99    rotavirus  M         4         T          10         5        T
S4.bam    S4      163   rotavirus  M         0         T          10         5        T
S4.bam    S4      163   rotavirus  M         0         T          10         5        T
```

 
