# FindGVCFsBlocks

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find common blocks of calleable regions from a set of gvcfs


## Usage

```
Usage: findgvcfsblocks [options] Files
  Options:
    --min-size, --block-size
      min block size. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 0
    -c, --chrom, --chromosome, --contig
      limit to that contig
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -T
      temporary directory

```


## Keywords

 * gvcf
 * gatk
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew findgvcfsblocks
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210806

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gvcf/FindGVCFsBlocks.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gvcf/FindGVCFsBlocks.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findgvcfsblocks** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

find regions for running GATK CombineGVCFs in parallel.

## Input

input is a set of path to the indexed g.vcf files or a picard-style interval file generated with a previous invocation of findgvcfsblocks.jar with one sample.
or it's a file with the '.list' suffix containing the path to the g.vcf files/interval files

g.vcf files must be indexed if option `-c` is used.

## Output

output is a picard-style **Interval** file containing the calleable GVCFs blocks.

## Example

```
$ java -jar dist/findgvcfsblocks.jar --min-size 100 --chrom RF11 S1.g.vcf.gz S2.g.vcf.gz S3.g.vcf.gz 
@HD	VN:1.6
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@CO	findgvcfsblocks. compilation:20210807160340 githash:b442941 htsjdk:2.24.1 date:20210807160354. cmd:--min-size 100 --chrom RF11 S1.g.vcf.gz S2.g.vcf.gz S3.g.vcf.gz
RF11	1	95	+	.
RF11	96	182	+	.
RF11	183	237	+	.
RF11	238	428	+	.
RF11	429	528	+	.
RF11	529	628	+	.
RF11	629	666	+	.
(...)
```

## Nextflow

```
(...)
Channel.fromPath(params.reference+".fai").
	splitCsv(header: false,sep:'\t',strip:true).
	filter{T->T[0].matches("(chr)?[0-9XY]+")}.
	map{T->[T[0]]}.
	set{each_contig}

process gvcflists {
executor "local"
output:
	path("gvcfs.list") into (gvcfs_list1,gvcfs_list2)
script:
"""
test -s "${params.gvcfs}"

SQRT=`awk 'END{z=sqrt(NR); print (z==int(z)?z:int(z)+1);}' "${params.gvcfs}"`

split -a 9 --additional-suffix=.list --lines=\${SQRT} "${params.gvcfs}" chunck.

find \${PWD} -type f -name "chunck.*.list"  > gvcfs.list

"""
}


gvcfs_list2.splitCsv(header: false,sep:'\t',strip:true).
	map{T->T[0]}.
	set{one_gvcf_list}


process findBlocks {
tag "${file(vcfs).name} ${contig}"
memory "5g"
input:
	tuple vcfs,contig from one_gvcf_list.combine(each_contig)
output:
	tuple contig,path("block.interval_list.gz") into blocks_interval
script:
"""
module load jvarkit/0.0.0

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/findgvcfsblocks.jar \
	--contig "${contig}" \
	-T . \
	-o "block.interval_list.gz" \
	"${vcfs}"

test -s block.interval_list.gz
"""
}


process findCommonBlocks {
tag "${contig} N=${L.size()}"
memory "10g"
input:
	tuple contig,L from blocks_interval.groupTuple()
output:
	path("block.interval_list") into blocks_common
script:
"""
module load jvarkit/0.0.0

cat << EOF > jeter.list
${L.join("\n")}
EOF

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/findgvcfsblocks.jar \
	--contig "${contig}" \
	-T . \
	--block-size "${params.blocksize}" \
	-o "block.interval_list" \
	jeter.list

rm jeter.list

test -s block.interval_list
"""
}


blocks_common.splitCsv(header: false,sep:'\t',strip:true).
	filter{T->!T[0].startsWith("@")}.
	map{T->T[0]+":"+T[1]+"-"+T[2]}.
	set{intervals}

process genotypeRegion {
tag "${region}"
cache "lenient"
errorStrategy "finish"
afterScript "rm -rf TMP"
memory "15g"
input:
	val region from intervals
	path gvcfs from gvcfs_list1
output:
	path("genotyped.vcf.gz")  into genotyped
	path("genotyped.vcf.gz.tbi")
script:
(...)

```



