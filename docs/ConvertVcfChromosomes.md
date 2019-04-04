# ConvertVcfChromosomes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert the names of the chromosomes in a VCF file


## DEPRECATED

use `bcftools annotate` with `--rename-chrs file`

## Usage

```
Usage: vcfrenamechr [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -f, --mapping, -m
      load a custom name mapping. Format (chrom-source\tchrom-dest\n)+
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -c, -convert
      What should I do when  a converstion is not found
      Default: RAISE_EXCEPTION
      Possible Values: [RAISE_EXCEPTION, SKIP, RETURN_ORIGINAL]

```


## Keywords

 * vcf
 * contig
 * chromosome
 * convert


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfrenamechr
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertVcfChromosomes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertVcfChromosomes.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/ConvertVcfChromosomesTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/ConvertVcfChromosomesTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfrenamechr** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Deprecated

use `bcftools annotate` with `--rename-chrs file`

> rename chromosomes according to the map in file, with "old_name new_name\n" pairs separated by whitespaces, each on a separate line. 

## Example

```bash
$ curl  "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  java -jar dist/vcfrenamechr.jar -i -C -f src/main/resources/chromnames/hg19_to_g1kv37.tsv |\
cut -f 1-5

(...)
##contig=<ID=1,length=249250621>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=GL000202.1,length=40103>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=GL000203.1,length=37498>
##contig=<ID=GL000204.1,length=81310>
##contig=<ID=GL000205.1,length=174588>
##contig=<ID=GL000206.1,length=41001>
##contig=<ID=18,length=78077248>
##contig=<ID=GL000207.1,length=4262>
##contig=<ID=19,length=59128983>
##contig=<ID=GL000208.1,length=92689>
##contig=<ID=GL000209.1,length=159169>
##contig=<ID=GL000191.1,length=106433>
##contig=<ID=GL000192.1,length=547496>
##contig=<ID=2,length=243199373>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=GL000210.1,length=27682>
##contig=<ID=22,length=51304566>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=GL000193.1,length=189789>
##contig=<ID=GL000194.1,length=191469>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=GL000195.1,length=182896>
##contig=<ID=8,length=146364022>
##contig=<ID=GL000196.1,length=38914>
##contig=<ID=GL000197.1,length=37175>
##contig=<ID=9,length=141213431>
##contig=<ID=GL000198.1,length=90085>
##contig=<ID=GL000199.1,length=169874>
##contig=<ID=GL000200.1,length=187035>
##contig=<ID=GL000201.1,length=36148>
##contig=<ID=GL000211.1,length=166566>
##contig=<ID=GL000212.1,length=186858>
##contig=<ID=GL000213.1,length=164239>
##contig=<ID=GL000214.1,length=137718>
##contig=<ID=GL000215.1,length=172545>
##contig=<ID=GL000216.1,length=172294>
##contig=<ID=GL000217.1,length=172149>
##contig=<ID=GL000218.1,length=161147>
##contig=<ID=GL000219.1,length=179198>
##contig=<ID=GL000220.1,length=161802>
##contig=<ID=GL000221.1,length=155397>
##contig=<ID=GL000222.1,length=186861>
##contig=<ID=GL000223.1,length=180455>
##contig=<ID=GL000224.1,length=179693>
##contig=<ID=GL000225.1,length=211173>
##contig=<ID=GL000226.1,length=15008>
##contig=<ID=GL000227.1,length=128374>
##contig=<ID=GL000228.1,length=129120>
##contig=<ID=GL000229.1,length=19913>
##contig=<ID=GL000230.1,length=43691>
##contig=<ID=GL000231.1,length=27386>
##contig=<ID=GL000232.1,length=40652>
##contig=<ID=GL000233.1,length=45941>
##contig=<ID=GL000234.1,length=40531>
##contig=<ID=GL000235.1,length=34474>
##contig=<ID=GL000236.1,length=41934>
##contig=<ID=GL000237.1,length=45867>
##contig=<ID=GL000238.1,length=39939>
##contig=<ID=GL000239.1,length=33824>
##contig=<ID=GL000240.1,length=41933>
##contig=<ID=GL000241.1,length=42152>
##contig=<ID=GL000242.1,length=43523>
##contig=<ID=GL000243.1,length=43341>
##contig=<ID=GL000244.1,length=39929>
##contig=<ID=GL000245.1,length=36651>
##contig=<ID=GL000246.1,length=38154>
##contig=<ID=GL000247.1,length=36422>
##contig=<ID=GL000248.1,length=39786>
##contig=<ID=GL000249.1,length=38502>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##reference=file:///m/cphg-quinlan/cphg-quinlan/shared/genomes/hg19/bwa/gatk/hg19_gatk.fa
#CHROM	POS	ID	REF	ALT
1	145273345	.	T	C
1	156011444	.	T	C
5	64982321	.	T	C
10	1142208	.	T	C
10	126678092	.	G	A
10	135210791	.	T	C
13	48873835	.	G	A
20	36779424	.	G	A
X	17819377	.	T	C


```

## See also

* https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/g1kv37_to_hg19.tsv
* https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/hg19_to_g1kv37.tsv
* [[VcfSampleRename]]
* [[BamRenameChromosomes]]
* http://plindenbaum.blogspot.fr/2013/07/g1kv37-vs-hg19.html


