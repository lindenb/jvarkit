# VcfAlleleDepthFraction

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

filter VCF for strange FORMAT:AD fraction


## Usage

```
Usage: vcfadfraction [options] Files
  Options:
    -dp, --dp
      Only consider Genotypes having DP> 'x'
      Default: -1
    -filter, --filter
      Variant FILTER
      Default: AD_RATIO
    -f, --filtered
      ignore FILTER-ed **GENOTYPES**
      Default: false
    -gtf, --gtf
      Genotype FILTER
      Default: AD_RATIO
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -het, --het
      AD ratio for **HET** genotypes. HET genotype should have x <= 
      AD[1]/(AD[0]+AD[1])<= (1-x)
      Default: 0.2
    -hom, --hom
      AD ratio for **HOM_REF** or **HOM_VAR** genotypes. HOM_REF genotype 
      should have x <= AD[1]/(AD[0]+AD[1]). HOM_VAR genotype should have  
      AD[1]/(AD[0]+AD[1]) >= (1-x).
      Default: 0.05
    -maxFilteredGenotypes, --maxFilteredGenotypes
      Set Variant FILTER if number of BAD genotype is greater than 'x'. 
      Negative is ignore.
      Default: -1
    -maxFractionFilteredGenotypes, --maxFractionFilteredGenotypes
      Set Variant FILTER if percent of BAD genotype is greater than 'x'. 
      Negative is ignore.
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * allele-balance
 * depth


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfadfraction
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190723

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleDepthFraction.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleDepthFraction.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfadfraction** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ wget -O - -q "https://github.com/KennethJHan/Bioinformatics_Programming_101/raw/f40a5daa01cb9b37232d7d7a576100a10f867e80/GATK_BtPractice/SRR000982.filtered.variants.vcf" |  java -jar dist/vcfadfraction.jar | java -jar dist/vcf2table.jar | grep AD_RATIO -C 10

 Genotype Types
 +---------+-------+-----+
 | Type    | Count | %   |
 +---------+-------+-----+
 | HOM_VAR | 1     | 100 |
 +---------+-------+-----+
 Genotypes
 +-----------+---------+-----+----+----------+----+-----+-----+---------------+---------+
 | Sample    | Type    | AD  | DP | FT       | GQ | GT  | PGT | PID           | PL      |
 +-----------+---------+-----+----+----------+----+-----+-----+---------------+---------+
 | SRR000982 | HOM_VAR | 1,4 | 5  | AD_RATIO | 5  | 1/1 | 1|1 | 196625750_A_T | 163,5,0 |
 +-----------+---------+-----+----+----------+----+-----+-----+---------------+---------+
<<GRCh37 chr3:196625765/A (n. 90)
>>GRCh37 chr3:197119835/G (n. 91)
 Variant
 +-------+-----------+
 | Key   | Value     |
 +-------+-----------+
 | CHROM | chr3      |
 | POS   | 197119835 |
 | end   | 197119835 |
--
 Genotype Types
 +------+-------+-----+
 | Type | Count | %   |
 +------+-------+-----+
 | HET  | 1     | 100 |
 +------+-------+-----+
 Genotypes
 +-----------+------+------+----+----------+----+-----+----------+
 | Sample    | Type | AD   | DP | FT       | GQ | GT  | PL       |
 +-----------+------+------+----+----------+----+-----+----------+
 | SRR000982 | HET  | 3,16 | 19 | AD_RATIO | 56 | 0/1 | 472,0,56 |
 +-----------+------+------+----+----------+----+-----+----------+
<<GRCh37 chr10:42385236/A (n. 167)
>>GRCh37 chr10:42385255/G (n. 168)
```


