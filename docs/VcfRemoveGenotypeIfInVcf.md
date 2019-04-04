# VcfRemoveGenotypeIfInVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Reset Genotypes in VCF (./.) if they've been found in another VCF indexed with tabix


## Usage

```
Usage: vcfresetvcf [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -t, --tabix
      Tabix indexed VCF file
    --version
      print version and exit
    -x
      remove variant if there is no called genotype
      Default: false

```


## Keywords

 * vcf
 * genotype


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfresetvcf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveGenotypeIfInVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveGenotypeIfInVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfresetvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$ gunzip -c input.vcf.gz | grep -E '1308871'
1	1308871	.	A	T	6.2	.	.	GT:PL:DP:GQ	1/1:35,3,0:1:4



$ gunzip -c input.vcf.gz | grep -E '(^#|1308871)' |\
  java -jar dist/vcfresetvcf.jar -x ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz |\
  grep -v '^#'
1	1308871	.	A	T	6.20	.	.	GT	./.

```
