# VcfRemoveGenotypeJs

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Reset Genotype in VCF using a javascript expression


## Usage

```
Usage: vcfremovegenotypejs [options] Files
  Options:
    -e, --expression
       (js expression). Optional.
    -filter, --filter
      if not empty, don't delete the genotype but filter it.
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -homref, --homref
      Replace variant with homref instead of nocall
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -g, --removeCtxNoGenotype
      Remove variants having no called genotype or all home Ref.
      Default: false
    -f, --script
       (js file). Optional.
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfremovegenotypejs
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveGenotypeJs.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveGenotypeJs.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfremovegenotypejs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Example

The script injects in the context:
 *  header a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html
 *  variant a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
 *  genotype a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html

if the returned value is false, the genotype is set to no call or to hom-ref.


```
$ cat  ~/src/gatk-ui/testdata/mutations.vcf | java -jar dist/vcfresetgenotypejs.jar -homref -e '!genotype.isHomRef()' |grep -v "##"

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  S1      S2      S3      S4
rotavirus       51      .       A       G       22.55   .       AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SGB=-61.9012;VDB=3.53678e-05  GT:PL   0/0:0,255,134   0/0:0,255,127   0/0:0,255,137   1/1:70,255,0
rotavirus       91      .       A       T       5.45    .       AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SGB=-369.163;VDB=0.313337 GT:PL   0/0:0,255,133   0/1:40,0,31     0/0:0,255,134   0/0:0,255,82
rotavirus       130     .       T       C       4.12    .       AC1=1;AF1=0.124933;BQB=1;DP=1349;DP4=1139,0,204,0;FQ=4.48321;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.762964;SGB=-335.275;VDB=0.00084636  GT:PL0/1:38,0,35      0/0:0,255,132   0/0:0,255,132   0/0:0,255,79


```


 */
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
 
 The script injects in the context:

  * **header** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html)
  * **variant** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html)
  * **genotype** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html)



if the returned value is false, the genotype is set to no call or to hom-ref.

```bash
$ cat  ~/src/gatk-ui/testdata/mutations.vcf | java -jar dist/vcfresetgenotypejs.jar -homref -e '!genotype.isHomRef()' |grep -v "##"

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SGB=-61.9012;VDB=3.53678e-05	GT:PL	0/0:0,255,134	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SGB=-369.163;VDB=0.313337	GT:PL	0/0:0,255,133	0/1:40,0,31	0/0:0,255,134	0/0:0,255,82
rotavirus	130	.	T	C	4.12	.	AC1=1;AF1=0.124933;BQB=1;DP=1349;DP4=1139,0,204,0;FQ=4.48321;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.762964;SGB=-335.275;VDB=0.00084636	GT:PL	0/1:38,0,35	0/0:0,255,132	0/0:0,255,132	0/0:0,255,79
```



 
