# HowManyBamDict

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

finds if there's are some differences in the sequence dictionaries.


## Usage

```
Usage: howmanybamdict [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * dict
 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/468541](https://www.biostars.org/p/468541)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew howmanybamdict
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20131108

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/HowManyBamDict.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/HowManyBamDict.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **howmanybamdict** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
t$ find /home/lindenb/src/jvarkit/src/test/resources -name "*.vcf.gz" -o -name "*.bam" -o -name "*.cram" | java -jar dist/howmanybamdict.jar | cut -c 1-${COLUMNS}
9a5c58c2c91e731135b27ed14974523a	.	85	3101976562	1=249250621;2=243199373;3=198022430;4=191154276;5=180915260;6=17111
9a5c58c2c91e731135b27ed14974523a	/home/lindenb/src/jvarkit/src/test/resources/ExAC.r1.sites.vep.vcf.gz
4677ece43eea2b029d0d33fe130ea6c7	.	86	3137454505	chr1=249250621;chr2=243199373;chr3=198022430;chr4=191154276;chr5=18
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.B00I9CJ.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	.	11	18490	RF01=3302;RF02=2687;RF03=2592;RF04=2362;RF05=1579;RF06=1356;RF07=1074;RF
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S2.vcf.gz
9a5c58c2c91e731135b27ed14974523a	/home/lindenb/src/jvarkit/src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S3.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S1.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S4.vcf.gz
635de5cb51973d45844fa713ac0b7719	.	2	85	ref=45;ref2=40	/home/lindenb/src/jvarkit/src/test/resources/toy.vcf.gz
635de5cb51973d45844fa713ac0b7719	/home/lindenb/src/jvarkit/src/test/resources/toy.vcf.gz
635de5cb51973d45844fa713ac0b7719	/home/lindenb/src/jvarkit/src/test/resources/toy.bam
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S4.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/roxan.hs37d5.csq.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.ann.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.vcf.gz
df8200dd2a49e25bc98df5f2c45ac36a	.	2	85	ref=45;ref2=40	/home/lindenb/src/jvarkit/src/test/resources/toy.cram
df8200dd2a49e25bc98df5f2c45ac36a	/home/lindenb/src/jvarkit/src/test/resources/toy.cram
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.freebayes.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S1.bam
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S2.bam
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S5.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.B00GWIU.vcf.gz
5de1aef8e9e3be34ae7696779a4791c7	.	3366	3217346917	chr1=248956422;chr2=242193529;chr3=198295559;chr4=190214555;chr5=
5de1aef8e9e3be34ae7696779a4791c7	/home/lindenb/src/jvarkit/src/test/resources/FAB23716.nanopore.bam
e07f3b093938833945fa357c4b37bdf9	.	292	3100014256	chr1=248956422;chr2=242193529;chr3=198295559;chr4=190214555;chr5=1
e07f3b093938833945fa357c4b37bdf9	/home/lindenb/src/jvarkit/src/test/resources/ENCFF331CGL.rnaseq.b38.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/retrocopy01.bwa.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.B00GWGD.vcf.gz
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.D000Q1R.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S5.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S3.bam
f8d942cb3fc6ebef618a0b0ba3f4ef99	.	24	3095677412	1=249250621;2=243199373;3=198022430;4=191154276;5=180915260;6=17111
f8d942cb3fc6ebef618a0b0ba3f4ef99	/home/lindenb/src/jvarkit/src/test/resources/gnomad_v2_sv.sites.vcf.gz
4327878bf3a3073edf6e77fda48033fe	.	86	3137454505	1=249250621;2=243199373;3=198022430;4=191154276;5=180915260;6=17111
4327878bf3a3073edf6e77fda48033fe	/home/lindenb/src/jvarkit/src/test/resources/HG02260.transloc.chr9.14.bam
9a5c58c2c91e731135b27ed14974523a	/home/lindenb/src/jvarkit/src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz
(...)
```

