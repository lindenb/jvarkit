# VcfGnomadPext

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Peek annotations from gnomadpext


## Usage

```
Usage: vcfgnomadpext [options] Files
  Options:
    --bufferSize
      When we're looking for variant in Gnomad, load the variants for 'N' 
      bases instead of doing a random access for each variant. A distance 
      specified as a positive integer.Comma are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 10000
  * -d, --database, --pext
      Pext file. tab delimited :(chrom\tpos\tref\talt\ttx_annotation). Bgziped 
      and indexed with tabix.
    -filtered, --filtered
      Skip Filtered User Variants
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tissues
      Restrict to those tissues.
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * gnomad


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgnomadpext
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190220

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomadpext/VcfGnomadPext.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomadpext/VcfGnomadPext.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomadpext** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Context

https://www.biorxiv.org/content/10.1101/554444v1

> "Transcript expression-aware annotation improves rare variant discovery and interpretation" 
> Here, we develop a transcript-level annotation metric, the proportion expressed across transcripts (pext), which summarizes isoform quantifications for variants. We calculate this metric using 11,706 tissue samples from the Genotype Tissue Expression project (GTEx) and show that it clearly differentiates between weakly and highly evolutionarily conserved exons, a proxy for functional importance. 

## Example

```
# bgzip if needed and index the database with tabix
$ bgzip data.tsv
$ tabix -f -b 2 -e 2 -s 1 -c 'c' data.tsv.gz

# annotate
$ java -jar dist/vcfgnomadpext.jar -d data.tsv.gz input.vcf.gz

(...)
1	135803	.	C	G	.	.	AC=8;AF=0.500;AN=16;DP=113;GNOMAD_PEXT=ensg:ENSG00000237683|csq:splice_acceptor_variant|symbol:AL627309.1|lof:LC|lof_flag:NON_CAN_SPLICE|Spleen:1.0|Brain_FrontalCortex_BA9_:1.0|SmallIntestine_TerminalIleum:1.0|Skin_SunExposed_Lowerleg_:1.0|Artery_Coronary:1.0|Brain_Hippocampus:1.0|Esophagus_Muscularis:1.0|Brain_Nucleusaccumbens_basalganglia_:1.0|Artery_Tibial:1.0|Brain_Hypothalamus:1.0|Adipose_Visceral_Omentum_:1.0|Nerve_Tibial:1.0|Brain_CerebellarHemisphere:1.0|Breast_MammaryTissue:1.0|Liver:1.0|Skin_NotSunExposed_Suprapubic_:1.0|AdrenalGland:1.0|Pancreas:1.0|Lung:1.0|Pituitary:1.0|Muscle_Skeletal:1.0|Colon_Transverse:1.0|Artery_Aorta:1.0|Heart_AtrialAppendage:1.0|Adipose_Subcutaneous:1.0|Esophagus_Mucosa:1.0|Heart_LeftVentricle:1.0|Brain_Cerebellum:1.0|Brain_Cortex:1.0|Thyroid:1.0|Stomach:1.0|WholeBlood:1.0|Brain_Anteriorcingulatecortex_BA24_:1.0|Brain_Putamen_basalganglia_:1.0|Brain_Caudate_basalganglia_:1.0|Colon_Sigmoid:1.0|Esophagus_GastroesophagealJunction:1.0|Brain_Amygdala:1.0|mean_proportion:1.0	GT:DP	1/10/0:39	1/1	1/0:25	0/0:11	1/0	0/1:38	0/1
1	135804	.	G	A	.	.	AC=6;AF=0.375;AN=16;DP=82;GNOMAD_PEXT=ensg:ENSG00000237683|csq:splice_acceptor_variant|symbol:AL627309.1|lof:LC|lof_flag:NON_CAN_SPLICE|Spleen:1.0|Brain_FrontalCortex_BA9_:1.0|SmallIntestine_TerminalIleum:1.0|Skin_SunExposed_Lowerleg_:1.0|Artery_Coronary:1.0|Brain_Hippocampus:1.0|Esophagus_Muscularis:1.0|Brain_Nucleusaccumbens_basalganglia_:1.0|Artery_Tibial:1.0|Brain_Hypothalamus:1.0|Adipose_Visceral_Omentum_:1.0|Nerve_Tibial:1.0|Brain_CerebellarHemisphere:1.0|Breast_MammaryTissue:1.0|Liver:1.0|Skin_NotSunExposed_Suprapubic_:1.0|AdrenalGland:1.0|Pancreas:1.0|Lung:1.0|Pituitary:1.0|Muscle_Skeletal:1.0|Colon_Transverse:1.0|Artery_Aorta:1.0|Heart_AtrialAppendage:1.0|Adipose_Subcutaneous:1.0|Esophagus_Mucosa:1.0|Heart_LeftVentricle:1.0|Brain_Cerebellum:1.0|Brain_Cortex:1.0|Thyroid:1.0|Stomach:1.0|WholeBlood:1.0|Brain_Anteriorcingulatecortex_BA24_:1.0|Brain_Putamen_basalganglia_:1.0|Brain_Caudate_basalganglia_:1.0|Colon_Sigmoid:1.0|Esophagus_GastroesophagealJunction:1.0|Brain_Amygdala:1.0|mean_proportion:1.0	GT:DP	0/00/0	0/1:15	0/0:5	1/1:15	0/1:17	0/0:8	1/1:22
1	137619	.	G	A	.	.	AC=7;AF=0.438;AN=16;DP=193;GNOMAD_PEXT=ensg:ENSG00000237683|csq:splice_donor_variant|symbol:AL627309.1|lof:LC|lof_flag:NON_CAN_SPLICE|Spleen:1.0|Brain_FrontalCortex_BA9_:1.0|SmallIntestine_TerminalIleum:1.0|Skin_SunExposed_Lowerleg_:1.0|Artery_Coronary:1.0|Brain_Hippocampus:1.0|Esophagus_Muscularis:1.0|Brain_Nucleusaccumbens_basalganglia_:1.0|Artery_Tibial:1.0|Brain_Hypothalamus:1.0|Adipose_Visceral_Omentum_:1.0|Nerve_Tibial:1.0|Brain_CerebellarHemisphere:1.0|Breast_MammaryTissue:1.0|Liver:1.0|Skin_NotSunExposed_Suprapubic_:1.0|AdrenalGland:1.0|Pancreas:1.0|Lung:1.0|Pituitary:1.0|Muscle_Skeletal:1.0|Colon_Transverse:1.0|Artery_Aorta:1.0|Heart_AtrialAppendage:1.0|Adipose_Subcutaneous:1.0|Esophagus_Mucosa:1.0|Heart_LeftVentricle:1.0|Brain_Cerebellum:1.0|Brain_Cortex:1.0|Thyroid:1.0|Stomach:1.0|WholeBlood:1.0|Brain_Anteriorcingulatecortex_BA24_:1.0|Brain_Putamen_basalganglia_:1.0|Brain_Caudate_basalganglia_:1.0|Colon_Sigmoid:1.0|Esophagus_GastroesophagealJunction:1.0|Brain_Amygdala:1.0|mean_proportion:1.0	GT:DP	0/1:25	0/1:10	0/1:46	0/0:32	0/0	1/0	1/1:31	0/1:49
```



