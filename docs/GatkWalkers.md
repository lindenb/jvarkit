Jvarkit contains some experimental Broad-Insitute GATK custom walkers ( https://software.broadinstitute.org/gatk/download/ ).

The walkers have been tested with **GATK 3.7**

# Compiling

* javac 1.8 and GNU Make are required.
* Download & extract GATK https://software.broadinstitute.org/gatk/download/
* In the root of jvarkit, create a file `local.mk` with a property `gatk.jar` to the full path of **GenomeAnalysisTK.jar** . e.g: 

```makefile
gatk.jar=/dir1/dir2/GenomeAnalysisTK.jar
```

* run the Makefile to compile the custom walkers:

```
$ make gatkwalkers
```

this should create a `mygatk.jar` in the `dist` directory

# Using the Walkers

The GATK is invoked as follow:

```bash
$ java -cp /dir1/dir2/GenomeAnalysisTK.jar:/path/to/mygatk.jar org.broadinstitute.gatk.engine.CommandLineGATK -T ...
```

# Walkers

## EigenVariants

Variant Annotator for the data of https://xioniti01.u.hpc.mssm.edu/v1.1/ : Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance.           

```
Arguments for EigenVariants:
 -eigen,--eigenDirectory <eigenDirectory>   The Eigen directory containing the tabix indexed files  *.tab.gz.
 -V,--variant <variant>                     Input VCF file
 -tabixPrefix,--tabixPrefix <tabixPrefix>   Override prefix of tabix file in the eigen directory. Leave null for 
                                            default.
 -o,--out <out>                             File to which variants should be written
```


## GroupByGenotypes              

Reads a VCF file and creates a genotype summary table

```
Arguments for GroupByGenotypes:
 -V,--variant <variant>                           Input VCF file
 -mgq,--minGenotypeQuality <minGenotypeQuality>   Minimum genotype quality to put a assign a category
 -chrom,--chrom                                   Group by Chromosome/Contig
 -ID,--ID                                         Group by ID
 -variantType,--variantType                       Group by VariantType
 -genotypeType,--genotypeType                     Group by GenotypeType
 -filter,--filter                                 Group by FILTER
 -gfilter,--gfilter                               Group by GENOTYPE FILTER
 -impact,--impact                                 Group by ANN/IMPACT
 -onlysingletons,--onlysingletons                 only consider singletons (one sample affected/variant)
 -o,--out <out>                                   File to which result should be written
```

## GroupByVariants               

Reads a VCF file and creates a variant summary table


```
Arguments for GroupByVariants:
 -V,--variant <variant>                  Input VCF file
 -mq,--minQuality <minQuality>           Group by Quality. Set the treshold for Minimum Quality
 -chrom,--chrom                          Group by Chromosome/Contig
 -ID,--ID                                Group by ID
 -variantType,--variantType              Group by VariantType
 -filter,--filter                        Group by FILTER
 -impact,--impact                        Group by ANN/IMPACT
 -biotype,--biotype                      Group by ANN/biotype
 -nalts,--nalts                          Group by number of ALTS
 -affected,--affected                    Group by number of Samples called and not HOMREF
 -called,--called                        Group by number of Samples called
 -maxSamples,--maxSamples <maxSamples>   if the number of samples affected is greater than --maxSamples use the label 
                                         "GT_MAX_SAMPLES"
 -tsv,--tsv                              Group by Transition/Transversion
 -allelesize,--allelesize                Group by Max(allele.size)
 -o,--out <out>                          File to which result should be written
```

## WindowVariants   

Annotate Variants using a sliding window

```
Arguments for WindowVariants:
 -V,--variant <variant>                            Input VCF file
 -select,--selectexpressions <selectexpressions>   Optional Jexl expression to use when selecting the data
 -shift,--windowShift <windowShift>                Window shift (in bp.)
 -wsize,--windowSize <windowSize>                  Window Size (in bp.)
 -wname,--windowName <windowName>                  INFO Attribute name that will be added
 -o,--out <out>                                    File to which variants should be written
```


e.g:

```
 (...) -T WindowVariants -V input.vcf -R ref.fasta -select 'vc.hasID()'

##INFO=<ID=WINDOW,Number=.,Type=String,Description="Window : start|end|number-of-matching-variants|number-of-non-matching-variants">
(...)
2	35565407	.	G	A	.	PASS	WINDOW=35565400|35565550|0|1,35565350|35565500|2|2,35565300|35565450|2|2
2	35565471	rs11234	C	T	.	PASS	WINDOW=35565450|35565600|1|0
2	35565628	.	T	G	.	PASS	WINDOW=35565600|35565750|1|2

```

