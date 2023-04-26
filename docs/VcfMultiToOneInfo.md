# VcfMultiToOneInfo

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

'one variant with INFO with N values' to 'N variants with one INFO'


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfmulti2oneinfo  [options] Files

Usage: vcfmulti2oneinfo [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -i, --info
      The INFO tag
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf



## Creation Date

20260106

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmulti2oneinfo/VcfMultiToOneInfo.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmulti2oneinfo/VcfMultiToOneInfo.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmulti2oneinfo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


before:

```
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|frameshift_variant&start_lost|HIGH|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|p.Met1fs||1/297|1/98||,CAA|disruptive_inframe_deletion|MODERATE|Gene_20_616|Gene_20_616|transcript|AAG15311.1|protein_coding|1/1|c.57_59delAAA|p.Lys19del|57/597|57/597|19/198||,CAA|upstream_gene_variant|MODIFIER|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|||||2|;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
```

run:

```
java -jar dist/jvarkit.jar  vcfmulti2oneinfo -i ANN src/test/resources/rotavirus_rf.ann.vcf.gz
```

after:

```
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|frameshift_variant&start_lost|HIGH|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|p.Met1fs||1/297|1/98||;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|disruptive_inframe_deletion|MODERATE|Gene_20_616|Gene_20_616|transcript|AAG15311.1|protein_coding|1/1|c.57_59delAAA|p.Lys19del|57/597|57/597|19/198||;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|upstream_gene_variant|MODIFIER|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|||||2|;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
```


