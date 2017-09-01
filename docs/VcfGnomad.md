# VcfGnomad

Peek annotations from gnomad


## Usage

```
Usage: vcfgnomad [options] Files
  Options:
    -ac, --alleleconcordance
      ALL Alt allele must be found in gnomad before setting a FILTER
      Default: false
    --bufferSize
      When we're looking for variant in Exac, load the variants for 'N' bases 
      instead of doing a random access for each variant
      Default: 100000
    -filtered, --filtered
      Skip Filtered User Variants
      Default: false
    -filteredGnomad, --filteredGnomad
      [20170706] Skip Filtered GNOMAD Variants
      Default: false
    -gf, --gnomadFilter
      if defined, add this FILTER when the variant is found in nomad
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -m, --manifest
      manifest file descibing how to map a contig to an URI . 3 columns: 1) 
      exome|genome 2) contig 3) path or URL.
    --noAlleleCount
      do Not Insert AC /Allele Count
      Default: false
    --noAlleleFreq
      do Not Insert AF /Allele Freq.
      Default: false
    --noAlleleNumber
      do Not Insert AN /Allele Number
      Default: false
    -noMultiAltGnomad, --noMultiAltGnomad
      [20170706] Skip Multi Allelic GNOMAD Variants
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --streaming
      [20170707] Don't use tabix random-access (which are ok for small inputs) 
      but you a streaming process (better to annotate a large WGS file). 
      Assume dictionaries are sorted the same way.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * gnomad


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfgnomad
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomad.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomad.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomad** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Manifest
 
 the manifest is a tab delimited file containing 3 columns. It's used to map a contig to a URI
 
   * 1st column is a keyword 'exome' or 'genome'
   * 2d column is a contig name e.g: '1' .  Use '*' for 'any' chromosome
   * 3d column is a URL or file path where to find the data
 
 
## Example:
 
 ```
  curl -s "https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz" |\
     gunzip -c | head -n 400 |\
     java  -jar ~/src/jvarkit-git/dist/vcfgnomad.jar -ac -gf IN_GNOMAD 

 (...)
 1	13595	.	AGT	A	379.68	AC0;IN_GNOMAD;RF	AB_HIST_ALL=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_MEDIAN=1.44068e-01;AC=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_Female=0;AC_Male=0;AC_NFE=0;AC_OTH=0;AC_POPMAX=.;AC_SAS=0;AC_raw=1;AF=0.00000e+00;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_Female=0.00000e+00;AF_Male=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;AF_POPMAX=.;AF_SAS=0.00000e+00;AF_raw=9.99900e-06;AN=50778;AN_AFR=4986;AN_AMR=10892;AN_ASJ=1274;AN_EAS=7560;AN_FIN=694;AN_Female=24940;AN_Male=25838;AN_NFE=17556;AN_OTH=1486;AN_POPMAX=.;AN_SAS=6330;AN_raw=100010;AS_FilterStatus=RF|AC0;AS_RF=1.49748e-01;BaseQRankSum=-4.60000e-01;CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034|YES|||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.561_562delTG||558-559||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.847_848delTG||844-845||||||1||1||deletion|1|HGNC|37102|YES||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene|||||||||||1|807|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.840_841delTG||837-838||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.678_679delTG||675-676||||||1||1||deletion|1|HGNC|37102|||||||||||||3|||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene|||||||||||1|814|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene|||||||||||1|766|-1||deletion|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site|||||||||||1||||deletion|1||||||||||||||||||||||||||||||||||||||||||||;ClippingRankSum=5.63000e-01;DP=2519792;DP_HIST_ALL=20921|3680|466|85|62|97|652|4365|4551|3656|2891|2039|1464|1114|954|811|688|497|352|310;DP_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;DP_MEDIAN=118;DREF_MEDIAN=3.98107e-38;FS=1.59250e+01;GC=25389,0,0;GC_AFR=2493,0,0;GC_AMR=5446,0,0;GC_ASJ=637,0,0;GC_EAS=3780,0,0;GC_FIN=347,0,0;GC_Female=12470,0,0;GC_Male=12919,0,0;GC_NFE=8778,0,0;GC_OTH=743,0,0;GC_SAS=3165,0,0;GC_raw=50004,1,0;GQ_HIST_ALL=11211|8535|2038|2055|803|203|195|95|28|49|65|37|115|64|88|117|164|34|237|23872;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;GQ_MEDIAN=99;Hom=0;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_Female=0;Hom_Male=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;Hom_raw=0;InbreedingCoeff=-4.37000e-02;MQ=3.15600e+01;MQRankSum=-8.97000e-01;POPMAX=.;QD=3.22000e+00;ReadPosRankSum=-1.23200e+00;SOR=1.09000e-01;VQSLOD=-1.83100e+00;VQSR_NEGATIVE_TRAIN_SITE;VQSR_culprit=QD;gnomad.exome.AC_AFR=0;gnomad.exome.AC_AMR=0;gnomad.exome.AC_ASJ=0;gnomad.exome.AC_EAS=0;gnomad.exome.AC_FIN=0;gnomad.exome.AC_Female=0;gnomad.exome.AC_Male=0;gnomad.exome.AC_NFE=0;gnomad.exome.AC_OTH=0;gnomad.exome.AC_raw=1;gnomad.exome.AN_AFR=4986;gnomad.exome.AN_AMR=10892;gnomad.exome.AN_ASJ=1274;gnomad.exome.AN_EAS=7560;gnomad.exome.AN_FIN=694;gnomad.exome.AN_Female=24940;gnomad.exome.AN_Male=25838;gnomad.exome.AN_NFE=17556;gnomad.exome.AN_OTH=1486;gnomad.exome.AN_raw=100010;gnomad.genome.AC_AFR=0;gnomad.genome.AC_AMR=0;gnomad.genome.AC_ASJ=0;gnomad.genome.AC_EAS=0;gnomad.genome.AC_FIN=0;gnomad.genome.AC_Female=0;gnomad.genome.AC_Male=0;gnomad.genome.AC_NFE=0;gnomad.genome.AC_OTH=0;gnomad.genome.AC_raw=1;gnomad.genome.AN_AFR=8680;gnomad.genome.AN_AMR=794;gnomad.genome.AN_ASJ=224;gnomad.genome.AN_EAS=1592;gnomad.genome.AN_FIN=3490;gnomad.genome.AN_Female=13274;gnomad.genome.AN_Male=16168;gnomad.genome.AN_NFE=13754;gnomad.genome.AN_OTH=908;gnomad.genome.AN_raw=30500

 
 ```

## Note to self: Another alternative with VariantAnnotator,

but I think it slower...

(javascript / Makefile generation)

```javascript
out.print(" ${java.exe} -jar ${gatk.jar} -R $(REF) -L $(addsuffix .tmp.vcf,$@) -T VariantAnnotator --variant $(addsuffix .tmp.vcf,$@) -o $(addsuffix .tmp2.vcf,$@) --resourceAlleleConcordance ");

out.print(" --resource:gnomad_exome /commun/data/pubdb/broadinstitute.org/gnomad/release-170228/vcf/exome/gnomad.exomes.r2.0.1.sites.vcf.gz ");
out.print("$(foreach A,${GFIELDS}, -E gnomad_exome.${A} ) ");

var genome="/commun/data/pubdb/broadinstitute.org/gnomad/release-170228/vcf/genome/gnomad.genomes.r2.0.1.sites."+chrom+".vcf.gz";

out.print("$(if $(realpath "+genome+"), --resource:gnomad_genome  "+genome+"  $(foreach A,${GFIELDS}, -E gnomad_genome.${A} ) )");
```



## Generating jar helper for knime
  
(for the people in my lab)
  
generate big jar

```
$ cd jvarkit
$ rm -rf tmp && mkdir tmp && echo '1.jar:2.jar:...N.jar:vcfgnomad.jar' | tr ":" "\n" | sort | uniq | while read F; do unzip -o $F -d tmp ; done && jar cvf vcfgnomad4knime.jar -C tmp . && rm -rf tmp
```
  
Open KNIME

we're going to create the following workflow : http://imgur.com/a/QcrKW

* create a new Node `java Snippet`
* in the tab 'additional libraries', add 'vcfgnomad4knime.jar'.
* in the tab 'java snippet'. Declare the following inputs: `c_CHROM,c_POS,c_REF,c_ALT`, the output string `GNOMAD`.

And insert the following code:

```java
// Your custom imports:
import  com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad.KnimeAdapter;
// Enter your code here:

System.setProperty("http.proxyHost","cache.ha.univ-nantes.fr");
System.setProperty("https.proxyHost","cache.ha.univ-nantes.fr");
System.setProperty("http.proxyPort","3128");
System.setProperty("https.proxyPort","3128");

final KnimeAdapter app= new KnimeAdapter();
if(app.instanceMain(new String[]{"-ac",c_CHROM,String.valueOf(c_POS),c_REF,c_ALT})==0)
	{
	out_GNOMAD = app.getOutputString();
	}
else
	{
	out_GNOMAD = ".";
	}

// expression end

```
 


