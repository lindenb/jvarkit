**JfxNgs** is a java-based user interface used to display data from Variant Call Format ( VCF) files, or BAM files.

#Screenshots

https://twitter.com/yokofakun/status/831524118541561857

https://twitter.com/yokofakun/status/830480940052930560

https://twitter.com/yokofakun/status/830022958643036162

https://twitter.com/yokofakun/status/829360900150816769

https://twitter.com/yokofakun/status/829270290626727936

https://twitter.com/yokofakun/status/829264066371780612

https://twitter.com/yokofakun/status/828963166860214273

https://twitter.com/yokofakun/status/828626410814238720

# Supported Files

* local VCF files must be indexed with tabix or tribble
* remote VCF files must be indexed with tabix (there must be an associated '.tbi' files)
* local and remote BAM files must be indexed (there must be an associated '.bai' files)
* VCF and BAM file must have a valid header and a sequence dictionary
* servers hosting the remote VCFs and BAMs must support `Accept-Range:` ( https://en.wikipedia.org/wiki/Byte_serving )

# WebStart version

A java webstart version is currenlty available at: http://redonlab.univ-nantes.fr/public_html/jnlp/jfxngs/ , or you can launch it from the command-line:

```
javaws "http://redonlab.univ-nantes.fr/public_html/jnlp/jfxngs/JfxNgs.jnlp"
```

Security: For now I don't know how to configure my application to be always accepted by webstart.

You must add `http://redonlab.univ-nantes.fr` to the 'Exception Site List?' of the java Control-Panel . See https://www.java.com/en/download/faq/exception_sitelist.xml
The program will store temporary files in your TMPDIR (for example, remote index files are downloaded on your computer), preferences and need to read the local VCF/BAM files.


## The main screen

* The main screen is the first opened window.
* In the 'File' menu you can open a BAM file, a VCF file, or specify a remote URL to open a remote BAM or a remote VCF file.

Example remote BAM url: OPTIONS XHR https://data.broadinstitute.org/igvdata/BodyMap/hg19/IlluminaHiSeq2000_BodySites/brain_merged/accepted_hits.bam

* There is a text entry to specify a genomic location. This text field can be used to change the genomic locations of all the windows.
* closing the main screen will exit the application.

## The VCF window

* Local VCF files must be indexed with tribble or tabix
* a **pedigree** with the suffix "*.ped" can be associated to a VCF if they have the same basename. e.g:

```
input.vcf.gz
input.vcf.gz.tbi
input.ped
```


* Remote VCF files must be indexed with tabix.
* The File menu contains items to load new VCF/bam files as well as items to save the tables. An item "Save variants as" will filter the variants from the source file and save them in another vcf.
* The spinner 'Limit' is used to specify the maximun number of variants **before** filtration that will be read from the VCF input.
* The 'Goto' field is used to select variant of a genomic location. Syntax is 'chrom' or 'chrom:start' or 'chrom:start-end' .
* The 'IGV' button will move IGV in the region of the selected variant.
* The INFO tab is a table containing the definitions of the INFO fields found in the VCF header
* The FORMAT tab is a table containing the definitions of the FORMAT fields found in the VCF header
* The FILTER tab is a table containing the definitions of the FILTER fields found in the VCF header
* The 'Dict' tab is a table containing the reference dictionary, the 'contig' lines found in the VCF header.
* The 'JS' tab, is the javascript area. It's a powerful way to filter the variants. The text-area can be used to specify a javascript script that will be used to filter the variants. The script will inject 'variant' , the current variant as a  **VariantContext** ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html ) , 'header' a  VCFheader ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html ). It also inject a variable 'tools' which, at the time of writing, have the following methods:

```java
class VcfTool {
	public List<AnnPredictionParser.AnnPrediction> getAnnPredictions(final VariantContext ctx);
	public List<VepPredictionParser.VepPrediction> getVepPredictions(final VariantContext ctx);
	public boolean isMendelianIncompatibility(final Genotype child,final Genotype parent);
	public boolean isMendelianIncompatibility(final Genotype child,final Genotype father,final Genotype mother);
	}
```

The script should return a boolean to accept (TRUE) or reject (FALSE) the variant. For example, the following script will Remove all sites with more than 5% of no-call.

```javascript
var n=0.0;
for(var i=0;i< variant.getNSamples();i++ )
	{
	if(variant.getGenotype(i).isNoCall()) n++;
	}
variant.getNSamples()>0 && n/variant.getNSamples() < 0.05
```

There is a text field to set a FILTER name. If this text area is not empty, the variant won't be discarded but the column FILTER will be set for this variant.
The 'Save...' will save the current script on disk
The 'Open...' will load a script.
The 'Validate' button will validate the syntax of the javascript program.

* The variant table display the VCF rows.
* The filter table contains the FILTERs for the currently selected variant.
* The genotype table contains the genotypes for the currently selected variant. Two checkboxes can be used to hide "Hom-Ref" and "No-Call" variants.
* A pie-chart (yes, I know you hate pie-charts) is used to display the types of genotypes.
* The INFO table displays the data found in the INFO column. A text-field can be used to filter the info using a regular expression.
* The 'ANN' tab displays the annotations (called by SnpEff for example).
* The 'Stats' menu will generate various statistics about the VCF file. It can be applied to the current variant or to the whole VCF file.
* The 'Snippets' menu contains some javascript examples.

## The BAM window

* Local and BAM files must be indexed (there must be an associated *.bai file).
* The File menu contains items to load new VCF/bam files, as well as items to save the tables. An item "Save variants as" will filter the reads from the BAM file and save them in another new BAM.
* The spinner 'Limit' is used to specify the maximun number of short-reads **before** filtration that will be read from the BAM input.
* The 'Flags' button can be used to select reads on the SAM flag.
* The 'Goto' field is used to select variant of a genomic location. Syntax is 'unmapped' or 'chrom' or 'chrom:start' or 'chrom:start-end' .
* The 'IGV' button will move IGV in the region of the selected variant.
* The 'Header' contains the BAM header.
* The 'ReadGroup' tab is a table containing the 'Read-Groups' found in the BAM header.
* The 'PG' tab is a table containing the 'Program-Group' found in the BAM header.
* The 'Dict' tab is a table containing the reference dictionary, the 'contig' lines found in the BAM header.
* The 'JS' tab, is the javascript area. It's a powerful way to filter the short-reads.
The text-area can be used to specify a javascript script that will be used to filter the short-reads.
The script will inject 'record' , the current short-read as a  **SAMRecord** ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html ) , 'header' a  SAMFileHeader ( http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html ). It also inject a variable 'tools' which, at the time of writing, have the following methods:

```java
class BamTool {
	public String reverseComplement(final String sequenceData);
	}
```

The script should return a boolean to accept (TRUE) or reject (FALSE) the variant. For example , the following script will remove all the soft/hard clipped variants.

```javascript
record.getReadUnmappedFlag()==false && (
	record.getCigarString().contains("S") ||
	record.getCigarString().contains("H")
	);
```

The 'Save...' will save the current script on disk
The 'Open...' will load a script.
The 'Validate' button will validate the syntax of the javascript program.

* The 'Canvas' tab displays a mini-genome browser.
* The 'Reads' table display the BAM rows.
* The Flags table contains the sam-fags for the currently selected read.
* The key-value table contains meta-data of the selected read.
* A cigar table displays the positions of all the bases of the selected read.
* The 'Stats' menu will generate various statistics about the BAM file. It can be applied to the current short-read or to the whole BAM file.
* The 'Snippets' menu contains some javascript examples.

## Bialcidae

Both  BAM and VCF windows have an interface to Bioalcidae: https://github.com/lindenb/jvarkit/wiki/BioAlcidae .
 

## Compiling the standalone version

Compilation, requires Oracle java SDK 8, curl, GNU make. After cloning the github repo, `cd`in the jvarkit folder and type

```
make jfxngs
```

This should create a executable jar file in `dist/jfxngs.jar`

## Running the standalone version

```
java -jar dist/jfxngs.jar

```

you can also specify the URLs/File on the command line


```
java -jar dist/jfxngs.jar "http://path/to/file.vcf.gz" "/path/to/file.bam"

```

## Author

Pierre Lindenbaum PhD Insitut du Thorax , UMR1087 , Nantes, France. @yokofakun
