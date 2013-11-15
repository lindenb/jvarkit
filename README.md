JVARKIT
=======

Java utilities for Next Generation Sequencing

Pierre Lindenbaum PhD

http://plindenbaum.blogspot.com

[@yokofakun](https://twitter.com/yokofakun)

## Download and install

see [Download and Install](https://github.com/lindenb/jvarkit/wiki/Compilation)

##Tools

###Main
<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SplitBam">SplitBam</a></th><td>Split a BAM by chromosome group. Creates EMPTY bams if no reads was found for a given group. </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamJS">SamJS</a></th><td>Filtering a SAM/BAM with javascript (rhino).</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFilterJS">VCFFilterJS</a></th><td>Filtering a VCF with javascript (rhino)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SortVCFOnRef">SortVCFOnRef<a></th><td>Sort a VCF using the order of the chromosomes in a REFerence index.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Illuminadir">Illuminadir</a></th><td>Create a structured (**JSON** or **XML**) representation of a directory containing some Illumina FASTQs.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats04">BamStats04<a></th><td>Coverage statistics for a BED file. It uses the Cigar string instead of the start/end to compute the coverage</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats01">BamStats01<a></th><td>Statistics about the reads in a BAM.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFBed">VCFBed<a></th><td>Annotate a VCF with the content of a BED file indexed with tabix.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFPolyX">VCFPolyX<a></th><td>Number of repeated REF bases around POS.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFBigWig">VCFBigWig<a></th><td>Annotate a VCF with the data of a bigwig file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFTabixml">VCFTabixml<a></th><td>Annotate a value from a vcf+xml file.4th column of the BED indexed with TABIX is a XML string.</td></tr>
</table>

###Less used, but useful
<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFPredictions">VCFPredictions<a></th><td>Basic variant prediction using UCSC knownGenes.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FindCorruptedFiles">FindCorruptedFiles<a></th><td>Reads filename from stdin and prints corrupted NGS files (VCF/BAM/FASTQ). </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCF2XML">VCF2XML</a></th><td>Transforms a VCF to XML. </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFAnnoBam">VCFAnnoBam<a></th><td>Annotate a VCF with the Coverage statistics of a BAM file + BED file of capture. It uses the Cigar string instead of the start/end to get the voverage</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFTrio">VCFTrio<a></th><td>Check for mendelian incompatibilities in a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamGrep">SamGrep<a></th><td>Search reads in a BAM</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFixIndels">VCFFixIndels<a></th><td>Fix samtools INDELS for @SolenaLS</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/NgsFilesSummary">NgsFilesSummary<a></th><td>Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..).</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/NoZeroVariationVCF">NoZeroVariationVCF<a></th><td>creates a VCF containing one fake variation if the input is empty.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/HowManyBamDict">HowManyBamDict<a></th><td>for @abinouze : quickly find the number of distinct BAM Dictionaries from a set of BAM files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/ExtendBed">ExtendBed<a></th><td>Extends a BED file by 'X' bases.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/CmpBams">CmpBams<a></th><td>Compare two or more BAMs.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/IlluminaFastqStats">IlluminaFastqStats<a></th><td>Statistics on Illumina Fastqs</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Bam2Raster">Bam2Raster<a></th><td>Save a BAM alignment as a PNG image.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfRebase">VcfRebase<a></th><td>Finds restriction sites overlapping variants in a VCF file</td></tr>
</table>


###One shots, Answers to biostar,...
<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BlastMapAnnots">BlastMapAnnots<a></th><td>Maps uniprot/genbank annotations on a blast result. See http://www.biostars.org/p/76056</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfViewGui">VcfViewGui<a></th><td>Simple java-Swing-based VCF viewer.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar81455">Biostar81455<a></th><td>Defining precisely the genomic context based on a position http://www.biostars.org/p/81455/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/MapUniProtFeatures">MapUniProtFeatures<a></th><td>map Uniprot features on reference genome.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar86363">Biostar86363<a></th><td>Set genotype of specific sample/genotype comb to unknown in multisample vcf file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FixVCF">FixVCF<a></th><td>Fix a VCF HEADER when I forgot to declare a FILTER or an INFO field in the HEADER</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar78400">Biostar78400<a></th><td>Add the read group info to the sam file on a per lane basis</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar78285">Biostar78285<a></th><td>Extract regions of genome that have 0 coverage See http://www.biostars.org/p/78285/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar77288">Biostar77288<a></th><td>Low resolution sequence alignment visualization http://www.biostars.org/p/77288/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar77828">Biostar77828<a></th><td>Divide the human genome among X cores, taking into account gaps See http://www.biostars.org/p/77828/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar76892">Biostar76892<a></th><td>Fix strand of two paired reads close but on the same strand http://www.biostars.org/p/76892/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFCompareGT">VCFCompareGT<a></th><td>VCF : compare genotypes of two or more callers for the same samples.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SAM4WebLogo">SAM4WebLogo<a></th><td>Creates an Input file for BAM + WebLogo.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SAM2Tsv">SAM2Tsv<a></th><td>Tabular view of each base of the reads vs the reference.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar84786">Biostar84786<a></th><td>Table transposition</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCF2SQL">VCF2SQL<a></th><td>Generate the SQL code to insert a VCF into a database</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Bam4DeseqIntervals">Bam4DeseqIntervals<a></th><td>creates a table for DESEQ with the number of reads within a sliding window for multiple BAMS</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFStripAnnotations">VCFStripAnnotations<a></th><td>Removes one or more field from the INFO column from a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFGeneOntology">VCFGeneOntology<a></th><td>Finds the GO terms for VCF annotated with SNPEFF or VEP</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFilterGO">VCFFilterGO<a></th><td>Set the VCF FILTERs on VCF files annotated with SNPEFF or VCP testing wether a Gene belong or not to the descendants of a GO term.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar86480">Biostar86480<a></th><td>Genomic restriction finder See http://www.biostars.org/p/86480/</td></tr>
</table>
