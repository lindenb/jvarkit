JVARKIT
=======

Java utilities for Bioinformatics

[![Build Status](https://travis-ci.org/lindenb/jvarkit.svg)](https://travis-ci.org/lindenb/jvarkit)

## Warning

<div style="background-color:red;">Since 2015-12-10, I'm slowly moving to a <b>XML</b>-based description of my tools and I'm now using <b>java8</b> .
I do my best to change all those tools. The documentation in the wiki migh be out of date. If you think you have found an error leave a message at https://github.com/lindenb/jvarkit/issues .
Furthermore, I'm now using the "Apache Commons CLI library" for parsing the command line. For the arguments takings as input more than one parameter, you might have to
add a double dash '--' to separate with the input files. 
</div>


## Author

Pierre Lindenbaum PhD

http://plindenbaum.blogspot.com

[@yokofakun](https://twitter.com/yokofakun)

## Cite

see [Cite](https://github.com/lindenb/jvarkit/wiki/Citing)

## Download and install

See [Download and Install](https://github.com/lindenb/jvarkit/wiki/Compilation)

##Tools


<table>
<tr><th>Tool</th><th>Description</th></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SplitBam">SplitBam</a></th><td>Split a BAM by chromosome group. Creates EMPTY bams if no reads was found for a given group. </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamJS">SamJS</a></th><td>Filtering a SAM/BAM with javascript (rhino).</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFFilterJS">VCFFilterJS</a></th><td>Filtering a VCF with javascript (rhino)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SortVCFOnRef">SortVCFOnRef<a></th><td>Sort a VCF using the order of the chromosomes in a REFerence index.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Illuminadir">Illuminadir</a></th><td>Create a structured (**JSON** or **XML**) representation of a directory containing some Illumina FASTQs.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats04">BamStats04<a></th><td>Coverage statistics for a BED file. It uses the Cigar string instead of the start/end to compute the coverage</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats05">BamStats05<a></th><td>same as BamStats04 but group by gene</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats01">BamStats01<a></th><td>Statistics about the reads in a BAM.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFBed">VCFBed<a></th><td>Annotate a VCF with the content of a BED file indexed with tabix.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFPolyX">VCFPolyX<a></th><td>Number of repeated REF bases around POS.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFBigWig">VCFBigWig<a></th><td>Annotate a VCF with the data of a bigwig file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFTabixml">VCFTabixml<a></th><td>Annotate a value from a vcf+xml file.4th column of the BED indexed with TABIX is a XML string.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/GroupByGene">GroupByGene<a></th><td>Group VCF data by gene/transcript.</td></tr>
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
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqRevComp">FastqRevComp<a></th><td>Reverse complement a FASTQ file for mate-pair alignment</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/PicardMetricsToXML">PicardMetricsToXML<a></th><td>Convert picards metrics file to XML.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Bam2Wig">Bam2Wig</thd><td>Bam to Wiggle converter</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/TViewWeb">TViewWeb</thd><td>CGI/Web based version of samtools tview</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfRegistryWeb">VcfRegistryWeb</thd><td>CGI/Web tool printing all variants at a given position for a collection VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BlastMapAnnots">BlastMapAnnots<a></th><td>Maps uniprot/genbank annotations on a blast result. See http://www.biostars.org/p/76056</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfViewGui">VcfViewGui<a></th><td>Simple java-Swing-based VCF viewer.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamViewGui">BamViewGui<a></th><td>Simple java-Swing-based BAM viewer.</td></tr>
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
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFStripAnnotations">VCFStripAnnotations<a></th><td>Removes one or more field from the INFO column from a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFGeneOntology">VCFGeneOntology<a></th><td>Finds and filters the GO terms for VCF annotated with SNPEFF or VEP</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar86480">Biostar86480<a></th><td>Genomic restriction finder See http://www.biostars.org/p/86480/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamToFastq">BamToFastq<a></th><td>Shrink your FASTQ.bz2 files by 40+% using this one weird tip  by ordering them by alignment to reference</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/PadEmptyFastq">PadEmptyFastq<a></th><td>Pad empty fastq sequence/qual with N/#</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamFixCigar">SamFixCigar<a></th><td>Replace 'M'(match) in SAM cigar by 'X' or '='</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FixVcfFormat">FixVcfFormat<a></th><td>Fix PL format in VCF. Problem is described in http://gatkforums.broadinstitute.org/discussion/3453</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfToRdf">VcfToRdf<a></th><td>Convert a VCF to RDF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFShuffle">VcfShuffle<a></th><td>Shuffle a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/DownSampleVcf">DownSampleVcf<a></th><td>Down sample a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfHead">VcfHead<a></th><td>Print the first variants of a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfTail">VcfTail<a></th><td>Print the last variants of a VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfCutSamples">VcfCutSamples<a></th><td>Select/Exclude some samples from a VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfStats">VcfStats<a></th><td>Generate some statistics from a VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfSampleRename">VcfSampleRename<a></th><td>Rename Samples in a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcffilterSequenceOntology">VcffilterSequenceOntology<a></th><td>Filter a VCF on Seqence Ontology (SO).</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar59647">Biostar59647<a></th><td>position of mismatches per read from a sam/bam file (XML) See http://www.biostars.org/p/59647/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfRenameChromosomes">VcfRenameChromosomes<a></th><td>Rename chromosomes in a VCF (eg. convert hg19/ucsc to grch37/ensembl)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamRenameChromosomes">BamRenameChromosomes<a></th><td>Rename chromosomes in a BAM (eg. convert hg19/ucsc to grch37/ensembl)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BedRenameChromosomes">BedRenameChromosomes<a></th><td>Rename chromosomes in a BED (eg. convert hg19/ucsc to grch37/ensembl)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BlastnToSnp">BlastnToSnp<a></th><td>Map variations from a BLASTN-XML file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Blast2Sam">Blast2Sam<a></th><td>Convert a BLASTN-XML input to SAM</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfMapUniprot">VcfMapUniprot<a></th><td>Map uniprot features on VCF annotated with VEP or SNPEff.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfCompare">VcfCompare<a></th><td>Compare two VCF files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfBiomart">VcfBiomart<a></th><td>Annotate a VCF with the data from Biomart.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfLiftOver">VcfLiftOver<a></th><td>LiftOver a VCF file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BedLiftOver">BedLiftOver<a></th><td>LiftOver a BED file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfConcat">VcfConcat<a></th><td>Concatenate VCF files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/MergeSplittedBlast">MergeSplittedBlast<a></th><td>Merge Blast hit from a splitted database</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FindMyVirus">FindMyVirus<a></th><td>Virus+host cell : split BAM into categories.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar90204">Biostar90204<a></th><td>linux split equivalent for BAM file .</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfJaspar">VcfJaspar<a></th><td>Finds JASPAR profiles in VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/GenomicJaspar">GenomicJaspar<a></th><td>Finds JASPAR profiles in Fasta</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfTreePack">VcfTreePack<a></th><td>Create a TreeMap from one or more VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamTreePack">BamTreePack<a></th><td>Create a TreeMap from one or more Bam.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqRecordTreePack">FastqRecordTreePack<a></th><td>Create a TreeMap from one or more Fastq files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/WorldMapGenome">WorldMapGenome<a></th><td>Map bed file to Genome + geographic data.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/AddLinearIndexToBed">AddLinearIndexToBed<a></th><td>Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFComm">VCFComm<a></th><td>Compare mulitple VCF files, ouput a new VCF file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfIn">VcfIn<a></th><td>Prints variants that are contained/not contained into another VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar92368">Biostar92368<a></th><td>Binary interactions depth See also http://www.biostars.org/p/92368</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFStopCodon">VCFStopCodon<a></th><td>TODO</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqGrep">FastqGrep<a></th><td>Finds reads in fastq files</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfCadd">VcfCadd<a></th><td>Annotate a VCF with Combined Annotation Dependent Depletion (CADD) data.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SortVCFOnInfo">SortVCFOnInfo<a></th><td>sort a VCF using a field in the INFO column</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamChangeReference">SamChangeReference<a></th><td>TODO</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamExtractClip">SamExtractClip<a></th><td>TODO</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/GCAndDepth">GCAndDepth<a></th><td>Extracts GC% and depth for multiple bam using a sliding window.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/MsaToVcf">MsaToVcf<a></th><td>Getting a VCF file from a CLUSTAW or FASTA alignment</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/CompareBamAndBuild">CompareBamAndBuild<a></th><td>Compare two  BAM files mapped on two different builds. Requires a liftover chain file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/KnownGenesToBed">KnownGenesToBed<a></th><td>Convert UCSC KnownGene to BED.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar95652">Biostar95652<a></th><td>Drawing a schematic genomic context tree. See also http://www.biostars.org/p/95652/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamToPsl">SamToPsl<a></th><td>Convert SAM/BAM to PSL  or BED12 .</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BWAMemNOp">BWAMemNOp<a></th><td>merge the SA:Z:* attributes of a read mapped with bwa-mem and prints a read containing a cigar string with 'N' (Skipped region from the REF).</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqEntropy">FastqEntropy<a></th><td>Compute the Entropy of a Fastq file (distribution of the length(gzipped(sequence)))</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/NgsFilesScanner">NgsFilesScanner<a></th><td>Build a persistent database of NGS file. Dump as XML.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SigFrame">SigFrame<a></th><td>GUI displaying CGH data</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar103303">Biostar103303<a></th><td>Calculate Percent Spliced In (PSI)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFComparePredictions">VCFComparePredictions<a></th><td>Compare the variant predictions of VCFs</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BackLocate">BackLocate<a></th><td>Map a  position in a protein back to the genomic coordinates.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FindAVariation">FindAVariation<a></th><td>Search for variations in a set of VCF files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/AlleleFrequencyCalculator">AlleleFrequencyCalculator<a></th><td>VCF: Alelle Frequency Calculator</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BuildWikipediaOntology">BuildWikipediaOntology<a></th><td>Build a simple RDFS/XML ontology from Wikipedia Categories.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/AlmostSortedVcf">AlmostSortedVcf<a></th><td>Sort an 'almost' sorted VCF using an in-memory buffer.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar105754">Biostar105754<a></th><td>bigwig: peak distance from specific genomic BED region</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfRegulomeDB">VcfRegulomeDB<a></th><td>Annotate a VCF with the RegulomeDB data (http://regulome.stanford.edu/)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar106668">Biostar106668<a></th><td>unmark duplicates (deprecated)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BatchIGVPictures">BatchIGVPictures<a></th><td>GUI: Batch pictures with IGV</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/PubmedDump">PubmedDump<a></th><td>Dump pubmed data as XML.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamIndexReadNames">BamIndexReadNames<a></th><td>Build a dictionary of read names to be searched with BamQueryReadNames.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamQueryReadNames">BamQueryReadNames<a></th><td>Query a Bam file indexed with BamIndexReadNames.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqShuffle">FastqShuffle<a></th><td>Shuffle Fastq files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqSplitInterleaved">FastqSplitInterleaved<a></th><td>Split interleaved Fastq files</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/PubmedFilterJS">PubmedFilterJS<a></th><td>Filters pubmed XML using javascript.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/ReferenceToVCF">ReferenceToVCF<a></th><td>Creates a VCF containing all possible substitutions in a Reference Genome..</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfEnsemblReg">VcfEnsemblReg<a></th><td>Annotate a VCF with the UCSC genome hub tracks for Ensembl Regulation.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FastqJS">FastqJS<a></th><td>Filters a FASTQ file using javascript.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Bam2SVG">Bam2SVG<a></th><td>Convert a BAM to SVG</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/LiftOverToSVG">LiftOverToSVG</a></th><td>Convert UCSC LiftOver chain files to animated SVG</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFMerge">VCFMerge</a></th><td>Combines VCF files.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FixVcfMissingGenotypes">FixVcfMissingGenotypes</a></th><td>Use BAM to fill missing genotypes in merged VCFs</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/NcbiTaxonomyToXml">NcbiTaxonomyToXml</a></th><td> Dump NCBI taxonomy tree as a hierarchical XML document</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamCmpCoverage">BamCmpCoverage</a></th><td> Creates the figure of a comparative view of the depths sample vs sample</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/FindAllCoveragesAtPosition">FindAllCoveragesAtPosition</a></th><td>Find depth at specific position in a list of BAM files</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfMultiToOne">VcfMultiToOne</a></th><td>Convert VCF with multiple samples to a VCF with one SAMPLE</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Evs2Xml">Evs2Xml</a></th><td>Download data from Exome Variant Server as XML.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfRemoveGenotypeIfInVcf">VcfRemoveGenotypeIfInVcf</a></th><td>Reset Genotypes in VCF if they've been found in another VCF indexed with tabix</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar130456">Biostar130456</a></th><td>Generate one VCF file for each sample from a multi-samples VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/UniprotFilterJS">UniprotFilterJS</a></th><td>Filter Uniprot XML with a javascript expression.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SkipXmlElements">SkipXmlElements</a></th><td>Filter XML elements with a javascript expression.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/MiniCaller">MiniCaller</a></th><td>Simple and Stupid Variant Caller designed for @AdrienLeger2</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfCompareCallersOneSample">VcfCompareCallersOneSample</a></th><td>For my colleague Julien. Compare VCF allers with VCF with one sample.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamRetrieveSeqAndQual">SamRetrieveSeqAndQual</a></th><td> Is there a tool to add seq and qual to BAM? for @sjackman </td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfEnsemblVepRest">VcfEnsemblVepRest</a></th><td>Annotate a VCF with Ensembl REST API.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfCompareCallers">VcfCompareCallers</a></th><td>Compare two VCFs and print common/exclusive information for each sample/genotype</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamStats02">BamStats02</a></th><td>Generate and explore statistics about the reads in a BAM (Sample/File/Flags/chroms/MAPQ)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BamTile">BamTile</a></th><td>Bam tiling Path.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/XContaminations">XContaminations</a></th><td>for @AdrienLeger2 : test for cross contamination between samples in same flowcell/runlane.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFJoinVcfJS">VCFJoinVcfJS</a></th><td>Join two VCF files using javascript.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar139647">Biostar139647</a></th><td>Convert Clustal/Fasta alignment to SAM/BAM</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BioAlcidae">BioAlcidae</a></th><td>Reformat bioinformatics files using javascript/rhino (~ awk)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFBedSetFilter">VCFBedSetFilter</a></th><td>Set FILTER for VCF having intersection with BED</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VCFReplaceTag">VCFReplaceTag</a></th><td>Replace the key in INFO/FORMAT/FILTER</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfIndexTabix">VcfIndexTabix</a></th><td>sort, Compress (bgz) a VCF and create tabix index on the fly.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfPeekVcf">VcfPeekVcf</a></th><td>Peek INFO Tag and ID from another VCF</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfGetVariantByIndex">VcfGetVariantByIndex</a></th><td>Access a (plain or tabix-indexed) VCF file by the i-th index.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele">VcfMultiToOneAllele</a></th><td>VCF: "one variant with N ALT alleles" to "N variants with one ALT"</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneInfo">VcfMultiToOneInfo</a></th><td>VCF: "one variant having INFO with N values" to "N variants having INFO with one value"</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BedIndexTabix">BedIndexTabix</a></th><td>Index and sort a BED on the fly with Tabix</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfToHilbert">VcfToHilbert</a></th><td>Plot a Hilbert Curve from a VCF file.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar145820">Biostar145820</a></th><td>Shuffl Bam/Subsample BAM to fixed number of alignments</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/PcrClipReads">PcrClipReads</a></th><td>Soft clip BAM files based on PCR target regions https://www.biostars.org/p/147136/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/ExtendReferenceWithReads">ExtendReferenceWithReads</a></th><td>Extending ends of REF sequence with the help of reads in BAM https://www.biostars.org/p/148089/</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/PcrSliceReads">PcrSliceReads</a></th><td>Mark PCR reads to their PCR amplicon https://www.biostars.org/p/149687/"</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SamJmx">SamJmx</a></th><td>Monitor/interrupt/break a BAM/SAM stream with java JMX</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfJmx">VcfJmx</a></th><td>Monitor/interrupt/break a VcfJmx stream with java JMX</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Gtf2Xml">Gtf2Xml</a></th><td>convert gff to XML  in order to be processed with XSLT</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/SortSamRefName">SortSamRefName</a></th><td>Sort a SAM/BAM on REF/contig and then on read/query name</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar154220">Biostar154220</a></th><td>Cap BAM to a given coverage. see https://www.biostars.org/p/154220</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/VcfToBam">VcfToBam</a></th><td>create a BAM from a VCF.</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar165777">Biostar165777</a></th><td>Split a XML file (e.g: blast)</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/BlastFilterJS">BlastFilterJS</a></th><td>Filters a XML Blast Output with a javascript expression</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar170742">Biostar170742</a></th><td>SAM to AXT converter</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar172515">Biostar172515</a></th><td>Convert Bam Index to XML</td></tr>
<tr><th><a href="https://github.com/lindenb/jvarkit/wiki/Biostar173114">Biostar173114</a></th><td>make a bam file smaller by removing unwanted information see also https://www.biostars.org/p/173114</td></tr>
</table>

