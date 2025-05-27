JVARKIT
=======

Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.
Version     : 2e415c3ad
Compilation : 20250527180525
Github      : https://github.com/lindenb/jvarkit
Issues      : https://github.com/lindenb/jvarkit/issues

## Usage

```
  java -jar jvarkit.jar [options]
```
or
```
  java -jar jvarkit.jar <command name> (other arguments)
```

## Options

 + --help show this screen
 + --help-all show all commands, including the private ones.
 + --version print version

## Compilation Installation

Please, read [how to run and install jvarkit](JvarkitCentral.md)

## Tools

### BAM Visualization

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [bam2raster](Bam2Raster.md) | BAM to raster graphics |  |  |
| [bam2svg](BamToSVG.md) | BAM to Scalar Vector Graphics (SVG) | 20141013 | 20210728 |
| [biostar139647](Biostar139647.md) | Convert alignment in Fasta/Clustal format to SAM/BAM file |  |  |
| [biostar145820](Biostar145820.md) | subsample/shuffle BAM to fixed number of alignments. | 20150615 | 20211005 |
| [lowresbam2raster](LowResBam2Raster.md) | Low Resolution BAM to raster graphics | 20170523 | 20211126 |
| [mkminibam](MakeMiniBam.md) | Creates an archive of small bams with only a few regions. | 20190410 | 20221019 |
| [plotsashimi](PlotSashimi.md) | Print Sashimi plots from Bam | 20191117 | 20191104 |
| [prettysam](PrettySam.md) | Pretty SAM alignments | 20171215 | 20211105 |
| [sv2svg](SvToSVG.md) | BAM to SVG. Used to display structural variations. | 20181115 | 20230505 |
| [wgscoverageplotter](WGSCoveragePlotter.md) | Whole genome coverage plotter | 20201125 | 20230505 |

### CNV/SV

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [bammatrix](BamMatrix.md) | Bam matrix, inspired from 10x/loupe | 20190620 | 20211206 |
| [cnvpanelofnormal](CNVPaneOfNormal.md) | Call CNV from panel of normal computed with 'basecoverage' | 20241123 | 20241123 |
| [cnvtview](CnvTView.md) | Text visualization of bam DEPTH for multiple regions in a terminal | 20181018 | 20210412 |
| [coveragegrid](CoverageGrid.md) | Display an image of depth to display any anomaly an intervals+bams as a grid image | 20241009 | 20241021 |
| [coverageplotter](CoveragePlotter.md) | Display an image of depth to display any anomaly an intervals+bams | 20200605 | 20241009 |
| [indexcov2vcf](IndexCovToVcf.md) | convert indexcov data to vcf | 20200528 | 20400313 |
| [samfindclippedregions](SamFindClippedRegions.md) | Fins clipped position in one or more bam. | 20140228 | 20220329 |
| [swingindexcov](SwingIndexCov.md) | indexcov visualization | 2020511 | 2020512 |
| [vcfstrech2svg](VcfStrechToSvg.md) | another VCF to SVG | 20210304 | 20210309 |
| [wescnvsvg](WesCnvSvg.md) | SVG visualization of bam DEPTH for multiple regions | 20180726 | 20210726 |

### Functional prediction

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [backlocate](BackLocate.md) | Mapping a mutation on a protein back to the genome. | 20140619 | 20190820 |
| [groupbygene](GroupByGene.md) | Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff | 20140531 | 20220529 |

### BED Manipulation

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [bedcluster](BedCluster.md) | Clusters a BED file into a set of BED files. | 20200130 | 2050129 |
| [bedclustername](BedClusterName.md) | Clusters a BED file into a set of BED files using the 4th column of the bed name. | 2050428 | 2050428 |
| [bedmergecnv](BedMergeCnv.md) | Merge continuous sorted bed records if they overlap a fraction of their lengths. | 20200330 | 20200603 |
| [bednonoverlappingset](BedNonOverlappingSet.md) | Split a Bed file into non-overlapping data set. | 20180607 | 20200408 |
| [bedrenamechr](BedRenameChromosomes.md) | Convert the names of the chromosomes in a Bed file | 20190503 | 20240515 |
| [setfile2bed](SetFileToBed.md) | Convert setfile to bed | 20210125 | 20240724 |
| [setfilecluster](SetFileCluster.md) | Cluster records of setfiles into files containing a sum to basepaires close to 'x' bp | 20210125 | 20240724 |
| [setfilefrombed](SetFileFromBed.md) | Convert bed chrom/start/end/name sorted on 4th column to set file | 20210125 | 20240724 |
| [setfiletools](SetFileTools.md) | Utilities for the setfile format | 20210125 | 20220426 |

### Biostars

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [biostar103303](Biostar103303.md) | Calculate Percent Spliced In (PSI). |  |  |
| [biostar105754](Biostar105754.md) | bigwig : peak distance from specific genomic region | 20140708 | 20220110 |
| [biostar165777](Biostar165777.md) | Split a XML file | 20151114 | 20151114 |
| [biostar170742](Biostar170742.md) | convert sam format to axt Format | 20151228 | 20210412 |
| [biostar172515](Biostar172515.md) | Convert BAI to XML |  |  |
| [biostar173114](Biostar173114.md) | make a bam file smaller by removing unwanted information see also https://www.biostars.org/p/173114/ |  |  |
| [biostar175929](Biostar175929.md) | Construct a combination set of fasta sequences from a vcf | 20160208 | 20211012 |
| [biostar178713](Biostar178713.md) | split bed file into several bed files where each region is separated of any other by N bases | 20160226 | 20200818 |
| [biostar214299](Biostar214299.md) | Extract allele specific reads from bamfiles | 20160930 | 20220420 |
| [biostar234081](Biostar234081.md) | convert extended CIGAR to regular CIGAR ('X','=' -> 'M') | 20170130 | 20200409 |
| [biostar234230](Biostar234230.md) | Sliding Window : discriminate partial and fully contained fragments (from a bam file) |  | 20190417 |
| [biostar251649](Biostar251649.md) | Annotating the flanking bases of SNPs in a VCF file | 20170508 | 20200213 |
| [biostar322664](Biostar322664.md) | Extract PE Reads (with their mates) supporting variants in vcf file | 20180625 | 20250313 |
| [biostar332826](Biostar332826.md) | Fast Extraction of Variants from a list of IDs | 20180817 | 20210412 |
| [biostar336589](Biostar336589.md) | displays circular map as SVG from BED and REF file | 20180907 | 20210818 |
| [biostar352930](Biostar352930.md) | Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information. |  |  |
| [biostar398854](Biostar398854.md) | Extract every CDS sequences from a VCF file | 20190916 | 20240418 |
| [biostar404363](Biostar404363.md) | introduce artificial mutation SNV in bam | 20191023 | 20191024 |
| [biostar480685](Biostar480685.md) | paired-end bam clip bases outside insert range | 20201223 | 20200220 |
| [biostar489074](Biostar489074.md) | call variants for every paired overlaping read | 20200205 | 20210412 |
| [biostar497922](Biostar497922.md) | Split VCF into separate VCFs by SNP count | 20210319 | 20210319 |
| [biostar59647](Biostar59647.md) | SAM/BAM to XML | 20131112 | 20190926 |
| [biostar76892](Biostar76892.md) | fix strand of two paired reads close but on the same strand. |  |  |
| [biostar77288](Biostar77288.md) | Low resolution sequence alignment visualization |  | 20240729 |
| [biostar77828](Biostar77828.md) | Divide the human genome among X cores, taking into account gaps |  |  |
| [biostar78285](Biostar78285.md) | Extract BAMs coverage as a VCF file. |  |  |
| [biostar81455](Biostar81455.md) | Defining precisely the exonic genomic context based on a position . | 20130918 | 20200603 |
| [biostar84452](Biostar84452.md) | remove clipped bases from a BAM file |  | 20250408 |
| [biostar84786](Biostar84786.md) | Matrix transposition |  | 20250408 |
| [biostar86363](Biostar86363.md) | Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/ |  |  |
| [biostar86480](Biostar86480.md) | Genomic restriction finder | 20131114 | 20220426 |
| [biostar90204](Biostar90204.md) | Bam version of linux split. |  | 20250408 |
| [biostar9462889](Biostar9462889.md) | Extracting reads from a regular expression in a bam file | 20210402 | 20210402 |
| [biostar9469733](Biostar9469733.md) | Extract reads mapped within chosen intronic region from BAM file | 20210511 | 20210511 |
| [biostar9501110](Biostar9501110.md) | Keep reads including/excluding variants from VCF | 20211210 | 20211213 |
| [biostar9556602](Biostar9556602.md) | Filtering of tricky overlapping sites in VCF |  |  |

### Deprecated/barely used

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [addlinearindextobed](AddLinearIndexToBed.md) | Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart. | 20140201 | 20250115 |
| [bam2sql](BamToSql.md) | Convert a SAM/BAM to sqlite statements | 20160414 | 20160414 |
| [bam2xml](Bam2Xml.md) | converts a BAM to XML | 20130506 | 20210315 |

### Pubmed

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [pubmed404](Pubmed404.md) | Test if URL in the pubmed abstracts are reacheable. | 20181210 | 20200204 |
| [pubmedcodinglang](PubmedCodingLanguages.md) | Programming language use distribution from recent programs / articles | 20170404 | 20200223 |
| [pubmeddump](PubmedDump.md) | Dump XML results from pubmed/Eutils | 20140805 | 20200204 |
| [pubmedgender](PubmedGender.md) | Add gender-related attributes in the Author tag of pubmed xml. |  |  |
| [pubmedgraph](PubmedGraph.md) | Creates a Gephi-gexf graph of references-cites for a given PMID | 20150605 | 20200220 |

### GTF/GFF Manipulation

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [gtf2bed](GtfToBed.md) | Convert GTF/GFF3 to BED. | 20220629 | 20220630 |
| [gtf2xml](Gtf2Xml.md) | Convert GTF/GFF to XML | 20150811 | 20230512 |

### Utilities

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [goutils](GoUtils.md) | Gene Ontology Utils. Retrieves terms from Gene Ontology | 20180130 | 20240523 |
| [multiqcpostproc](MultiqcPostProcessor.md) | Enhances multiqc output by reading the data folder and producing new plots (eg. boxplot per population. | 20240708 | 20240730 |
| [ncbitaxonomy2xml](NcbiTaxonomyToXml.md) | Dump NCBI taxonomy tree as a hierarchical XML document or as a table | 20120320 | 20240320 |
| [oboutils](OboUtils.md) | OBO Ontology Utils. | 20230105 | 20230105 |
| [ukbiobanksamples](UKBiobankSelectSamples.md) | Select samples from ukbiobank | 20210705 | 20220322 |
| [uniprot2svg](UniprotToSvg.md) | plot uniprot to SVG | 20220608 | 20220922 |
| [xsltstream](XsltStream.md) | XSLT transformation for large XML files. xslt is only applied on a given subset of nodes. |  | 20190222 |

### Unclassfied

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [applyvelocity](ApplyVelocity.md) | Execute apache velocity macros | 20241023 | 20241023 |
| [bamclip2insertion](BamClipToInsertion.md) | Convert SOFT clip to Insertion if other read confirm it |  |  |
| [bamcmpcoverage](BamCmpCoverage.md) | Creates the figure of a comparative view of the depths sample vs sample. Memory consideration: the tool alloc an array of bits which size is: (MIN(maxdepth-mindepth,pixel_width_for_one_sample) * count_samples)^2 |  |  |
| [bamliftover](BamLiftOver.md) | Lift-over a BAM file. |  | 20250114 |
| [bamslicebed](BamSliceBed.md) | For @wouter_decoster : slice (long reads) overlapping the records of a BED file | 20191030 | 20210615 |
| [bamtile](BamTile.md) | Answer to @sjackman : Is there a bedtools command to determine a minimal tiling path? A minimal set of features that cover a maximum of the target. |  | 20191010 |
| [barcodegenerator](BarcodeGenerator.md) | Barcode generator for EricCharp | 20230629 | 20230629 |
| [batchigvpictures](BatchIGVPictures.md) | Takes IGV pictures in batch. Save as HTML+png image | 20140731 | 20220524 |
| [bedindextabix](BedIndexTabix.md) | Index and sort a Bed on the fly with Tabix (deprecated). | 20150708 | 20240724 |
| [bedliftover](BedLiftOver.md) | LiftOver a BED file | 20140311 | 20250526 |
| [bedremovebed](BedRemoveBed.md) | Remove bed file from each record of input bed file. Output is a SETFILE | 20221210 | 20221210 |
| [bigwigmerge](BigwigMerge.md) | merge several Bigwig files using different descriptive statistics (mean, median, etc..) | 20240417 | 20240417 |
| [bigwigtview](BigWigTView.md) | view bigwig file coverage in a terminal | 20240704 | 20240704 |
| [bioalcidae](BioAlcidae.md) | javascript version of awk for bioinformatics |  | 20250328 |
| [biostar160470](Biostar160470.md) | Getting untranslated nucleotide sequences on tblastn standalone |  | 20240701 |
| [biostar3654](Biostar3654.md) | show blast alignment with annotations |  |  |
| [biostar95652](Biostar95652.md) | Drawing a schematic genomic context tree. |  |  |
| [biostar9608448](Biostar9608448.md) | Convert long reads to short paired reads | 20250130 | 20250130 |
| [blast2sam](BlastToSam.md) | Convert a **BLASTN-XML** input to SAM |  |  |
| [blastfilterjs](BlastFilterJS.md) | Filters a BlastOutput with a javascript expression. The script injects each <Hit> as the variable 'blasthit'. The user script should return 'true' to keep the hit. |  |  |
| [blastmapannots](BlastMapAnnotations.md) | Maps uniprot/genbank annotations on a blast result. |  |  |
| [blastn2snp](BlastNToSnp.md) | print indel/mismatch in a blastn stream | 20131126 | 20240701 |
| [cmpbams](CompareBams.md) | Compare two or more BAM files | 20130506 | 20200221 |
| [cmpbams4](CompareBams4.md) | Compare two query-name sorted BAM files. Print a tab-delimited report | 20161206 | 20250115 |
| [cmpbamsandbuild](CompareBamAndBuild.md) | Compare two  BAM files mapped on two different builds. Requires a liftover chain file | 20140307 | 20250115 |
| [cnvvalidatorserver](CNVValidatorServer.md) | Review files generated by coverageplotter | 20220818 | 20220826 |
| [commbams](CommBams.md) | Equivalent of unix 'comm' for bams sorted on queryname | 20170420 |  |
| [convertliftoverchain](ConvertLiftOverChain.md) | Convert the contigs in a liftover chain to match another REFerence. (eg. to remove chr prefix, unknown chromosomes etc...) | 20190409 | 20250114 |
| [copynumber01](CopyNumber01.md) | experimental CNV detection. | 20140201 | 20210318 |
| [coverageserver](CoverageServer.md) | Jetty Based http server serving Bam coverage. | 20200212 | 20200330 |
| [cytoband2svg](CytobandToSvg.md) | Creates a svg karyotype . |  |  |
| [dict2bed](DictToBed.md) | convert a SAM dictionary from vcf,sam,bam,dict, etc.. to bed. | 20240603 | 20240603 |
| [dict2xml](DictToXml.md) | convert a SAM dictionary from vcf,sam,bam,dict, etc.. to XML. | 20240824 | 20240824 |
| [dragenbnd2inv](DragenBndToInversion.md) | Converts Dragen BND to inversions | 20241016 | 20241016 |
| [evadumpfiles](EVADumpFiles.md) | Dump files locations from European Variation Archive | 20230314 | 20230314 |
| [extendrefwithreads](ExtendReferenceWithReads.md) | Extending ends of sequences with the help of reads |  | 20190926 |
| [fasta2vcf](ReferenceToVCF.md) | Creates a VCF containing all the possible substitutions from a Reference Genome. | 20140910 | 20240711 |
| [fastqrevcomp](FastqRevComp.md) | produces a reverse-complement fastq (for mate pair alignment see http://seqanswers.com/forums/showthread.php?t=5085 ) |  |  |
| [fastqshuffle](FastqShuffle.md) | Shuffle Fastq files | 20140901 | 20240129 |
| [findhtsfiledict](FindHtsFileDictionary.md) | Scan a set of HTS files (VCF, BAM, CRAM, BCF, etc...), return a tab delimited file (path-of-file,path/url-to-fasta) | 20190912 | 20240824 |
| [fixvcfmissinggenotypes](FixVcfMissingGenotypes.md) | After a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then a missing genotype is said hom-ref. | 20141109 | 20200525 |
| [gatkhc](GatkHaplotypeCaller.md) | Wrapper for GATK HaplotypeCaller | 20240625 | 20240625 |
| [gff2fasta](Gff3ToFasta.md) | extract fasta from gtf | 20241016 | 20241017 |
| [gff2kg](Gff2KnownGene.md) | Convert GFF3 format to UCSC knownGene format. | 20160404 | 20250328 |
| [gff3upstreamorf](Gff3UpstreamOrf.md) | Takes a ucsc genpred file, scan the 5' UTRs and generate a GFF3 containing upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs | 20220724 | 20250131 |
| [gtexrs2qtl](GtexRsToQTL.md) | extract gtex eqtl data from a list of RS | 20230215 | 20240225 |
| [gtf2gff](GtfToGff.md) | Convert GTF to gff | 20220703 | 20250328 |
| [gtfliftover](GtfLiftOver.md) | LiftOver GTF file. | 20190823 | 20250115 |
| [gtfretrocopy](GtfRetroCopy.md) | Scan retrocopies by comparing the gtf/intron and the deletions in a VCF | 20190813 | 20191104 |
| [haplogroupcasectrl](HaploGroupCaseControl.md) | Run Fisher test for Haplogroup input. | 20240610 | 20240610 |
| [howmanybamdict](HowManyBamDict.md) | finds if there's are some differences in the sequence dictionaries. | 20131108 | 20201021 |
| [htsfreemarker](HtsVelocity.md) | Apply apache velocity to VCF/BAM/JSON files. | 20230616 | 20230616 |
| [illuminadir](IlluminaDirectory.md) | Create a structured (**JSON** or **XML**) representation of a directory containing some Illumina FASTQs. | 20131021 | 20180717 |
| [jbrowse2](JBrowse2Server.md) | create a run a local instance of jbrowse2 | 20250404 | 20250404 |
| [kg2bed](KnownGenesToBed.md) | converts UCSC knownGenes file to BED. | 20140311 | 20230815 |
| [kg2fa](KnownGeneToFasta.md) | convert ucsc genpred to fasta | 20190213 | 20230815 |
| [kg2gff](KgToGff.md) | Convert UCSC genpred/knowngene file to gff3 or gtf | 20210106 | 20250324 |
| [knownretrocopy](KnownRetroCopy.md) | Annotate VCF structural variants that could be intron from retrocopies. | 20190815 | 20230817 |
| [liftover2svg](LiftOverToSVG.md) | Convert LiftOver chain files to animated SVG |  |  |
| [manhattan](Manhattan.md) | Manhattan plot SVG picture from different sources. | 20220525 | 20240525 |
| [mapuniprot](MapUniProtFeatures.md) | map uniprot features on reference genome. |  | 20250331 |
| [mergeblastxml](MergeBlastXml.md) | merge XML blast results (same Iteration/Iteration_query-def in multiple xml files |  |  |
| [mergesplittedblast](MergeSplittedBlast.md) | merge blast Hits from splitted BLAST database |  |  |
| [ngsfilessummary](NgsFilesSummary.md) | Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..). Useful to get a summary of your samples. | 20140430 | 20240324 |
| [optimizefisher](OptimizeFisher.md) | Optimize fisher test on VCF using genetic algo | 20221013 | 20240207 |
| [pcrclipreads](PcrClipReads.md) | Soft clip bam files based on PCR target regions | 20150618 | 20210322 |
| [pcrslicereads](PcrSliceReads.md) | Mark PCR reads to their PCR amplicon | 20150707 | 20240724 |
| [plotbcftoolsstats](PlotBcftoolsStats.md) | Plot bcftools stats output | 20210622 | 20210622 |
| [pubmedauthorgraph](PubmedAuthorGraph.md) | Creates a graph from Pubmed and Authors |  |  |
| [pubmedmap](PubmedMap.md) | Use Pubmed Author's Affiliation to map the authors in the world. | 20160426 |  |
| [pubmedorcidgraph](PubmedOrcidGraph.md) | Creates a graph from Pubmed and Authors' Orcid identifiers | 20160520 | 20210712 |
| [qqplotter](QQPlotter.md) | plot QQplot | 20250324 | 20250324 |
| [reduceblast](ReduceBlast.md) | Reduce the size of XML blast, by removing iterations that have no Hit |  |  |
| [regeniebedannot](RegenieBedAnnot.md) | Create annotation files for regenie using sliding annotations | 20250311 | 202050515 |
| [regeniefunctionalannot](RegenieFunctionalAnnot.md) | Create annotation files for regenie using snpEff annotations | 20250311 | 20250320 |
| [regeniemakeannot](RegenieMakeAnnot.md) | Create annotation files for regenie from a TSV input file | 20250311 | 20250320 |
| [regenieslidingannot](RegenieSlidingAnnot.md) | Create annotation files for regenie using sliding annotations | 20250311 | 20250320 |
| [repairfastq](RepairFastq.md) | Join single end reads to paired end | 20240128 | 20240128 |
| [rnaseqpolya](RNASeqPolyA.md) | find poly-A tail in RNASeq data | 20210913 | 20210914 |
| [sam2json](SamToJson.md) | Convert a SAM input to JSON | 20210402 | 20210315 |
| [sam4weblogo](SAM4WebLogo.md) | Sequence logo for different alleles or generated from SAM/BAM | 20130524 | 20250326 |
| [samaddpi](SamAddPI.md) | Add predicted median insert size 'PI' to SAM Read groups (RG). |  |  |
| [samedict](SameDict.md) | check if all HTS files share the same dictionary | 20240724 | 20240724 |
| [samfixcigar](SamFixCigar.md) | Fix Cigar String in SAM replacing 'M' by 'X' or '=' | 20131126 | 20210223 |
| [samjdk](SamJdk.md) | Filters a BAM using a java expression compiled in memory. | 20170807 | 20191119 |
| [samslop](SamSlop.md) | extends sam by 'x' bases using the reference sequence | 20160119 | 20250522 |
| [scanlabguru](ScanLabGuru.md) | scan the files stored in labguru | 20240325 | 20240325 |
| [scansv](ScanStructuralVariants.md) | Scan structural variants for case/controls data | 20190815 | 20240916 |
| [shiftbam](ShiftBam.md) | shit all coordinates of a bam | 20241001 | 20241001 |
| [shiftvcf](ShiftVcf.md) | shit all coordinates of a VCF | 20241002 | 20241002 |
| [sortvcfoninfo](SortVcfOnInfo.md) | Sort a VCF a field in the INFO column | 20140218 | 20201204 |
| [sv2fasta](StructuralVariantToFasta.md) | convert VCF of structural variant(s) to fasta for pggb | 20230403 | 20230403 |
| [swingplinkselectcluster](SwingPLinkSelectCluster.md) | Swing-based Plink/MDS sample selector | 20231123 | 20240606 |
| [swingregenie](RegenieSwing.md) | view regenie output | 20250324 | 20250324 |
| [translategff3](TranslateGff3.md) | translates the output of bcftools consensus | 20241003 | 20241003 |
| [tssenrich](TSSEnrichment.md) | Transcription Start Site (TSS) Enrichment Score calculation | 20240130 | 20240206 |
| [ukbbdump](UkbiobankDump.md) | Dump data https://afb.ukbiobank.ac.uk/ ukbiobank server. The server might not like too many requests. Use a your own risk. Doesn't work with jdk17 (?!) | 20250425 | 20250425 |
| [uniprotfilterjs](UniprotFilterJS.md) | Filters Uniprot DUMP+ XML with a javascript  (java rhino) expression. Context contain 'entry' an uniprot entry and 'index', the index in the XML file. |  | 20250328 |
| [variantsinwindow](VariantsInWindow.md) | Annotate Number of Variants overlaping a sliding window. |  |  |
| [vcf2bam](VcfToBam.md) | vcf to bam | 20150612 | 20211022 |
| [vcf2hilbert](VcfToHilbert.md) | Plot a Hilbert Curve from a VCF file as SVG | 20171201 | 20240517 |
| [vcf2xml](Vcf2Xml.md) | Convert VCF to XML |  | 20230822 |
| [vcfburdencnv](VcfBurdenCNV.md) | Burden on CNV (experimental) | 20250404 | 20250408 |
| [vcfburdenfisherh](VcfBurdenFisherH.md) | Fisher Case /Controls per Variant | 20160418 | 20200713 |
| [vcfburdenslidingwindow](VcfBurdenSlidingWindow.md) | apply fisher test on VCF using a sliding window | 20190920 | 20231213 |
| [vcffilterbyliftover](VcfFilterByLiftOver.md) | Add FILTER(s) to a variant when it is known to map elsewhere after liftover. | 20190418 | 20210603 |
| [vcfgatkeval](VcfGatkEval.md) | Eval/Plot gatk INFO tags for filtering | 20230424 | 20240321 |
| [vcfgroupbypop](VcfGroupByPopulation.md) | create INFO data by population | 20190319 | 20230712 |
| [vcfpeekvcf](VcfPeekVcf.md) | Get the INFO from a VCF and use it for another VCF | 20150521 | 20240405 |
| [vcfsamplesprs](VcfSamplesPRS.md) | another program for @AntoineRimbert | 20220915 | 20220915 |
| [vcfscanupstreamorf](VcfScanUpstreamOrf.md) | Scan BAM for upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs | 20190218 | 20200804 |
| [vcfserver](VcfServer.md) | Web Server displaying VCF file. A web interface for vcf2table | 20171027 | 20220517 |
| [vcfspliceai](VcfSpliceAI.md) | Annotate VCF with spiceai web service | 20201107 | 20201107 |
| [vcftbi2bed](VcfTbiToBed.md) | extracts BED for each contig in a tabix-indexed VCF peeking first of last variant for each chromosome. | 20230214 | 20230214 |
| [vcfukbb](VcfUkbiobank.md) | annotates an VCF with the https://afb.ukbiobank.ac.uk/ ukbiobank server. The server might not like too many requests. Use a your own risk. Doesn't work with jdk17 (?!) | 20250424 | 20250424 |
| [wib2bedgraph](WibToBedGraph.md) | Extract Wib files to bedgraph or wig | 20230819 | 20230819 |
| [xcontaminations](XContaminations.md) | For @AdrienLeger2 : cross contamination between samples by looking at the homozygous genotypes. |  |  |

### VCF Manipulation

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [bcftoolsmergebest](BCFToolsMergeBest.md) | Scan a VCF file generated by 'bcftools merge --force-samples', identify duplicate samples, keep the best | 20240604 | 20240604 |
| [bed2vcf](BedToVcf.md) | Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele | 20240604 | 20240604 |
| [bioalcidaejdk](BioAlcidaeJdk.md) | java-based version of awk for bioinformatics | 20170712 | 20210412 |
| [biostar130456](Biostar130456.md) | Split individual VCF files from multisamples VCF file | 20150210 | 20200603 |
| [breakdancer2vcf](BreakdancerToVcf.md) | Convert output of breakdancer to VCF | 20200511 | 20241122 |
| [builddbsnp](BuildDbsnp.md) | Build a DBSNP file from different sources for GATK | 20200904 | 2021070726 |
| [findavariation](FindAVariation.md) | Finds a specific mutation in a list of VCF files | 20140623 | 20200217 |
| [findgvcfsblocks](FindGVCFsBlocks.md) | Find common blocks of calleable regions from a set of gvcfs | 20210806 | 20220401 |
| [mantamerger](MantaMerger.md) | Merge Vcf from Manta VCF. | 20190916 | 20230320 |
| [minicaller](MiniCaller.md) | Simple and Stupid Variant Caller designed for @AdrienLeger2 | 201500306 | 20250327 |
| [svcasescontrols](SVCasesControls.md) | Find SV present in cases but not in controls. | 20240513 | 20240516 |
| [swingvcfjexl](SwingVcfJexlFilter.md) | Filter VCF using Java Swing UI and JEXL/Javascript expression | 20220413 | 20220414 |
| [swingvcfview](SwingVcfView.md) | VCFviewer using Java Swing UI | 20210503 | 20210503 |
| [vcf2intervals](VcfToIntervals.md) | split a vcf to interval or bed for parallelization | 20211112 | 20221128 |
| [vcf2r](VcfToRScript.md) | Convert VCF to R so it can be used for burden testing |  | 20240607 |
| [vcf2rdf](VcfToRdf.md) | convert VCF to RDF (N3 notation) | 20191213 | 20250520 |
| [vcf2table](VcfToTable.md) | convert a vcf to a table, to ease display in the terminal | 20170511 | 20250523 |
| [vcfallelebalance](VcfAlleleBalance.md) | Insert missing allele balance annotation using FORMAT:AD | 20180829 | 20200805 |
| [vcfancestralalleles](VcfAncestralAllele.md) | Annotate a VCF with it's ancestral allele. Data from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README | 20180418 | 20220126 |
| [vcfbigbed](VcfBigBed.md) | Annotate a VCF with values from a bigbed file | 20220107 | 20220107 |
| [vcfbigwig](VCFBigWig.md) | Annotate a VCF with values from a bigwig file | 20200506 | 20230819 |
| [vcfbraiding](VcfBraiding.md) | visualization for variants and attributes using https://visdunneright.github.io/sequence_braiding/docs/ . | 20201021 | 20201022 |
| [vcfburdenmaf](VcfBurdenMAF.md) | MAF for Cases / Controls | 20160418 | 202000713 |
| [vcfcadd](VcfCadd.md) | Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892.PubMed PMID: 24487276. | 20220119 | 20240524 |
| [vcfcombinetwosnvs](VCFCombineTwoSnvs.md) | Detect Mutations than are the consequences of two distinct variants. This kind of variant might be ignored/skipped from classical variant consequence predictor. Idea from @SolenaLS and then @AntoineRimbert | 20160215 | 20200425 |
| [vcfcomposite](VCFComposite.md) | (in developpement) Finds Variants involved in a Het Compound Disease | 20170331 | 20200210 |
| [vcfconcat](VcfConcat.md) | Concatenate VCFs with same sample. See also bcftools concat | 20131230 | 20240426 |
| [vcfdistancevariants](VcfDistanceBetweenVariants.md) | Annotate variants with the distance between previous and next variant. | 20190410 | 20230510 |
| [vcffiltergenes](VcFilterGenes.md) | Filter VEP/SnpEff Output from a list of genes. | 20160322 | 20230505 |
| [vcffiltergtf](VcfFilterGtf.md) | Filter VCF on GTF | 20230703 | 20230704 |
| [vcffilterjdk](VcfFilterJdk.md) | Filtering VCF with dynamically-compiled java expressions | 20170705 | 20220830 |
| [vcffilterso](VcfFilterSequenceOntology.md) | Filter a VCF file annotated with SNPEff or VEP with terms from Sequence-Ontology. Reasoning : Children of user's SO-terms will be also used. | 20170331 | 20200924 |
| [vcfflatten](VCFFlatten.md) | Flatten variants to one variant | 20230222 | 20230222 |
| [vcfgenesplitter](VcfGeneSplitter.md) | Split VCF+VEP by gene/transcript. | 20160310 | 202220531 |
| [vcfgnomad](VcfGnomad.md) | Peek annotations from gnomad | 20170407 | 20231103 |
| [vcfgnomadsv](VcfGnomadSV.md) | Peek annotations from gnomad structural variants | 20190814 | 20211109 |
| [vcfgrantham](VcfGrantham.md) | add grantham score from annotated VCF variant | 20230503 | 20230503 |
| [vcfhead](VcfHead.md) | print the first variants of a vcf | 202050510 | 20200518 |
| [vcfliftover](VcfLiftOver.md) | Lift-over a VCF file | 20240114 | 20210603 |
| [vcfmovefilterstoinfo](VcfMoveFiltersToInfo.md) | Move any FILTER to the INFO column. reset FILTER to PASS | 20161025 | 20220323 |
| [vcfmulti2one](VcfMultiToOne.md) | Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column. Never used. | 20150312 | 20241125 |
| [vcfmulti2oneinfo](VcfMultiToOneInfo.md) | 'one variant with INFO with N values' to 'N variants with one INFO' | 20260106 | 20230524 |
| [vcfnearest](VCFNearest.md) | find nearest feature near a variant | 20250523 | 20250523 |
| [vcfpar](VcfPseudoAutosomalRegion.md) | Flag human sexual regions excluding PAR. | 20200908 | 20200908 |
| [vcfpeekaf](VcfPeekAf.md) | Peek the AF from another VCF | 20200624 | 20200904 |
| [vcfphased01](VcfPhased01.md) | X10 Phased SVG to Scalar Vector Graphics (SVG) | 20190710 | 20230818 |
| [vcfpolyx](VCFPolyX.md) | Number of repeated REF bases around POS. | 20200930 | 20230526 |
| [vcfrebase](VcfRebase.md) | Restriction sites overlaping variations in a vcf | 20131115 | 20200624 |
| [vcfregulomedb](VcfRegulomeDB.md) | Annotate a VCF with the Regulome2 data (https://regulomedb.org/) | 20140709 | 20230512 |
| [vcfsetdict](VcfSetSequenceDictionary.md) | Set the `##contig` lines in a VCF header on the fly | 20140105 | 20210201 |
| [vcfshuffle](VCFShuffle.md) | Shuffle a VCF | 20131210 | 20250523 |
| [vcfsplitnvariants](VcfSplitNVariants.md) | Split VCF to 'N' VCF files | 202221122 | 202221201 |
| [vcfsplitvep](VCFSplitVEP.md) | Split CSQ vep annotations | 20250517 | 20250517 |
| [vcfspringfilter](VcfSpringFilter.md) | Uses the java spring Framework to build complex vcf filters | 20230526 | 20230526 |
| [vcfstats](VcfStats.md) | Produce VCF statitics | 20131212 | 20230707 |
| [vcfsvannotator](VCFSVAnnotator.md) | SV Variant Effect prediction using gtf, gnomad, etc | 20190815 | 20230509 |
| [vcftail](VcfTail.md) | print the last variants of a vcf | 20131210 | 20200518 |
| [vcftrio](VCFTrios.md) | Find mendelian incompatibilitie / denovo variants in a VCF | 20130705 | 20200624 |

### Retrocopy

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [scanretrocopy](ScanRetroCopy.md) | Scan BAM for retrocopies | 20190125 | 20230818 |
| [starretrocopy](StarRetroCopy.md) | Scan retrocopies from the star-aligner/bwa output | 20190710 | 20191008 |

### BAM Manipulation

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [bam2haplotypes](BamToHaplotypes.md) | Reconstruct SNP haplotypes from reads | 20211015 | 20211020 |
| [bamleftalign](BamLeftAlign.md) | Left Align Reads around deletions | 20250327 | 20250327 |
| [bamphased01](BamPhased01.md) | Extract Reads from a SAM/BAM file supporting at least two variants in a VCF file. | 20210218 | 20210218 |
| [bamrenamechr](ConvertBamChromosomes.md) | Convert the names of the chromosomes in a BAM file | 20131217 | 20191210 |
| [bamstats04](BamStats04.md) | Coverage statistics for a BED file. | 20130513 | 20191003 |
| [bamstats05](BamStats05.md) | Coverage statistics for a BED file, group by gene | 20151012 | 20210317 |
| [bamwithoutbai](BamWithoutBai.md) | Query a Remote BAM without bai | 20191213 | 20191217 |
| [basecoverage](BaseCoverage.md) | 'Depth of Coverage' per base. | 20220420 | 20241122 |
| [bioalcidaejdk](BioAlcidaeJdk.md) | java-based version of awk for bioinformatics | 20170712 | 20210412 |
| [biostar154220](Biostar154220.md) | Cap BAM to a given coverage | 20150812 | 20210312 |
| [biostar9566948](Biostar9566948.md) | Trim Reads So Only First Base Remains | 20230621 | 20230621 |
| [findallcoverageatposition](FindAllCoverageAtPosition.md) | Find depth at specific position in a list of BAM files. My colleague Estelle asked: in all the BAM we sequenced, can you give me the depth at a given position ? | 20141128 | 20250327 |
| [sam2tsv](Sam2Tsv.md) | Prints the SAM alignments as a TAB delimited file. | 20170712 | 20250327 |
| [samextractclip](SamExtractClip.md) | Extract Soft Clipped Sequences from a SAM. Ouput is a FASTQ | 20140228 | 20240524 |
| [samgrep](SamGrep.md) | grep read-names in a bam file | 20130506 | 20210726 |
| [samrmdupnames](SamRemoveDuplicatedNames.md) | remove duplicated names in sorted BAM | 20240405 | 20221207 |
| [samviewwithmate](SamViewWithMate.md) | Extract reads within given region(s), and their mates | 20190207 | 20250219 |
| [sortsamrefname](SortSamRefName.md) | Sort a BAM on chromosome/contig and then on read/querty name | 20150812 | 20210312 |
| [swingbamcov](SwingBamCov.md) | Bam coverage viewer using Java Swing UI | 20210420 | 20220513 |
| [swingbamview](SwingBamView.md) | Read viewer using Java Swing UI | 20220503 | 20230427 |
| [texbam](TextBam.md) | Write text in a bam. Mostly for fun... | 20220708 | 20220708 |
| [tview](TViewCmd.md) | equivalent of samtools tview | 20130613 | 20241209 |


