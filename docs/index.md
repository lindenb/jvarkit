JVARKIT
=======

Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.
Version     : 1a9ec1562
Compilation : 20230127152508
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

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| [addlinearindextobed](AddLinearIndexToBed.md) | Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart. | 20140201 | 20230126 |
| [backlocate](BackLocate.md) | Mapping a mutation on a protein back to the genome. | 20140619 | 20190820 |
| [bam2raster](Bam2Raster.md) | BAM to raster graphics |  |  |
| [bam2sql](BamToSql.md) | Convert a SAM/BAM to sqlite statements |  |  |
| [bam2svg](BamToSVG.md) | BAM to Scalar Vector Graphics (SVG) | 20141013 | 20210728 |
| [bam2xml](Bam2Xml.md) | converts a BAM to XML | 20130506 | 20210315 |
| [bammatrix](BamMatrix.md) | Bam matrix, inspired from 10x/loupe | 20190620 | 20211206 |
| [bamphased01](BamPhased01.md) | Extract Reads from a SAM/BAM file supporting at least two variants in a VCF file. | 20210218 | 20210218 |
| [bamrenamechr](ConvertBamChromosomes.md) | Convert the names of the chromosomes in a BAM file | 20131217 | 20191210 |
| [bamstats05](BamStats05.md) | Coverage statistics for a BED file, group by gene | 20151012 | 20210317 |
| [bamwithoutbai](BamWithoutBai.md) | Query a Remote BAM without bai | 20191213 | 20191217 |
| [basecoverage](BaseCoverage.md) | 'Depth of Coverage' per base. | 20220420 | 20220420 |
| [bedcluster](BedCluster.md) | Clusters a BED file into a set of BED files. | 20200130 | 20220914 |
| [bedmergecnv](BedMergeCnv.md) | Merge Bed records if they overlap a fraction of their lengths. | 20200330 | 20200603 |
| [bednonoverlappingset](BedNonOverlappingSet.md) | Split a Bed file into non-overlapping data set. | 20180607 | 20200408 |
| [bedrenamechr](ConvertBedChromosomes.md) | Convert the names of the chromosomes in a Bed file |  | 20190503 |
| [bioalcidaejdk](BioAlcidaeJdk.md) | java-based version of awk for bioinformatics | 20170712 | 20210412 |
| [biostar103303](Biostar103303.md) | Calculate Percent Spliced In (PSI). |  |  |
| [biostar105754](Biostar105754.md) | bigwig : peak distance from specific genomic region | 20140708 | 20220110 |
| [biostar130456](Biostar130456.md) | Split individual VCF files from multisamples VCF file | 20150210 | 20200603 |
| [biostar139647](Biostar139647.md) | Convert alignment in Fasta/Clustal format to SAM/BAM file |  |  |
| [biostar145820](Biostar145820.md) | subsample/shuffle BAM to fixed number of alignments. | 20150615 | 20211005 |
| [biostar154220](Biostar154220.md) | Cap BAM to a given coverage | 20150812 | 20210312 |
| [biostar165777](Biostar165777.md) | Split a XML file |  |  |
| [biostar170742](Biostar170742.md) | convert sam format to axt Format | 20151228 | 20210412 |
| [biostar172515](Biostar172515.md) | Convert BAI to XML |  |  |
| [biostar173114](Biostar173114.md) | make a bam file smaller by removing unwanted information see also https://www.biostars.org/p/173114/ |  |  |
| [biostar175929](Biostar175929.md) | Construct a combination set of fasta sequences from a vcf | 20160208 | 20211012 |
| [biostar178713](Biostar178713.md) | split bed file into several bed files where each region is separated of any other by N bases | 20160226 | 20200818 |
| [biostar214299](Biostar214299.md) | Extract allele specific reads from bamfiles | 20160930 | 20220420 |
| [biostar234081](Biostar234081.md) | convert extended CIGAR to regular CIGAR ('X','=' -> 'M') | 20170130 | 20200409 |
| [biostar234230](Biostar234230.md) | Sliding Window : discriminate partial and fully contained fragments (from a bam file) |  | 20190417 |
| [biostar251649](Biostar251649.md) | Annotating the flanking bases of SNPs in a VCF file | 20170508 | 20200213 |
| [biostar322664](Biostar322664.md) | Extract PE Reads (with their mates) supporting variants in vcf file |  |  |
| [biostar332826](Biostar332826.md) | Fast Extraction of Variants from a list of IDs | 20180817 | 20210412 |
| [biostar336589](Biostar336589.md) | displays circular map as SVG from BED and REF file | 20180907 | 20210818 |
| [biostar352930](Biostar352930.md) | Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information. |  |  |
| [biostar398854](Biostar398854.md) | Extract every CDS sequences from a VCF file | 20190916 | 20190916 |
| [biostar404363](Biostar404363.md) | introduce artificial mutation SNV in bam | 20191023 | 20191024 |
| [biostar480685](Biostar480685.md) | paired-end bam clip bases outside insert range | 20201223 | 20200220 |
| [biostar489074](Biostar489074.md) | call variants for every paired overlaping read | 20200205 | 20210412 |
| [biostar497922](Biostar497922.md) | Split VCF into separate VCFs by SNP count | 20210319 | 20210319 |
| [biostar59647](Biostar59647.md) | SAM/BAM to XML | 20131112 | 20190926 |
| [biostar76892](Biostar76892.md) | fix strand of two paired reads close but on the same strand. |  |  |
| [biostar77288](Biostar77288.md) | Low resolution sequence alignment visualization |  |  |
| [biostar77828](Biostar77828.md) | Divide the human genome among X cores, taking into account gaps |  |  |
| [biostar78285](Biostar78285.md) | Extract BAMs coverage as a VCF file. |  |  |
| [biostar81455](Biostar81455.md) | Defining precisely the exonic genomic context based on a position . | 20130918 | 20200603 |
| [biostar84452](Biostar84452.md) | remove clipped bases from a BAM file |  |  |
| [biostar84786](Biostar84786.md) | Matrix transposition |  |  |
| [biostar86363](Biostar86363.md) | Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/ |  |  |
| [biostar86480](Biostar86480.md) | Genomic restriction finder | 20131114 | 20220426 |
| [biostar90204](Biostar90204.md) | Bam version of linux split. |  |  |
| [biostar9462889](Biostar9462889.md) | Extracting reads from a regular expression in a bam file | 20210402 | 20210402 |
| [biostar9469733](Biostar9469733.md) | Extract reads mapped within chosen intronic region from BAM file | 20210511 | 20210511 |
| [biostar9501110](Biostar9501110.md) | Keep reads including/excluding variants from VCF | 20211210 | 20211213 |
| [builddbsnp](BuildDbsnp.md) | Build a DBSNP file from different sources for GATK | 20200904 | 2021070726 |
| [cnvtview](CnvTView.md) | Text visualization of bam DEPTH for multiple regions in a terminal | 20181018 | 20210412 |
| [coverageplotter](CoveragePlotter.md) | Display an image of depth to display any anomaly an intervals+bams | 20200605 | 20221125 |
| [findallcoverageatposition](FindAllCoverageAtPosition.md) | Find depth at specific position in a list of BAM files. My colleague Estelle asked: in all the BAM we sequenced, can you give me the depth at a given position ? | 20141128 | 20210818 |
| [findavariation](FindAVariation.md) | Finds a specific mutation in a list of VCF files | 20140623 | 20200217 |
| [findgvcfsblocks](FindGVCFsBlocks.md) | Find common blocks of calleable regions from a set of gvcfs | 20210806 | 20220401 |
| [groupbygene](GroupByGene.md) | Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff | 20131209 | 20220529 |
| [gtf2bed](GtfToBed.md) | Convert GTF/GFF3 to BED. | 20220629 | 20220630 |
| [lowresbam2raster](LowResBam2Raster.md) | Low Resolution BAM to raster graphics | 20170523 | 20211126 |
| [minicaller](MiniCaller.md) | Simple and Stupid Variant Caller designed for @AdrienLeger2 | 201500306 | 20220705 |
| [mkminibam](MakeMiniBam.md) | Creates an archive of small bams with only a few regions. | 20190410 | 20221019 |
| [plotsashimi](PlotSashimi.md) | Print Sashimi plots from Bam | 20191117 | 20191104 |
| [prettysam](PrettySam.md) | Pretty SAM alignments | 20171215 | 20211105 |
| [pubmed404](Pubmed404.md) | Test if URL in the pubmed abstracts are reacheable. | 20181210 | 20200204 |
| [pubmedcodinglang](PubmedCodingLanguages.md) | Programming language use distribution from recent programs / articles | 20170404 | 20200223 |
| [pubmedgender](PubmedGender.md) | Add gender-related attributes in the Author tag of pubmed xml. |  |  |
| [pubmedgraph](PubmedGraph.md) | Creates a Gephi-gexf graph of references-cites for a given PMID | 20150605 | 20200220 |
| [sam2tsv](Sam2Tsv.md) | Prints the SAM alignments as a TAB delimited file. | 20170712 | 20210304 |
| [samgrep](SamGrep.md) | grep read-names in a bam file | 20130506 | 20210726 |
| [samrmdupnames](SamRemoveDuplicatedNames.md) | remove duplicated names in sorted BAM | 20221207 | 20221207 |
| [samviewwithmate](SamViewWithMate.md) | Extract reads within given region(s), and their mates | 20190207 | 20191004 |
| [scanretrocopy](ScanRetroCopy.md) | Scan BAM for retrocopies | 20190125 | 20190709 |
| [setfiletools](SetFileTools.md) | Utilities for the setfile format | 20210125 | 20220426 |
| [sortsamrefname](SortSamRefName.md) | Sort a BAM on chromosome/contig and then on read/querty name | 20150812 | 20210312 |
| [swingbamcov](SwingBamCov.md) | Bam coverage viewer using Java Swing UI | 20210420 | 20220513 |
| [swingbamview](SwingBamView.md) | Read viewer using Java Swing UI | 20220503 | 20230124 |
| [swingindexcov](SwingIndexCov.md) | indexcov visualization | 2020511 | 2020512 |
| [swingvcfjexl](SwingVcfJexlFilter.md) | Filter VCF using Java Swing UI and JEXL/Javascript expression | 20220413 | 20220414 |
| [swingvcfview](SwingVcfView.md) | VCFviewer using Java Swing UI | 20210503 | 20210503 |
| [ukbiobanksamples](UKBiobankSelectSamples.md) | Select samples from ukbiobank | 20210705 | 20220322 |
| [uniprot2svg](UniprotToSvg.md) | plot uniprot to SVG | 20220608 | 20220922 |
| [vcf2table](VcfToTable.md) | convert a vcf to a table, to ease display in the terminal | 20170511 | 20220507 |
| [vcfallelebalance](VcfAlleleBalance.md) | Insert missing allele balance annotation using FORMAT:AD | 20180829 | 20200805 |
| [vcfbigbed](VcfBigBed.md) | Annotate a VCF with values from a bigbed file | 20220107 | 20220107 |
| [vcfbigwig](VCFBigWig.md) | Annotate a VCF with values from a bigwig file | 20200506 | 20220110 |
| [vcffilterjdk](VcfFilterJdk.md) | Filtering VCF with dynamically-compiled java expressions | 20170705 | 20220830 |
| [vcffilterso](VcfFilterSequenceOntology.md) | Filter a VCF file annotated with SNPEff or VEP with terms from Sequence-Ontology. Reasoning : Children of user's SO-terms will be also used. | 20170331 | 20200924 |
| [vcfgenesplitter](VcfGeneSplitter.md) | Split VCF+VEP by gene/transcript. | 20160310 | 202220531 |
| [vcfgnomad](VcfGnomad.md) | Peek annotations from gnomad | 20170407 | 20200702 |
| [vcfhead](VcfHead.md) | print the first variants of a vcf | 20131210 | 20200518 |
| [vcfpar](VcfPseudoAutosomalRegion.md) | Flag human sexual regions excluding PAR. | 20200908 | 20200908 |
| [vcfpolyx](VCFPolyX.md) | Number of repeated REF bases around POS. | 20200930 | 20211102 |
| [vcfrebase](VcfRebase.md) | Restriction sites overlaping variations in a vcf | 20131115 | 20200624 |
| [vcfshuffle](VCFShuffle.md) | Shuffle a VCF | 20131210 | 20200818 |
| [vcfsplitnvariants](VcfSplitNVariants.md) | Split VCF to 'N' VCF files | 202221122 | 202221201 |
| [vcfstrech2svg](VcfStrechToSvg.md) | another VCF to SVG | 20210304 | 20210309 |
| [vcftail](VcfTail.md) | print the last variants of a vcf | 20131210 | 20200518 |
| [vcftrio](VCFTrios.md) | Find mendelian incompatibilitie / denovo variants in a VCF | 20130705 | 20200624 |
| [wescnvsvg](WesCnvSvg.md) | SVG visualization of bam DEPTH for multiple regions | 20180726 | 20210726 |
| [wgscoverageplotter](WGSCoveragePlotter.md) | Whole genome coverage plotter | 20201125 | 20210812 |
