JVARKIT
=======

Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.
Version     : 192fbc43
Compilation : 20230126211827
Github      : https://github.com/lindenb/jvarkit
Issues      : https://github.com/lindenb/jvarkit/issues

## Usage

  java -jar jvarkit.jar [options]
or
  java -jar jvarkit.jar <command name> (other arguments)

## Options

 + --help show this screen
 + --help-all show all commands, including the private ones.
 + --version print version

## Compilation Installation

Please, read [how to run and install jvarkit](JvarkitCentral.md)

## Tools

| Tool | Description | Creation | Update |
| ---: | :---------- | :------: | :----: |
| addlinearindextobed | Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart. | 20140201 | 20230126 |
| backlocate | Mapping a mutation on a protein back to the genome. | 20140619 | 20190820 |
| bam2raster | BAM to raster graphics |  |  |
| bam2svg | BAM to Scalar Vector Graphics (SVG) | 20141013 | 20210728 |
| bam2xml | converts a BAM to XML | 20130506 | 20210315 |
| bamrenamechr | Convert the names of the chromosomes in a BAM file | 20131217 | 20191210 |
| bamstats05 | Coverage statistics for a BED file, group by gene | 20151012 | 20210317 |
| basecoverage | 'Depth of Coverage' per base. | 20220420 | 20220420 |
| bedcluster | Clusters a BED file into a set of BED files. | 20200130 | 20220914 |
| bednonoverlappingset | Split a Bed file into non-overlapping data set. | 20180607 | 20200408 |
| bedrenamechr | Convert the names of the chromosomes in a Bed file |  | 20190503 |
| bioalcidaejdk | java-based version of awk for bioinformatics | 20170712 | 20210412 |
| biostar103303 | Calculate Percent Spliced In (PSI). |  |  |
| biostar105754 | bigwig : peak distance from specific genomic region | 20140708 | 20220110 |
| biostar130456 | Split individual VCF files from multisamples VCF file | 20150210 | 20200603 |
| biostar139647 | Convert alignment in Fasta/Clustal format to SAM/BAM file |  |  |
| biostar145820 | subsample/shuffle BAM to fixed number of alignments. | 20150615 | 20211005 |
| biostar154220 | Cap BAM to a given coverage | 20150812 | 20210312 |
| biostar165777 | Split a XML file |  |  |
| biostar170742 | convert sam format to axt Format | 20151228 | 20210412 |
| biostar172515 | Convert BAI to XML |  |  |
| biostar173114 | make a bam file smaller by removing unwanted information see also https://www.biostars.org/p/173114/ |  |  |
| biostar175929 | Construct a combination set of fasta sequences from a vcf | 20160208 | 20211012 |
| biostar178713 | split bed file into several bed files where each region is separated of any other by N bases | 20160226 | 20200818 |
| biostar214299 | Extract allele specific reads from bamfiles | 20160930 | 20220420 |
| biostar234081 | convert extended CIGAR to regular CIGAR ('X','=' -> 'M') | 20170130 | 20200409 |
| biostar234230 | Sliding Window : discriminate partial and fully contained fragments (from a bam file) |  | 20190417 |
| biostar251649 | Annotating the flanking bases of SNPs in a VCF file | 20170508 | 20200213 |
| biostar322664 | Extract PE Reads (with their mates) supporting variants in vcf file |  |  |
| biostar332826 | Fast Extraction of Variants from a list of IDs | 20180817 | 20210412 |
| biostar336589 | displays circular map as SVG from BED and REF file | 20180907 | 20210818 |
| biostar352930 | Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information. |  |  |
| biostar398854 | Extract every CDS sequences from a VCF file | 20190916 | 20190916 |
| biostar404363 | introduce artificial mutation SNV in bam | 20191023 | 20191024 |
| biostar480685 | paired-end bam clip bases outside insert range | 20201223 | 20200220 |
| biostar489074 | call variants for every paired overlaping read | 20200205 | 20210412 |
| biostar497922 | Split VCF into separate VCFs by SNP count | 20210319 | 20210319 |
| biostar59647 | SAM/BAM to XML | 20131112 | 20190926 |
| biostar76892 | fix strand of two paired reads close but on the same strand. |  |  |
| biostar77288 | Low resolution sequence alignment visualization |  |  |
| biostar77828 | Divide the human genome among X cores, taking into account gaps |  |  |
| biostar78285 | Extract BAMs coverage as a VCF file. |  |  |
| biostar81455 | Defining precisely the exonic genomic context based on a position . | 20130918 | 20200603 |
| biostar84452 | remove clipped bases from a BAM file |  |  |
| biostar84786 | Matrix transposition |  |  |
| biostar86363 | Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/ |  |  |
| biostar86480 | Genomic restriction finder | 20131114 | 20220426 |
| biostar90204 | Bam version of linux split. |  |  |
| biostar9462889 | Extracting reads from a regular expression in a bam file | 20210402 | 20210402 |
| biostar9469733 | Extract reads mapped within chosen intronic region from BAM file | 20210511 | 20210511 |
| biostar9501110 | Keep reads including/excluding variants from VCF | 20211210 | 20211213 |
| builddbsnp | Build a DBSNP file from different sources for GATK | 20200904 | 2021070726 |
| cnvtview | Text visualization of bam DEPTH for multiple regions in a terminal | 20181018 | 20210412 |
| coveragematrix | generate a VCF file from bam coverage | 20200618 | 20200618 |
| coverageplotter | Display an image of depth to display any anomaly an intervals+bams | 20200605 | 20221125 |
| findgvcfsblocks | Find common blocks of calleable regions from a set of gvcfs | 20210806 | 20220401 |
| groupbygene | Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff | 20131209 | 20220529 |
| gtf2bed | Convert GTF/GFF3 to BED. | 20220629 | 20220630 |
| lowresbam2raster | Low Resolution BAM to raster graphics | 20170523 | 20211126 |
| mkminibam | Creates an archive of small bams with only a few regions. | 20190410 | 20221019 |
| prettysam | Pretty SAM alignments | 20171215 | 20211105 |
| pubmed404 | Test if URL in the pubmed abstracts are reacheable. | 20181210 | 20200204 |
| pubmedcodinglang | Programming language use distribution from recent programs / articles | 20170404 | 20200223 |
| pubmedgender | Add gender-related attributes in the Author tag of pubmed xml. |  |  |
| pubmedgraph | Creates a Gephi-gexf graph of references-cites for a given PMID | 20150605 | 20200220 |
| sam2tsv | Prints the SAM alignments as a TAB delimited file. | 20170712 | 20210304 |
| samgrep | grep read-names in a bam file | 20130506 | 20210726 |
| samrmdupnames | remove duplicated names in sorted BAM | 20221207 | 20221207 |
| swingbamcov | Bam coverage viewer using Java Swing UI | 20210420 | 20220513 |
| swingbamview | Read viewer using Java Swing UI | 20220503 | 20230124 |
| swingindexcov | indexcov visualization | 2020511 | 2020512 |
| swingvcfjexl | Filter VCF using Java Swing UI and JEXL/Javascript expression | 20220413 | 20220414 |
| swingvcfview | VCFviewer using Java Swing UI | 20210503 | 20210503 |
| uniprot2svg | plot uniprot to SVG | 20220608 | 20220922 |
| vcf2table | convert a vcf to a table, to ease display in the terminal | 20170511 | 20220507 |
| vcfallelebalance | Insert missing allele balance annotation using FORMAT:AD | 20180829 | 20200805 |
| vcfbigbed | Annotate a VCF with values from a bigbed file | 20220107 | 20220107 |
| vcfbigwig | Annotate a VCF with values from a bigwig file | 20200506 | 20220110 |
| vcffilterjdk | Filtering VCF with dynamically-compiled java expressions | 20170705 | 20220830 |
| vcffilterso | Filter a VCF file annotated with SNPEff or VEP with terms from Sequence-Ontology. Reasoning : Children of user's SO-terms will be also used. | 20170331 | 20200924 |
| vcfgenesplitter | Split VCF+VEP by gene/transcript. | 20160310 | 202220531 |
| vcfgnomad | Peek annotations from gnomad | 20170407 | 20200702 |
| vcfhead | print the first variants of a vcf | 20131210 | 20200518 |
| vcfpar | Flag human sexual regions excluding PAR. | 20200908 | 20200908 |
| vcfpolyx | Number of repeated REF bases around POS. | 20200930 | 20211102 |
| vcfrebase | Restriction sites overlaping variations in a vcf | 20131115 | 20200624 |
| vcfrenamechr | Convert the names of the chromosomes in a VCF file |  | 20190411 |
| vcfshuffle | Shuffle a VCF | 20131210 | 20200818 |
| vcfsplitnvariants | Split VCF to 'N' VCF files | 202221122 | 202221201 |
| vcfstrech2svg | another VCF to SVG | 20210304 | 20210309 |
| vcftail | print the last variants of a vcf | 20131210 | 20200518 |
| wescnvsvg | SVG visualization of bam DEPTH for multiple regions | 20180726 | 20210726 |
| wgscoverageplotter | Whole genome coverage plotter | 20201125 | 20210812 |
