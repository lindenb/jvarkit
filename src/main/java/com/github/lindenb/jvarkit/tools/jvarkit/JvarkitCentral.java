/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.jvarkit;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.codec.language.Soundex;

import com.github.lindenb.jvarkit.tools.allelebalance.VcfAlleleBalance;
import com.github.lindenb.jvarkit.tools.backlocate.BackLocate;
import com.github.lindenb.jvarkit.tools.bam2graphics.Bam2Raster;
import com.github.lindenb.jvarkit.tools.bam2graphics.CoveragePlotter;
import com.github.lindenb.jvarkit.tools.bam2graphics.LowResBam2Raster;
import com.github.lindenb.jvarkit.tools.bam2graphics.WGSCoveragePlotter;
import com.github.lindenb.jvarkit.tools.bam2svg.BamToSVG;
import com.github.lindenb.jvarkit.tools.bam2svg.SvToSVG;
import com.github.lindenb.jvarkit.tools.bam2svg.WesCnvSvg;
import com.github.lindenb.jvarkit.tools.bam2wig.Bam2Wig;
import com.github.lindenb.jvarkit.tools.bam2xml.Bam2Xml;
import com.github.lindenb.jvarkit.tools.bamleftalign.BamLeftAlign;
import com.github.lindenb.jvarkit.tools.bamstats04.BamStats04;
import com.github.lindenb.jvarkit.tools.bamstats04.BamStats05;
import com.github.lindenb.jvarkit.tools.barcode.BarcodeGenerator;
import com.github.lindenb.jvarkit.tools.basecoverage.BaseCoverage;
import com.github.lindenb.jvarkit.tools.basecoverage.CNVPaneOfNormal;
import com.github.lindenb.jvarkit.tools.batchpicts.BatchIGVPictures;
import com.github.lindenb.jvarkit.tools.bcftools.PlotBcftoolsStats;
import com.github.lindenb.jvarkit.tools.bcftoolsmergebest.BCFToolsMergeBest;
import com.github.lindenb.jvarkit.tools.bed2vcf.BedToVcf;
import com.github.lindenb.jvarkit.tools.bedcluster.BedCluster;
import com.github.lindenb.jvarkit.tools.bedclustername.BedClusterName;
import com.github.lindenb.jvarkit.tools.bedrenamechr.BedRenameChromosomes;
import com.github.lindenb.jvarkit.tools.bedtools.BedMergeCnv;
import com.github.lindenb.jvarkit.tools.bedtools.BedNonOverlappingSet;
import com.github.lindenb.jvarkit.tools.bedtools.BedRemoveBed;
import com.github.lindenb.jvarkit.tools.bgen.bgen2vcf.BGenToVcf;
import com.github.lindenb.jvarkit.tools.bgen.bgenview.BGenView;
import com.github.lindenb.jvarkit.tools.bigwigmerge.BigwigMerge;
import com.github.lindenb.jvarkit.tools.bigwigtview.BigWigTView;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidae;
//import com.github.lindenb.jvarkit.tools.bio2rdf.BioToRDF;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk;
import com.github.lindenb.jvarkit.tools.biostar.Biostar103303;
import com.github.lindenb.jvarkit.tools.biostar.Biostar105754;
import com.github.lindenb.jvarkit.tools.biostar.Biostar130456;
import com.github.lindenb.jvarkit.tools.biostar.Biostar139647;
import com.github.lindenb.jvarkit.tools.biostar.Biostar145820;
import com.github.lindenb.jvarkit.tools.biostar.Biostar154220;
import com.github.lindenb.jvarkit.tools.biostar.Biostar160470;
import com.github.lindenb.jvarkit.tools.biostar.Biostar165777;
import com.github.lindenb.jvarkit.tools.biostar.Biostar170742;
import com.github.lindenb.jvarkit.tools.biostar.Biostar172515;
import com.github.lindenb.jvarkit.tools.biostar.Biostar173114;
import com.github.lindenb.jvarkit.tools.biostar.Biostar175929;
import com.github.lindenb.jvarkit.tools.biostar.Biostar178713;
import com.github.lindenb.jvarkit.tools.biostar.Biostar214299;
import com.github.lindenb.jvarkit.tools.biostar.Biostar234081;
import com.github.lindenb.jvarkit.tools.biostar.Biostar234230;
import com.github.lindenb.jvarkit.tools.biostar.Biostar251649;
import com.github.lindenb.jvarkit.tools.biostar.Biostar322664;
import com.github.lindenb.jvarkit.tools.biostar.Biostar332826;
import com.github.lindenb.jvarkit.tools.biostar.Biostar336589;
import com.github.lindenb.jvarkit.tools.biostar.Biostar352930;
import com.github.lindenb.jvarkit.tools.biostar.Biostar3654;
import com.github.lindenb.jvarkit.tools.biostar.Biostar398854;
import com.github.lindenb.jvarkit.tools.biostar.Biostar404363;
import com.github.lindenb.jvarkit.tools.biostar.Biostar480685;
import com.github.lindenb.jvarkit.tools.biostar.Biostar489074;
import com.github.lindenb.jvarkit.tools.biostar.Biostar497922;
import com.github.lindenb.jvarkit.tools.biostar.Biostar59647;
import com.github.lindenb.jvarkit.tools.biostar.Biostar76892;
import com.github.lindenb.jvarkit.tools.biostar.Biostar77288;
import com.github.lindenb.jvarkit.tools.biostar.Biostar77828;
import com.github.lindenb.jvarkit.tools.biostar.Biostar78285;
import com.github.lindenb.jvarkit.tools.biostar.Biostar81455;
import com.github.lindenb.jvarkit.tools.biostar.Biostar84452;
import com.github.lindenb.jvarkit.tools.biostar.Biostar84786;
import com.github.lindenb.jvarkit.tools.biostar.Biostar86363;
import com.github.lindenb.jvarkit.tools.biostar.Biostar86480;
import com.github.lindenb.jvarkit.tools.biostar.Biostar90204;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9462889;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9469733;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9501110;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9556602;
import com.github.lindenb.jvarkit.tools.biostar.Biostar95652;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9566948;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9608448;
import com.github.lindenb.jvarkit.tools.blast.BlastFilterJS;
import com.github.lindenb.jvarkit.tools.blast.BlastNToSnp;
import com.github.lindenb.jvarkit.tools.blast.MergeBlastXml;
import com.github.lindenb.jvarkit.tools.blast.MergeSplittedBlast;
import com.github.lindenb.jvarkit.tools.blast.ReduceBlast;
import com.github.lindenb.jvarkit.tools.blast2sam.BlastToSam;
import com.github.lindenb.jvarkit.tools.blastmapannots.BlastMapAnnotations;
import com.github.lindenb.jvarkit.tools.braiding.VcfBraiding;
import com.github.lindenb.jvarkit.tools.burden.OptimizeFisher;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherH;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenMAF;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenSlidingWindow;
import com.github.lindenb.jvarkit.tools.burden.VcfMoveFiltersToInfo;
import com.github.lindenb.jvarkit.tools.burdencnv.VcfBurdenCNV;
import com.github.lindenb.jvarkit.tools.cadd.VcfCadd;
import com.github.lindenb.jvarkit.tools.calling.MiniCaller;
import com.github.lindenb.jvarkit.tools.cmpbams.BamMatrix;
import com.github.lindenb.jvarkit.tools.cmpbams.CommBams;
import com.github.lindenb.jvarkit.tools.cmpbams.CompareBamAndBuild;
import com.github.lindenb.jvarkit.tools.cmpbams.CompareBams;
import com.github.lindenb.jvarkit.tools.cmpbams.CompareBams4;
import com.github.lindenb.jvarkit.tools.coveragegrid.CoverageGrid;
import com.github.lindenb.jvarkit.tools.coverageserver.CoverageServer;
import com.github.lindenb.jvarkit.tools.dbsnp.BuildDbsnp;
import com.github.lindenb.jvarkit.tools.dict2R.DictToR;
import com.github.lindenb.jvarkit.tools.dict2bed.DictToBed;
import com.github.lindenb.jvarkit.tools.dict2xml.DictToXml;
import com.github.lindenb.jvarkit.tools.drageninv.DragenBndToInversion;
import com.github.lindenb.jvarkit.tools.eva.EVADumpFiles;
import com.github.lindenb.jvarkit.tools.extendref.ExtendReferenceWithReads;
import com.github.lindenb.jvarkit.tools.fastq.BamToFastq;
import com.github.lindenb.jvarkit.tools.fastq.FastqShuffle;
import com.github.lindenb.jvarkit.tools.fastq.RepairFastq;
import com.github.lindenb.jvarkit.tools.findallcov.FindAllCoverageAtPosition;
import com.github.lindenb.jvarkit.tools.findhtsfiledict.FindHtsFileDictionary;
import com.github.lindenb.jvarkit.tools.fixvcfmissinggenotypes.FixVcfMissingGenotypes;
import com.github.lindenb.jvarkit.tools.gatk.GatkHaplotypeCaller;
import com.github.lindenb.jvarkit.tools.gff2fa.Gff3ToFasta;
import com.github.lindenb.jvarkit.tools.gff2kg.Gff2KnownGene;
import com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad;
import com.github.lindenb.jvarkit.tools.gnomad.VcfGnomadSV;
import com.github.lindenb.jvarkit.tools.go.GoUtils;
import com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene;
import com.github.lindenb.jvarkit.tools.gtex.GtexRsToQTL;
import com.github.lindenb.jvarkit.tools.gtf.Gtf2Xml;
import com.github.lindenb.jvarkit.tools.gtf.GtfLiftOver;
import com.github.lindenb.jvarkit.tools.gtf.GtfToBed;
import com.github.lindenb.jvarkit.tools.gtf.GtfToGff;
import com.github.lindenb.jvarkit.tools.gvcf.FindGVCFsBlocks;
import com.github.lindenb.jvarkit.tools.haplogroupcasectrl.HaploGroupCaseControl;
import com.github.lindenb.jvarkit.tools.hilbert.VcfToHilbert;
import com.github.lindenb.jvarkit.tools.htsvelocity.HtsVelocity;
import com.github.lindenb.jvarkit.tools.ibddb.IbdToVcf;
import com.github.lindenb.jvarkit.tools.jbrowse2.JBrowse2Server;
import com.github.lindenb.jvarkit.tools.kg2bed.KnownGenesToBed;
import com.github.lindenb.jvarkit.tools.kg2fa.KnownGeneToFasta;
import com.github.lindenb.jvarkit.tools.kg2gff.KgToGff;
import com.github.lindenb.jvarkit.tools.labguru.ScanLabGuru;
import com.github.lindenb.jvarkit.tools.liftover.BamLiftOver;
import com.github.lindenb.jvarkit.tools.liftover.BedLiftOver;
import com.github.lindenb.jvarkit.tools.liftover.ConvertLiftOverChain;
import com.github.lindenb.jvarkit.tools.liftover.LiftOverToSVG;
import com.github.lindenb.jvarkit.tools.liftover.VcfFilterByLiftOver;
import com.github.lindenb.jvarkit.tools.liftover.VcfLiftOver;
import com.github.lindenb.jvarkit.tools.manhattan.Manhattan;
import com.github.lindenb.jvarkit.tools.minibam.MakeMiniBam;
import com.github.lindenb.jvarkit.tools.misc.AddLinearIndexToBed;
import com.github.lindenb.jvarkit.tools.misc.AlleleFrequencyCalculator;
import com.github.lindenb.jvarkit.tools.misc.BamClipToInsertion;
import com.github.lindenb.jvarkit.tools.misc.BamCmpCoverage;
import com.github.lindenb.jvarkit.tools.misc.BamTile;
import com.github.lindenb.jvarkit.tools.misc.BamToSql;
import com.github.lindenb.jvarkit.tools.misc.BedIndexTabix;
import com.github.lindenb.jvarkit.tools.misc.ConvertBamChromosomes;
import com.github.lindenb.jvarkit.tools.misc.ConvertVcfChromosomes;
import com.github.lindenb.jvarkit.tools.misc.CytobandToSvg;
import com.github.lindenb.jvarkit.tools.misc.FastqGrep;
import com.github.lindenb.jvarkit.tools.misc.FastqRevComp;
import com.github.lindenb.jvarkit.tools.misc.FindAVariation;
import com.github.lindenb.jvarkit.tools.misc.HowManyBamDict;
import com.github.lindenb.jvarkit.tools.misc.IlluminaDirectory;
import com.github.lindenb.jvarkit.tools.misc.SamAddPI;
import com.github.lindenb.jvarkit.tools.misc.SamClipIndelFraction;
import com.github.lindenb.jvarkit.tools.misc.SamToJson;
import com.github.lindenb.jvarkit.tools.misc.SamToPsl;
import com.github.lindenb.jvarkit.tools.misc.SortSamRefName;
import com.github.lindenb.jvarkit.tools.misc.SplitVcf;
import com.github.lindenb.jvarkit.tools.misc.VCFShuffle;
import com.github.lindenb.jvarkit.tools.misc.VariantsInWindow;
import com.github.lindenb.jvarkit.tools.misc.VcfDistanceBetweenVariants;
import com.github.lindenb.jvarkit.tools.misc.VcfSetSequenceDictionary;
import com.github.lindenb.jvarkit.tools.misc.XsltStream;
import com.github.lindenb.jvarkit.tools.mosdepth.PlotMosdepth;
import com.github.lindenb.jvarkit.tools.msa2vcf.MsaToVcf;
import com.github.lindenb.jvarkit.tools.multiqc.MultiqcPostProcessor;
import com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesSummary;
import com.github.lindenb.jvarkit.tools.nobai.BamWithoutBai;
import com.github.lindenb.jvarkit.tools.obo.OboUtils;
import com.github.lindenb.jvarkit.tools.onekgenomes.VcfAncestralAllele;
import com.github.lindenb.jvarkit.tools.onesamplevcf.VcfMultiToOne;
import com.github.lindenb.jvarkit.tools.pcr.BamSliceBed;
import com.github.lindenb.jvarkit.tools.pcr.PcrClipReads;
import com.github.lindenb.jvarkit.tools.pcr.PcrSliceReads;
import com.github.lindenb.jvarkit.tools.phased.BamPhased01;
import com.github.lindenb.jvarkit.tools.phased.BamToHaplotypes;
import com.github.lindenb.jvarkit.tools.phased.BamToMNV;
import com.github.lindenb.jvarkit.tools.phased.VcfPhased01;
import com.github.lindenb.jvarkit.tools.plink.SwingPLinkSelectCluster;
import com.github.lindenb.jvarkit.tools.prs.VcfSamplesPRS;
import com.github.lindenb.jvarkit.tools.pubmed.Pubmed404;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedAuthorGraph;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedCodingLanguages;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedDump;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedGender;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedGraph;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedMap;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedOrcidGraph;
import com.github.lindenb.jvarkit.tools.qqplot.QQPlotter;
import com.github.lindenb.jvarkit.tools.redon.CopyNumber01;
//import com.github.lindenb.jvarkit.tools.rdfcombine.RDFCombine;
import com.github.lindenb.jvarkit.tools.ref2html.ReferenceToHtml;
import com.github.lindenb.jvarkit.tools.ref2vcf.ReferenceToVCF;
import com.github.lindenb.jvarkit.tools.regenie.RegenieBedAnnot;
import com.github.lindenb.jvarkit.tools.regenie.RegenieFunctionalAnnot;
import com.github.lindenb.jvarkit.tools.regenie.RegenieMakeAnnot;
import com.github.lindenb.jvarkit.tools.regenie.RegenieSlidingAnnot;
import com.github.lindenb.jvarkit.tools.regenie.RegenieSwing;
import com.github.lindenb.jvarkit.tools.retrocopy.GtfRetroCopy;
import com.github.lindenb.jvarkit.tools.retrocopy.KnownRetroCopy;
import com.github.lindenb.jvarkit.tools.retrocopy.ScanRetroCopy;
import com.github.lindenb.jvarkit.tools.retrocopy.StarRetroCopy;
import com.github.lindenb.jvarkit.tools.rnaseqpolya.RNASeqPolyA;
import com.github.lindenb.jvarkit.tools.sam2tsv.CnvTView;
import com.github.lindenb.jvarkit.tools.sam2tsv.PrettySam;
import com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv;
import com.github.lindenb.jvarkit.tools.sam4weblogo.SAM4WebLogo;
import com.github.lindenb.jvarkit.tools.samedict.SameDict;
import com.github.lindenb.jvarkit.tools.samfixcigar.SamFixCigar;
import com.github.lindenb.jvarkit.tools.samgrep.SamGrep;
import com.github.lindenb.jvarkit.tools.samjs.SamJdk;
import com.github.lindenb.jvarkit.tools.samrmdupnames.SamRemoveDuplicatedNames;
import com.github.lindenb.jvarkit.tools.samslop.SamSlop;
import com.github.lindenb.jvarkit.tools.sashimi.PlotSashimi;
import com.github.lindenb.jvarkit.tools.setfile.SetFileCluster;
import com.github.lindenb.jvarkit.tools.setfile.SetFileFromBed;
import com.github.lindenb.jvarkit.tools.setfile.SetFileToBed;
import com.github.lindenb.jvarkit.tools.setfile.SetFileTools;
import com.github.lindenb.jvarkit.tools.shiftbam.ShiftBam;
import com.github.lindenb.jvarkit.tools.shiftvcf.ShiftVcf;
import com.github.lindenb.jvarkit.tools.sortvcfonref.AlmostSortedVcf;
import com.github.lindenb.jvarkit.tools.sortvcfonref.SortVcfOnInfo;
import com.github.lindenb.jvarkit.tools.spliceai.VcfSpliceAI;
import com.github.lindenb.jvarkit.tools.structvar.CoverageMatrix;
import com.github.lindenb.jvarkit.tools.structvar.KnownDeletion;
import com.github.lindenb.jvarkit.tools.structvar.SVCasesControls;
import com.github.lindenb.jvarkit.tools.structvar.SamExtractClip;
import com.github.lindenb.jvarkit.tools.structvar.SamFindClippedRegions;
import com.github.lindenb.jvarkit.tools.structvar.ScanStructuralVariants;
import com.github.lindenb.jvarkit.tools.structvar.VcfStrechToSvg;
import com.github.lindenb.jvarkit.tools.structvar.breakdancer.BreakdancerToVcf;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.IndexCovToVcf;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.SwingIndexCov;
import com.github.lindenb.jvarkit.tools.structvar.manta.MantaMerger;
import com.github.lindenb.jvarkit.tools.sv2fasta.StructuralVariantToFasta;
import com.github.lindenb.jvarkit.tools.taxonomy.NcbiTaxonomyToXml;
import com.github.lindenb.jvarkit.tools.tbi2bed.VcfTbiToBed;
import com.github.lindenb.jvarkit.tools.testng.RunTestNG;
import com.github.lindenb.jvarkit.tools.textbam.TextBam;
import com.github.lindenb.jvarkit.tools.translategff3.TranslateGff3;
import com.github.lindenb.jvarkit.tools.tss.TSSEnrichment;
import com.github.lindenb.jvarkit.tools.tview.TViewCmd;
import com.github.lindenb.jvarkit.tools.ukbiobank.UKBiobankSelectSamples;
import com.github.lindenb.jvarkit.tools.ukbiobank.UkbiobankDump;
import com.github.lindenb.jvarkit.tools.uniprot.MapUniProtFeatures;
import com.github.lindenb.jvarkit.tools.uniprot.UniprotFilterJS;
import com.github.lindenb.jvarkit.tools.uniprot.UniprotToSvg;
import com.github.lindenb.jvarkit.tools.upstreamorf.Gff3UpstreamOrf;
import com.github.lindenb.jvarkit.tools.upstreamorf.VcfScanUpstreamOrf;
import com.github.lindenb.jvarkit.tools.validatorserver.CNVValidatorServer;
import com.github.lindenb.jvarkit.tools.vcf2bam.VcfToBam;
import com.github.lindenb.jvarkit.tools.vcf2intervals.VcfToIntervals;
import com.github.lindenb.jvarkit.tools.vcf2r.VcfToRScript;
import com.github.lindenb.jvarkit.tools.vcf2rdf.VcfToRdf;
import com.github.lindenb.jvarkit.tools.vcf2table.VcfToTable;
import com.github.lindenb.jvarkit.tools.vcf2xml.Vcf2Xml;
import com.github.lindenb.jvarkit.tools.vcfannot.VCFSVAnnotator;
import com.github.lindenb.jvarkit.tools.vcfannot.VCFCombineTwoSnvs;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VcfBigBed;
import com.github.lindenb.jvarkit.tools.vcfcomposite.VCFComposite;
import com.github.lindenb.jvarkit.tools.vcfconcat.VcfConcat;
import com.github.lindenb.jvarkit.tools.vcffiltergenes.VcFilterGenes;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcfflatten.VCFFlatten;
import com.github.lindenb.jvarkit.tools.vcfgatkeval.VcfGatkEval;
import com.github.lindenb.jvarkit.tools.vcfgrantham.VcfGrantham;
import com.github.lindenb.jvarkit.tools.vcfgroupbypop.VcfGroupByPopulation;
import com.github.lindenb.jvarkit.tools.vcfgtf.VcfFilterGtf;
import com.github.lindenb.jvarkit.tools.vcfhead.VcfHead;
import com.github.lindenb.jvarkit.tools.vcfmulti2oneallele.VcfMultiToOneAllele;
import com.github.lindenb.jvarkit.tools.vcfmulti2oneinfo.VcfMultiToOneInfo;
import com.github.lindenb.jvarkit.tools.vcfnearest.VCFNearest;
import com.github.lindenb.jvarkit.tools.vcfpar.VcfPseudoAutosomalRegion;
import com.github.lindenb.jvarkit.tools.vcfpolyx.VCFPolyX;
import com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase;
import com.github.lindenb.jvarkit.tools.vcfregulomedb.VcfRegulomeDB;
import com.github.lindenb.jvarkit.tools.vcfserver.VcfServer;
import com.github.lindenb.jvarkit.tools.vcfsplit.VcfSplitNVariants;
import com.github.lindenb.jvarkit.tools.vcfsplitgene.VcfGeneSplitter;
import com.github.lindenb.jvarkit.tools.vcfsplitvep.VCFSplitVEP;
import com.github.lindenb.jvarkit.tools.vcfspring.VcfSpringFilter;
import com.github.lindenb.jvarkit.tools.vcfstats.VcfStats;
import com.github.lindenb.jvarkit.tools.vcfstats.VcfStats2;
import com.github.lindenb.jvarkit.tools.vcftabixml.VCFTabixml;
import com.github.lindenb.jvarkit.tools.vcftail.VcfTail;
import com.github.lindenb.jvarkit.tools.vcftrios.VCFTrios;
import com.github.lindenb.jvarkit.tools.vcfukbb.VcfUkbiobank;
import com.github.lindenb.jvarkit.tools.vcfvcf.VcfPeekAf;
import com.github.lindenb.jvarkit.tools.vcfvcf.VcfPeekVcf;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamCov;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamView;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingVcfJexlFilter;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingVcfView;
import com.github.lindenb.jvarkit.tools.velocity.ApplyVelocity;
import com.github.lindenb.jvarkit.tools.viewmate.SamViewWithMate;
import com.github.lindenb.jvarkit.tools.wib.WibToBedGraph;
import com.github.lindenb.jvarkit.tools.xcontamination.XContaminations;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;

public class JvarkitCentral {
	private static final Logger LOG = Logger.of(JvarkitCentral.class);
	private class Command {
		private final Class<?> clazz;
		boolean hidden = false;
		Command(final Class<?> claz) {
			this.clazz = claz;
			}
		private Program program() {
				return this.clazz.getAnnotation(Program.class);
			}
		public boolean isHidden() {
			final Program p = program();
			if(p!=null && (p.jvarkit_hidden() || !p.generate_doc() || !p.deprecatedMsg().isEmpty())) return true;
			return hidden;
			}
		public String getMarkdownFile() {
			return this.clazz.getSimpleName()+".md";
			}
		public String getName() {
			final Program p = program();
			return p==null?this.clazz.getSimpleName().toLowerCase():p.name();
			}
		public String getDescription() {
			final Program p = program();
			return p==null?getName():p.description().trim();
			}
		public Set<String> getMenus() {
			final Program p = program();
			final String s= p==null?"Unclassified":p.menu();
			return Arrays.stream(s.split("[,;]")).map(S->S.trim()).filter(S->!S.isEmpty()).collect(Collectors.toSet());
			}
		
		public String getCreationDate() {
			final Program p = program();
			return p==null?".":p.creationDate();
			}
		public String getModificationDate() {
			final Program p = program();
			return p==null?".":p.modificationDate();
			}
		
		void execute(final String[] args) {
			try {
				final Method mainMethod = this.clazz.getMethod("main",String[].class);
				mainMethod.invoke(null, (Object)args);
				}
			catch(final Throwable err) {
				LOG.error(err);
				System.exit(-1);
				}
			}
		@SuppressWarnings("unused")
		Command setHidden() {
			this.hidden = true;
			return this;
			}
		
		
		
		void writeDoc(final File dir) throws IOException {
			final PrintStream old = System.out;
			final File filename = new File(dir,getMarkdownFile());
			final String[] cmd = new String[]{"--help","--helpFormat","markdown"};
			LOG.info((filename.exists()?"over":"")+"writing doc for "+getName()+" into "+filename);
			try (PrintStream out = new PrintStream(filename)) {
				System.setOut(out);
				if(Launcher.class.isAssignableFrom(this.clazz)) {
					Launcher instance= Launcher.class.cast(this.clazz.getConstructor().newInstance());
					final int  ret= instance.instanceMain(Arrays.asList(cmd));
					instance = null;
					if(ret!=0) {
						LOG.warn("return status!=0");
						}
					}
				else
					{
					execute(cmd);
					}
				}
			catch(Throwable err) {
				LOG.error("Cannot save doc for "+getName(), err);
				}
			finally
				{
				System.setOut(old);
				}
			}
		}
	private final Map<String,Command> commands = new TreeMap<>();
	private Command command(Class<?> clazz) {
		Command c = new Command(clazz);
		if(this.commands.containsKey(c.getName())) {
			return this.commands.get(c.getName());
			}
		this.commands.put(c.getName(),c);
		return c;
		}
	
	
	private void usageMain(final PrintStream out) {
		out.println("JVARKIT");
		out.println("=======");
		out.println();
		out.println("Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.");
		out.println("Version     : " + JVarkitVersion.getInstance().getGitHash());
		out.println("Compilation : " + JVarkitVersion.getInstance().getCompilationDate());
		out.println("Github      : https://github.com/lindenb/jvarkit");
		out.println("Issues      : https://github.com/lindenb/jvarkit/issues");
		out.println();
		out.println("## Usage");
		out.println();
		out.println("```\n  java -jar jvarkit.jar [options]\n```");
		out.println("or");
		out.println("```\n  java -jar jvarkit.jar <command name> (other arguments)\n```");
		out.println();
		out.println("## Options");
		out.println();
		out.println(" + --help show this screen");
		out.println(" + --help-all show all commands, including the private ones.");
		out.println(" + --version print version");
		out.println();
		}
	
	
	private void usage(PrintStream out,boolean showHidden) {
		usageMain(out);
		
		out.println("## Tools");
		out.println();
		final int maxlenname = this.commands.values().stream().
				filter(P->showHidden || !P.isHidden()).
				mapToInt(S->S.getName().length()).max().orElse(1);
		final int maxlendesc = this.commands.values().stream().
				filter(P->showHidden || !P.isHidden()).
				mapToInt(S->S.getDescription().length()).max().orElse(1);

		out.print("| ");
		out.print(String.format("%"+maxlenname+"s","Name"));
		out.print(" | ");
		out.print(String.format("%-"+maxlendesc+"s","Description"));
		out.println(" |");
		out.print("|-");
		out.print(StringUtil.repeatCharNTimes('-', maxlenname));
		out.print("-+-");
		out.print(StringUtil.repeatCharNTimes('-', maxlendesc));
		out.print("-|");
		out.println();
		for(Command c: this.commands.values()) {
			if(!showHidden && c.isHidden()) continue;
			out.print("| ");
			out.print(String.format("%"+maxlenname+"s",c.getName()));
			out.print(" | ");
			out.print(String.format("%-"+maxlendesc+"s",c.getDescription()));
			out.print(" |");
			out.println();
			}
		out.println();
	}
	
	private String hyperlink(final String url) {
	return "["+url+"]("+url+")";
	}
	
	private void menuReadTheDoc(final PrintStream out,final String menu) {
		final List<Command> L = this.commands.values().stream().
				filter(C->C.getMenus().stream().anyMatch(X->X.equalsIgnoreCase(menu))).
				filter(C->!C.isHidden()).
				collect(Collectors.toList());
		if(L.isEmpty()) return;
		out.println("### " +menu+"\n");
		
		out.println("| Tool | Description | Creation | Update |");
		out.println("| ---: | :---------- | :------: | :----: |");
		for(Command c: L) {
			out.print("| ");
			out.print("["+c.getName()+"]("+c.getMarkdownFile()+")");
			out.print(" | ");
			out.print(c.getDescription());
			out.print(" | ");
			out.print(c.getCreationDate());
			out.print(" | ");
			out.print(c.getModificationDate());
			out.print(" |");
			out.println();
			}
		
		out.println();
		}

	
	void writeReadTheDoc(final File dir) throws IOException {
		IOUtil.assertDirectoryIsWritable(dir);
		final File docs = new File(dir,"docs");
		IOUtil.assertDirectoryIsWritable(docs);
		try (PrintStream out = new PrintStream(new File(docs,"index.md"))) {
			usageMain(out);
			out.println("## Compilation Installation\n");
			out.println("Please, read [how to run and install jvarkit](JvarkitCentral.md)\n");
			out.println("## Tools\n");
			
			for(String menu:  this.commands.values().stream().flatMap(T->T.getMenus().stream()).collect(Collectors.toSet())) {
				menuReadTheDoc(out,menu);
				}
			
			out.flush();
			out.println();
			}
		try (PrintStream out = new PrintStream(new File(docs,"JvarkitCentral.md"))) {
			usageMain(out);
			
			out.append("## Compilation\n");
			out.append("\n");
			out.append("### Requirements / Dependencies\n");
			out.append("\n");
			out.append("* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )\n");
			out.append("\n");
			out.append("\n");
			out.append("### Download a pre-compiled executable jar\n");
			out.append("\n");
			out.append("A pre-compiled java executable of jvarkit.jar is available at "+hyperlink("https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH"));
			out.append("\n");
			out.append("### Download and Compile\n");
			out.append("\n");
			out.append("```bash\n");
			out.append("$ git clone \"https://github.com/lindenb/jvarkit.git\"\n");
			out.append("$ cd jvarkit\n");
			out.append("$ ./gradlew jvarkit\n");
			out.append("```\n");
			out.append("\n");
			out.append("The java jar file will be installed in the `dist` directory.\n");
			out.append("\n");
			
			out.append("## Contribute\n");
			out.append("\n");
			out.append("- Issue Tracker: "+hyperlink("http://github.com/lindenb/jvarkit/issues")+"\n");
			out.append("- Source Code: "+hyperlink("http://github.com/lindenb/jvarkit")+"\n");
			out.append("\n");
			out.append("## License\n");
			out.append("\n");
			out.append("The project is licensed under the MIT license.\n");
			out.append("\n");
			out.append("## Citing\n");
			out.append("\n");
			out.append("Should you cite **jvarkit** ? "+hyperlink("https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md")+"\n");
			out.append("\n");
			out.append("The current reference is:\n");
			out.append("\n");
			out.append(hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
			out.append("\n");
			out.append("> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.\n");
			out.append("> "+hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
			out.append("\n");

			out.flush();
			}
		try (PrintStream out = new PrintStream(new File(dir,"mkdocs.yml"))) {
				out.println("site_name: \"Jvarkit\"");
				out.println("repo_url: \"https://github.com/lindenb/jvarkit\"");
				out.println("repo_name: \"GitHub\"");
				out.println("docs_dir: docs");
				out.println("");
				out.println("theme:");
				out.println("    name: readthedocs");
				out.println("    highlightjs: true");
				out.println("    html_show_sourcelink: false");
				out.println("    prev_next_buttons_location: both");
				out.println("");
				out.println("nav:");
				out.println("        - \"jvarkit\" : JvarkitCentral.md");
				for(Command c:this.commands.values()) {
					out.println("        - \""+ c.getName() +"\" : "+ c.getMarkdownFile());
					//out.println("        - "VcfTail": VcfTail.md");
					}
				out.println();
				out.flush();
			}
		for(Command c:this.commands.values()) {
			try {
				c.writeDoc(docs);
				}
			catch(IOException err) {
				LOG.warn(err);
				}
			}
		}
	
	private void run(final String[] args) {
		command(AlleleFrequencyCalculator.class);
		command(AlmostSortedVcf.class);
		command(ApplyVelocity.class);
		command(AddLinearIndexToBed.class);
		command(TSSEnrichment.class);
		command(BedMergeCnv.class);
		command(BarcodeGenerator.class);
		command(BackLocate.class);
		command(BamCmpCoverage.class);
		command(BamToSVG.class);
		command(BamToFastq.class);
		command(BamToHaplotypes.class);
		command(BamClipToInsertion.class);
		command(BamWithoutBai.class);
		command(BamToSql.class);
		command(Bam2Xml.class);
		command(BamPhased01.class);
		command(BamLiftOver.class);
		command(BamSliceBed.class);
		command(BamTile.class);
		command(BamToMNV.class);
		command(BamMatrix.class);
		command(Bam2Wig.class);
		command(Bam2Raster.class);
		command(BamStats04.class);
		command(BamStats05.class);
		command(BatchIGVPictures.class);
		command(SortSamRefName.class);
		command(BaseCoverage.class);
		command(BamLeftAlign.class);
		command(BlastToSam.class);
		command(BlastNToSnp.class);
		command(BedCluster.class);
		command(BedClusterName.class);
		command(BedIndexTabix.class);
		command(BedNonOverlappingSet.class);
		command(BedRemoveBed.class);
		command(BedLiftOver.class);
		command(BedToVcf.class);
		command(BGenToVcf.class);
		command(BGenView.class);
		command(BioAlcidae.class);
		command(BioAlcidaeJdk.class);
		command(BigwigMerge.class);
		command(BigWigTView.class);
		command(BCFToolsMergeBest.class);
		//command(BioToRDF.class);
		command(Biostar103303.class);
		command(Biostar105754.class);
		command(Biostar154220.class);
		command(Biostar130456.class);
		command(Biostar139647.class);
		command(Biostar145820.class);
		command(Biostar160470.class);
		command(Biostar165777.class);
		command(Biostar170742.class);
		command(Biostar172515.class);
		command(Biostar173114.class);
		command(Biostar175929.class);
		command(Biostar178713.class);
		command(Biostar214299.class);
		command(Biostar234081.class);
		command(Biostar234230.class);
		command(Biostar251649.class);
		command(Biostar322664.class);
		command(Biostar332826.class);
		command(Biostar336589.class);
		command(Biostar352930.class);
		command(Biostar3654.class);
		command(Biostar398854.class);
		command(Biostar404363.class);
		command(Biostar480685.class);
		command(Biostar489074.class);
		command(Biostar497922.class);
		command(Biostar59647.class);
		command(Biostar76892.class);
		command(Biostar77288.class);
		command(Biostar77828.class);
		command(Biostar78285.class);
		command(Biostar81455.class);
		command(Biostar84452.class);
		command(Biostar84786.class);
		command(Biostar86363.class);
		command(Biostar86480.class);
		command(Biostar90204.class);
		command(Biostar95652.class);
		command(Biostar9462889.class);
		command(Biostar9469733.class);
		command(Biostar9501110.class);
		command(Biostar9556602.class);
		command(Biostar9566948.class);
		command(Biostar9608448.class);
		command(BlastMapAnnotations.class);
		command(BlastFilterJS.class);
		command(BreakdancerToVcf.class);
		command(DictToBed.class);
		command(DictToXml.class);
		command(DictToR.class);
		command(DragenBndToInversion.class);
		command(CompareBamAndBuild.class);
		command(ConvertLiftOverChain.class);
		command(CommBams.class);
		command(CompareBams.class);
		command(CompareBams4.class);
		command(CNVPaneOfNormal.class);
		command(CoverageMatrix.class);
		command(CytobandToSvg.class);
		command(BuildDbsnp.class);
		command(CnvTView.class);
		command(ConvertBamChromosomes.class);
		command(CopyNumber01.class);
		command(BedRenameChromosomes.class);
		command(ConvertVcfChromosomes.class);
		command(CoveragePlotter.class);
		command(CoverageServer.class);
		command(CoverageGrid.class);
		command(EVADumpFiles.class);
		command(ExtendReferenceWithReads.class);
		command(FastqShuffle.class);
		command(FastqRevComp.class);
		command(FastqGrep.class);
		command(FixVcfMissingGenotypes.class);
		command(FindGVCFsBlocks.class);
		command(FindAVariation.class);
		command(FindAllCoverageAtPosition.class);
		command(FindHtsFileDictionary.class);
		command(GatkHaplotypeCaller.class);
		command(GtfToBed.class);
		command(GtfToGff.class);
		command(GtfLiftOver.class);
		command(GtexRsToQTL.class);
		command(GtfRetroCopy.class);
		command(Gff3UpstreamOrf.class);
		command(Gff3ToFasta.class);
		command(GoUtils.class);
		command(GroupByGene.class);
		command(Gtf2Xml.class);
		command(Gff2KnownGene.class);
		command(HaploGroupCaseControl.class);
		command(HowManyBamDict.class);
		command(HtsVelocity.class);
		command(IbdToVcf.class);
		command(IlluminaDirectory.class);
		command(IndexCovToVcf.class);
		command(JBrowse2Server.class);
		command(KgToGff.class);
		command(KnownRetroCopy.class);
		command(KnownGenesToBed.class);
		command(KnownGeneToFasta.class);
		command(KnownDeletion.class);
		command(LiftOverToSVG.class);
		command(LowResBam2Raster.class);
		command(MantaMerger.class);
		command(MiniCaller.class);
		command(MsaToVcf.class);
		command(Manhattan.class);
		command(MultiqcPostProcessor.class);
		command(MakeMiniBam.class);
		command(MapUniProtFeatures.class);
		command(MergeBlastXml.class);
		command(MergeSplittedBlast.class);
		command(NcbiTaxonomyToXml.class);
		command(NgsFilesSummary.class);
		command(OptimizeFisher.class);
		command(OboUtils.class);
		command(PlotBcftoolsStats.class);
		command(PlotMosdepth.class);
		command(PubmedDump.class);
		command(Pubmed404.class);
		command(PubmedCodingLanguages.class);
		command(PubmedGender.class);
		command(PubmedGraph.class);
		command(PubmedAuthorGraph.class);
		command(PubmedOrcidGraph.class);
		command(PubmedMap.class);
		command(PlotSashimi.class);
		command(PrettySam.class);
		command(PcrClipReads.class);
		command(PcrSliceReads.class);
		command(QQPlotter.class);
		command(ReferenceToHtml.class);
		command(ReduceBlast.class);
		command(RegenieMakeAnnot.class);
		command(RegenieBedAnnot.class);
		command(RegenieFunctionalAnnot.class);
		command(RegenieSlidingAnnot.class);
		command(RegenieSwing.class);
		//command(RDFCombine.class);
		command(RepairFastq.class);
		command(ReferenceToVCF.class);
		command(RNASeqPolyA.class);
		command(RunTestNG.class);
		command(VCFBigWig.class);
		command(VcfBigBed.class);
		command(VcfFilterGtf.class);
		command(VcfFilterByLiftOver.class);
		command(VcfPeekVcf.class);
		command(VCFPolyX.class);
		command(VcfSpliceAI.class);
		command(VcfToTable.class);
		command(VcfToRScript.class);
		command(VcfStats.class);
		command(VcfStats2.class);
		command(VcfGnomad.class);
		command(VcfSetSequenceDictionary.class);
		command(VcfFilterSequenceOntology.class);
		command(VcfToBam.class);
		command(Sam2Tsv.class);
		command(SamGrep.class);
		command(SamViewWithMate.class);
		command(SamJdk.class);
		command(SameDict.class);
		command(SamClipIndelFraction.class);
		command(SamToJson.class);
		command(SamSlop.class);
		command(SamRemoveDuplicatedNames.class);
		command(ScanRetroCopy.class);
		command(ScanLabGuru.class);
		command(SortVcfOnInfo.class);
		command(StarRetroCopy.class);
		command(SetFileTools.class);
		command(SetFileFromBed.class);
		command(SetFileToBed.class);
		command(SetFileCluster.class);
		command(ShiftBam.class);
		command(ShiftVcf.class);
		command(SplitVcf.class);
		command(StructuralVariantToFasta.class);
		command(SVCasesControls.class);
		command(SAM4WebLogo.class);
		command(SamAddPI.class);
		command(SamToPsl.class);
		command(SamFindClippedRegions.class);
		command(SamFixCigar.class);
		command(ScanStructuralVariants.class);
		command(SwingVcfView.class);
		command(SwingBamCov.class);
		command(SwingBamView.class);
		command(SamExtractClip.class);
		command(SwingIndexCov.class);
		command(SwingVcfJexlFilter.class);
		command(SwingPLinkSelectCluster.class);
		command(SvToSVG.class);
		command(VCFSVAnnotator.class);
		command(TextBam.class);
		command(TViewCmd.class);
		command(TranslateGff3.class);
		command(UniprotToSvg.class);
		command(UniprotFilterJS.class);
		command(UKBiobankSelectSamples.class);
		command(UkbiobankDump.class);
		command(VariantsInWindow.class);
		command(VcfAncestralAllele.class);
		command(VcfMultiToOne.class);
		command(VcfBurdenCNV.class);
		command(VcFilterGenes.class);
		command(VcfBurdenSlidingWindow.class);
		command(Vcf2Xml.class);
		command(VcfPhased01.class);
		command(VcfConcat.class);
		command(VcfCadd.class);
		command(VCFComposite.class);
		command(VcfDistanceBetweenVariants.class);
		command(VcfFilterJdk.class);
		command(VcfAlleleBalance.class);
		command(VCFFlatten.class);
		command(VcfPseudoAutosomalRegion.class);
		command(VcfRebase.class);
		command(VcfGatkEval.class);
		command(VcfSpringFilter.class);
		command(VcfGnomadSV.class);
		command(VcfGrantham.class);
		command(VcfHead.class);
		command(VcfTail.class);
		command(VcfServer.class);
		command(VcfRegulomeDB.class);
		command(VcfScanUpstreamOrf.class);
		command(VCFTrios.class);
		command(VcfBurdenMAF.class);
		command(VcfLiftOver.class);
		command(VcfPeekAf.class);
		command(VcfGeneSplitter.class);
		command(VCFSplitVEP.class);
		command(VcfMultiToOneInfo.class);
		command(VcfMultiToOneAllele.class);
		command(VcfMoveFiltersToInfo.class);
		command(VcfTbiToBed.class);
		command(VcfSplitNVariants.class);
		command(VcfStrechToSvg.class);
		command(VCFCombineTwoSnvs.class);
		command(VcfToIntervals.class);
		command(VCFShuffle.class);
		command(VcfTbiToBed.class);
		command(VcfToRdf.class);
		command(VcfToHilbert.class);
		command(VcfGroupByPopulation.class);
		command(VcfBraiding.class);
		command(VcfSamplesPRS.class);
		command(VcfUkbiobank.class);
		command(VCFNearest.class);
		command(VCFTabixml.class);
		command(CNVValidatorServer.class);
		command(VcfBurdenFisherH.class);
		command(WibToBedGraph.class);
		command(WesCnvSvg.class);
		command(WGSCoveragePlotter.class);
		command(XsltStream.class);
		command(XContaminations.class);
		if(args.length==0) {
			usage(System.err,false);
			return;
		}
		else if(args[0].equals("--version") || args[0].equals("-v")) {
			System.out.println(JVarkitVersion.getInstance().getGitHash());
			return;
		}
		else if(args[0].equals("--help") || args[0].equals("-h")) {
			usage(System.out, false);
			return;
			}
		else if(args[0].equals("--help-all")) {
			usage(System.out, true);
			return;
			}
		else if(args.length==2 && args[0].equals("--generate-doc")) {
			final File dir = new File(args[1]);
			IOUtil.assertDirectoryIsWritable(dir);
			for(Command c:this.commands.values()) {
				try {
					c.writeDoc(dir);
					}
				catch(IOException err) {
					LOG.warn(err);
					}
				}
			return;
			}
		else if(args.length==2 && args[0].equals("--readthedocs")) {
			final File dir = new File(args[1]);
			try {
				writeReadTheDoc(dir);
				}
			catch(IOException err) {
				LOG.error(err);
				System.exit(-1);
				}
			return;
			}
		else if(args[0].startsWith("-")) {
			LOG.error("exepected a sub-jvarkit-program. but got "+args[0]+". type java -jar jvarkit.jar --help for more information.");
			System.exit(-1);
			}
		final String userCmd = args[0];
		final Command cmd = this.commands.get(userCmd);
		if(cmd!=null) {
			final String[] args2=new String[args.length-1];
			System.arraycopy(args, 1, args2, 0, args2.length);
			cmd.execute(args2);
			}
		else
			{
			final Soundex soundex = new Soundex();
			final Set<String> sounds = this.commands.keySet().
				stream().
				filter(S->{try { return soundex.difference(S, userCmd)>=4/* 4: high score */;} catch(Exception err) {return false;}} ).
				collect(Collectors.toSet());
			LOG.error("jvarkit command \""+userCmd+"\" not found.\n"
					+ "Available commands are: " + String.join(", ",this.commands.keySet())+"\n"+
					(sounds.isEmpty()?"":"\nDo you mean "+String.join(", ", sounds)+" ?\n")
					);
			System.exit(-1);
			}
		}
	
    public static void main(String[] args) {
        new JvarkitCentral().run(args);
    }
}
