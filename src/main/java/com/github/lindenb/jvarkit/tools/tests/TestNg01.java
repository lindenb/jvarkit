/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.tests;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.Random;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.imageio.ImageIO;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.xml.sax.helpers.DefaultHandler;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.backlocate.BackLocate;
import com.github.lindenb.jvarkit.tools.bam2graphics.Bam2Raster;
import com.github.lindenb.jvarkit.tools.bam2graphics.LowResBam2Raster;
import com.github.lindenb.jvarkit.tools.bam2wig.Bam2Wig;
import com.github.lindenb.jvarkit.tools.biostar.Biostar84452;
import com.github.lindenb.jvarkit.tools.burden.CaseControlCanvas;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFilterExac;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFilterGenes;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherH;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherV;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenMAF;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenRscriptV;
import com.github.lindenb.jvarkit.tools.burden.VcfFilterNotInPedigree;
import com.github.lindenb.jvarkit.tools.burden.VcfInjectPedigree;
import com.github.lindenb.jvarkit.tools.burden.VcfLoopOverGenes;
import com.github.lindenb.jvarkit.tools.burden.VcfMoveFiltersToInfo;
import com.github.lindenb.jvarkit.tools.calling.MiniCaller;
import com.github.lindenb.jvarkit.tools.fastq.FastqShuffle;
import com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad;
import com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene;
import com.github.lindenb.jvarkit.tools.hilbert.VcfToHilbert;
import com.github.lindenb.jvarkit.tools.misc.BamTile;
import com.github.lindenb.jvarkit.tools.misc.BamToSql;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.tools.misc.FindAVariation;
import com.github.lindenb.jvarkit.tools.misc.FindAllCoverageAtPosition;
import com.github.lindenb.jvarkit.tools.misc.FixVcfMissingGenotypes;
import com.github.lindenb.jvarkit.tools.misc.Gff2KnownGene;
import com.github.lindenb.jvarkit.tools.misc.GoUtils;
import com.github.lindenb.jvarkit.tools.misc.IlluminaReadName;
import com.github.lindenb.jvarkit.tools.misc.PadEmptyFastq;
import com.github.lindenb.jvarkit.tools.misc.VcfCreateDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfMultiToOneAllele;
import com.github.lindenb.jvarkit.tools.misc.VcfNoCallToHomRef;
import com.github.lindenb.jvarkit.tools.misc.VcfRemoveUnusedAlt;
import com.github.lindenb.jvarkit.tools.misc.VcfSetSequenceDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfToSvg;
import com.github.lindenb.jvarkit.tools.misc.VcfToTable;
import com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesSummary;
import com.github.lindenb.jvarkit.tools.onesamplevcf.VcfMultiToOne;
import com.github.lindenb.jvarkit.tools.sortvcfonref.SortVcfOnInfo;
import com.github.lindenb.jvarkit.tools.trap.TrapIndexer;
import com.github.lindenb.jvarkit.tools.trap.VcfTrap;
import com.github.lindenb.jvarkit.tools.vcf2rdf.VcfToRdf;
import com.github.lindenb.jvarkit.tools.vcf2sql.VcfToSql;
import com.github.lindenb.jvarkit.tools.vcf2xml.Vcf2Xml;
import com.github.lindenb.jvarkit.tools.vcfamalgation.VcfXmlAmalgamation;
import com.github.lindenb.jvarkit.tools.vcfbed.VCFBed;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig;
import com.github.lindenb.jvarkit.tools.vcfcmp.VcfCompareCallers;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcffixindels.VCFFixIndels;
import com.github.lindenb.jvarkit.tools.vcflist.VcfList;
import com.github.lindenb.jvarkit.tools.vcflist.VcfOffsetsIndexFactory;
import com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase;
import com.github.lindenb.jvarkit.tools.vcfstats.VcfStats;
import com.github.lindenb.jvarkit.tools.vcfstripannot.VCFStripAnnotations;
import com.github.lindenb.jvarkit.tools.vcftrios.VCFFamilies;
import com.github.lindenb.jvarkit.tools.vcftrios.VCFTrios;
import com.github.lindenb.jvarkit.util.Algorithms;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequenceReader;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceContig;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenome;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenomeFactory;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

class TestNg01 {
	static final File TEST_RESULTS_DIR= new File("test-results");
	static final File JETER_VCF = new File(TEST_RESULTS_DIR,"jeter.vcf");
	static final String TOY_FA="src/test/resources/toy.fa";
	static final String TOY_BED_GZ="src/test/resources/toy.bed.gz";
	static final String TOY_VCF_GZ="src/test/resources/toy.vcf.gz";
	static final String TOY_BAM="src/test/resources/toy.bam";
	static final String TOY_DICT="src/test/resources/toy.dict";
	static final String VCF01 = "src/test/resources/test_vcf01.vcf";
	static final String PED01 = "src/test/resources/test_vcf01.ped";
	static final String KNOWN_GENES01 = "src/test/resources/test_vcf01.knownGenes.txt.gz";
	
	private Properties properties = new Properties();
	
	
	@BeforeClass
    public void setup() throws IOException {
		TEST_RESULTS_DIR.mkdirs();
		final File propFile = new File(TEST_RESULTS_DIR,"properties.xml");
		if(propFile.exists()) {
			FileInputStream in = new FileInputStream(propFile);
			this.properties.loadFromXML(in);
			in.close();
			}
		}
	
	private static void assertIsXml(final File f) {
		boolean b=false; 
		try {
			Assert.assertTrue(f.exists(), "file "+f+" should exist");
			final SAXParser sax=SAXParserFactory.newInstance().newSAXParser();
			sax.parse(f, new DefaultHandler());
			b= true;
		}catch (final Exception e) {
			b= false;
		}
		Assert.assertTrue(b,"file "+f+" should be xml");
	}
	
	private static void assertTableIsConsitent(final File f,Predicate<String> ignoreLine) {
		final Pattern tab=Pattern.compile("[\t]");
		BufferedReader r=null;
		try {
			Assert.assertTrue(f.exists(), "file "+f+" should exist");
			r= IOUtils.openFileForBufferedReading(f);
			int nCols=-1;
			int nRow=0;
			String line;
			while((line=r.readLine())!=null)
				{
				if(ignoreLine!=null && ignoreLine.test(line)) continue;
				nRow++;
				final String tokens[]=tab.split(line);
				final int c= tokens.length;
				if(nRow==1)
					{
					nCols = c;
					}
				else 
					{
					if(nCols!=c)
						{
						for(int i=0;i< tokens.length;i++)
							{
							System.err.println("$"+(i+1)+":"+tokens[i]);
							}
						}
					Assert.assertEquals(nCols,c,"Line "+line+" expected "+nCols+" tokens but got "+c);
					}
				}
			}
		catch (final Exception e) {
			e.printStackTrace();
			Assert.fail();
			}
		finally
			{
			CloserUtil.close(r);
			}
	}
	
	
	private static void assertIsVcf(final File f) {
		boolean b; 
		try {
			streamVcf(f).
				flatMap(V->V.getGenotypes().stream()).
				flatMap(G->G.getAlleles().stream());
			b= true;
		}catch (final Exception e) {
			b= false;
		}
		Assert.assertTrue(b,"file "+f+" should be vcf");
	}

	
	private static void assertIsBam(final File f) {
		boolean b = false; 
		 try {
			 final SamReader r=SamReaderFactory.makeDefault().open(f);
			 r.iterator().stream().forEach(R->R.getReadName());
			 r.close();
			b = true;
		}catch (final Exception e) {
			b = false;
		}
		Assert.assertTrue(b,"file "+f+" should be bam.");
	}

	
	
	@DataProvider(name = "all_vcfs")
	public Object[][]  all_vcfs() {
		return new Object[][] {
			{TOY_VCF_GZ},
			{VCF01},
			{"src/test/resources/ExAC.r1.sites.vep.vcf.gz"},
			{"src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz"},
			{"src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz"}
			};
		}
	
	 static Stream<VariantContext> streamVcf(final File f) {
		final VCFFileReader r = new VCFFileReader(f,false);
		final CloseableIterator<VariantContext> iter = r.iterator();
		return StreamSupport.stream(new IterableAdapter<VariantContext>(iter).spliterator(), false).onClose(()->{iter.close();r.close();});
		}
	 static Stream<VariantContext> streamJeterVcf() {
		return streamVcf(JETER_VCF);
		}
	
	 static VCFHeader getVcfHeader(final File f) {
		final VCFFileReader r = new VCFFileReader(f,false);
		final VCFHeader h = r.getFileHeader();
		r.close();
		return h;
	 	}
    
    @Test(dataProvider="all_vcfs")
    public void testVcfCreateSequenceDictionary(final String vcfPath) throws IOException{
    	File output = new File(TEST_RESULTS_DIR,"jeter.dict");
    
        Assert.assertEquals(0,new VcfCreateDictionary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		vcfPath
        	}));
        Assert.assertTrue( output.exists());
    	}
    
    @Test
    public void testConcatSam() throws IOException{
    	File output = new File(TEST_RESULTS_DIR,"jeter.sam");
    
        Assert.assertEquals(0,new ConcatSam().instanceMain(new String[]{
        		"-o",output.getPath(),
        		TOY_BAM
        	}));
        Assert.assertTrue( output.exists());
    	}
    
    @Test
    public void testGroupByGene() throws IOException{
    	final File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new GroupByGene().instanceMain(new String[]{
        		"-o",output.getPath(),
        		VCF01
        	}));
        assertTableIsConsitent(output, null);
        Assert.assertTrue( output.exists());
    	}

    
    
    
    
    @Test(dataProvider="all_vcfs")
    public void testVcfFiltrerJdkVcf(final String vcfPath) throws IOException{
        Assert.assertEquals(0,new VcfFilterJdk().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-e","return variant.getStart()%2==0;",
        		vcfPath
        	}));
        Assert.assertTrue( JETER_VCF.exists());
        Assert.assertEquals(streamVcf(JETER_VCF).filter(V->V.getStart()%2!=0).count(),0L);
    	}
    
    @Test
    public void testBamToRaster() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
    	Assert.assertEquals(0,new Bam2Raster().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-R",TOY_FA,
        		"-r","ref:1-100",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}
    @Test
    public void testBamToSql() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.sql");
    	Assert.assertEquals(0,new BamToSql().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-R",TOY_FA,
        		"-r","ref:1-100",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}
    @Test
    public void testLowResBamToRaster() throws IOException{
	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
	Assert.assertEquals(0,new LowResBam2Raster().instanceMain(new String[]{
    		"-o",tmp.getPath(),
    		"-R",TOY_FA,
    		"-r","ref:1-100",
    		TOY_BAM
    	}));
    Assert.assertTrue( tmp.exists());
	}

    @Test
    public void testVcfInjectPed() throws IOException{ 
    	final File tmp=new File(TEST_RESULTS_DIR,"tmp.ped.vcf");
        Assert.assertEquals(0,new VcfInjectPedigree().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"--pedigree",PED01,
        		VCF01
        	}));
        assertIsVcf(tmp);
        Assert.assertTrue(tmp.exists());
        Assert.assertTrue(
        		getVcfHeader(tmp).
        		getOtherHeaderLines()
        		.stream().
        		anyMatch(H->H.getKey().equals("Sample"))
        		);
    	}
    
    
    @Test(dataProvider="all_vcfs")
    public void testVcfToHilbert(final String vcfPath) throws IOException{    
    	if(getVcfHeader(new File(vcfPath)).getNGenotypeSamples()==0) return;
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
    	Assert.assertEquals(0,new VcfToHilbert().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		vcfPath}));
    	Assert.assertTrue( tmp.exists());
    	}
    
    
    @Test(dataProvider="all_vcfs")
    public void testVcfFixIndels(final String vcfPath) throws IOException{    
    	Assert.assertEquals(0,new VCFFixIndels().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		vcfPath}));
    	}
    
    @Test(dataProvider="all_vcfs")
    public void testVcfRemoveUnusedAlt(final String vcfPath) throws IOException{   
    	File tmp = new File(TEST_RESULTS_DIR, "tmp.vcf");
    	Assert.assertEquals(0,new VcfRemoveUnusedAlt().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		vcfPath}));
    	assertIsVcf(tmp);
    	Assert.assertTrue( tmp.delete());
    	}

    
    @Test
    public void testVcfLoopOverGenes() throws IOException{    
    	final String prefix="vlog";
    	File genes = new File(TEST_RESULTS_DIR,"jeter.genes.txt");
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",genes.getPath(),
        		"-p",prefix,
        		"--snpEffNoIntergenic",
        		VCF01
        		}));
    	Assert.assertTrue( genes.exists());
    	
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",TEST_RESULTS_DIR.getPath(),
        		"-g",genes.getPath(),
        		VCF01
        		}));

    	
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",genes.getPath(),
        		"--splitMethod","VariantSlidingWindow",
        		"--variantsWinCount","50",
        		"--variantsWinShift","50",
        		"-p",prefix,
        		"--snpEffNoIntergenic",
        		VCF01
        		}));
    	
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",TEST_RESULTS_DIR.getPath(),
        		"-g",genes.getPath(),
        		VCF01
        		}));
    	Arrays.asList(TEST_RESULTS_DIR.listFiles()).stream().
    		filter(F->F.getName().startsWith(prefix) && F.getName().endsWith(".vcf")).
    		forEach(F->F.delete());
    	}
   
    @Test
    public void testBiostar84452() throws IOException{    
    	File tmp = new File(TEST_RESULTS_DIR, "tmp.bam");
    	Assert.assertEquals(0,new Biostar84452().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		TOY_BAM
        		}));
    	Assert.assertTrue(tmp.delete());
    	}

    @Test
    public void testVcfFilterSo() throws IOException{   
    	File output = new File(TEST_RESULTS_DIR, "jeter.filrerso.vcf");
    	final AnnPredictionParser parser = new AnnPredictionParserFactory().createDefaultParser();
    	final SequenceOntologyTree tree = SequenceOntologyTree.getInstance();
    	String acn = "SO:0001583";
    	final SequenceOntologyTree.Term term =tree.getTermByAcn(acn);
    	final Set<SequenceOntologyTree.Term> terms = term.getAllDescendants();
    	Assert.assertNotNull(term);
    	Assert.assertTrue(terms.size()>1);
    	Assert.assertTrue(terms.contains(term));
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-A",acn,
        		VCF01
        		}));
    	streamVcf(output).forEach(V->{
    		//System.err.println(V.getAttribute("ANN")+" vs "+ terms);
    		Assert.assertTrue(parser.getPredictions(V).stream().
    			flatMap(P->P.getSOTerms().stream()).
				anyMatch(T->terms.contains(T)));
    	});
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-A",acn,
        		"--rmatt",
        		"--invert",
        		VCF01
        		}));
    	streamVcf(output).forEach(V->{
    		Assert.assertFalse(parser.getPredictions(V).stream().
    			flatMap(P->P.getSOTerms().stream()).
				anyMatch(T->terms.contains(T)));
    	});

    	
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-A",acn,
        		"--rmatt",
        		VCF01
        		}));
    	Assert.assertTrue(streamVcf(output).findAny().isPresent());
    	
    	Assert.assertTrue(output.delete());
    	}
    @Test
    public void testVcfRebase() throws IOException{    
    	Assert.assertEquals(0,new VcfRebase().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-R",TOY_FA,
        		TOY_VCF_GZ}));
    	assertIsVcf(JETER_VCF);
    	}
    @Test
    public void testCaseCtrlCanvas() throws IOException{    
    	File tmp = new File(TEST_RESULTS_DIR, "tmp.png");

    	Assert.assertEquals(0,new CaseControlCanvas.Main().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"--pedigree",PED01,
        		VCF01}));
    	final BufferedImage img = ImageIO.read(tmp);
    	Assert.assertTrue(img.getWidth()>0);
    	Assert.assertTrue(img.getHeight()>0);
    	Assert.assertTrue(tmp.delete());
    	}
    
   
    
    @Test
    public void testVcfMoveFiltersToInfo() throws IOException{    

    	Assert.assertEquals(0,new VcfMoveFiltersToInfo().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		VCF01}));
    	assertIsVcf(JETER_VCF);
    	Assert.assertTrue( JETER_VCF.exists());
    	}
    
    @Test(groups={"burden"},dependsOnMethods={"testVcfInjectPed"})
    public void testVcfFilterNotInPedigree() throws IOException{    
    	final File input =new File(TEST_RESULTS_DIR,"tmp.ped.vcf");
    	Assert.assertEquals(0,new VcfFilterNotInPedigree().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		input.getPath()
        		}));
    	assertIsVcf(JETER_VCF);
    	Assert.assertTrue( JETER_VCF.exists());
    	}
    
    @Test(dependsOnMethods={"testVcfInjectPed"})
    public void testVcfFamilies() throws IOException{    
    	final File input =new File(TEST_RESULTS_DIR,"tmp.ped.vcf");
    	Assert.assertEquals(0,new VCFFamilies().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		input.getPath()
        		}));
    	assertIsVcf(JETER_VCF);
    	Assert.assertTrue( JETER_VCF.exists());
    	}
    
    

    
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfBurdenFilterExac() throws IOException{    
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");

	    	Assert.assertEquals(0,new VcfBurdenFilterExac().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		"--discardNotInExac",
	        		"--exac","src/test/resources/ExAC.r1.sites.vep.vcf.gz",
	        		input.getPath()
	        		}));
	    	assertIsVcf(JETER_VCF);
	    	Assert.assertTrue( JETER_VCF.exists());
	    	}
    @Test
    public void testVcfGnomad() throws IOException{    
    		final File manifest = new File(TEST_RESULTS_DIR,"gnomad.manifest");
    		PrintWriter pw = new PrintWriter(manifest);
    		pw.println("exome\t*\tsrc/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz");
    		pw.println("genome\t1\tsrc/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz");
    		pw.flush();
    		pw.close();
    		
	    	Assert.assertEquals(0,new VcfGnomad().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		"-m",manifest.getPath(),
	        		VCF01}));
	    	assertIsVcf(JETER_VCF);
	    	Assert.assertTrue( JETER_VCF.delete());
	    	Assert.assertTrue( manifest.delete());
	    	}
    @Test(dependsOnMethods={"testVcfInjectPed"})
    public void testVcfFisherH() throws IOException{
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
	    	Assert.assertEquals(0,new VcfBurdenFisherH().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		input.getPath()
	        		}));
	    	assertIsVcf(JETER_VCF);
	    	Assert.assertTrue( JETER_VCF.exists());
	    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfFisherV() throws IOException{    
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
	    	Assert.assertEquals(0,new VcfBurdenFisherV().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		input.getPath()
	        		}));
	    	assertIsVcf(JETER_VCF);
	    	Assert.assertTrue( JETER_VCF.exists());
	    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfBurdenMAF() throws IOException{    
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
	    	Assert.assertEquals(0,new VcfBurdenMAF().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		input.getPath()
	        		}));
	    	assertIsVcf(JETER_VCF);
	    	Assert.assertTrue( JETER_VCF.delete());
	    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfBurdenRscriptV() throws IOException{    
			final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
			final File output =new File(TEST_RESULTS_DIR,"tmp.R");
	    	Assert.assertEquals(0,new VcfBurdenRscriptV().instanceMain(new String[]{
	        		"-o",output.getPath(),
	        		input.getPath()
	        		}));
	    	Assert.assertTrue( output.exists());
	    	}
    @Test
    public void testVcfSortOnInfo() throws IOException{    
	    	Assert.assertEquals(0,new SortVcfOnInfo().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		"-T","DP",
	        		TOY_VCF_GZ
	        		}));
	    	assertIsVcf(JETER_VCF);
	    	Assert.assertTrue( JETER_VCF.delete());
	    	}
    
    @Test
    public void testBamToWig() throws IOException{    
		final File output =new File(TEST_RESULTS_DIR,"jeter.wig");
    	Assert.assertEquals(0,new Bam2Wig().instanceMain(new String[]{
        		"-o",output.getPath(),
        		TOY_BAM
        		}));
    	Assert.assertEquals(0,new Bam2Wig().instanceMain(new String[]{
    			"--region","ref:10-30",
    			"-w","3","-s","1","--percentile","MEDIAN",
        		"-o",output.getPath(),
        		TOY_BAM
        		}));
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testVcfCompareCallers() throws IOException{    
		final File output =new File(TEST_RESULTS_DIR,"jeter.zip");
    	Assert.assertEquals(0,new VcfCompareCallers().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-vcf1","VCF1",
        		"-vcf2","VCF2",
        		"--prefix","pfx",
        		TOY_VCF_GZ,TOY_VCF_GZ
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testNgsFilesSummary() throws IOException{   
		final File input =new File(TEST_RESULTS_DIR,"jeter.path.txt");
    	PrintWriter pw=new PrintWriter(input);
    	pw.println(TOY_VCF_GZ);
    	pw.println(TOY_BAM);
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.txt");
    	Assert.assertEquals(0,new NgsFilesSummary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		input.getPath()
        		}));
    	Assert.assertTrue( output.delete());
    	Assert.assertTrue( input.delete());
    	}
    @Test
    public void testFindAllCoverageAtPosition() throws IOException{   
		final File input =new File(TEST_RESULTS_DIR,"jeter.path.txt");
    	PrintWriter pw=new PrintWriter(input);
    	pw.println(TOY_BAM);
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.txt");
    	Assert.assertEquals(0,new FindAllCoverageAtPosition().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-p","ref2:14",
        		input.getPath()
        		}));
    	Assert.assertTrue( output.delete());
    	Assert.assertTrue( input.delete());
    	}
   
    @Test
    public void testMiniCaller() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new MiniCaller().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-R",TOY_FA,
        		TOY_BAM
        		}));
    	assertIsVcf(output);
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testVcfMultiToOne() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfMultiToOne().instanceMain(new String[]{
        		"-o",output.getPath(),
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testSortVcfOnInfo() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new SortVcfOnInfo().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-T","AC",
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testNoCallToHomRef() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfNoCallToHomRef().instanceMain(new String[]{
        		"-o",output.getPath(),
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testVcfAmalgation() throws IOException{   
    	final File tmpXml = File.createTempFile("_tmp", ".xml",TEST_RESULTS_DIR);
    	PrintWriter pw = new PrintWriter(tmpXml);
    	pw.println(
    		"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<config>\n"+
        		"<vcfhead><count>10</count></vcfhead>"+
        		"<vcffilterjdk><filter>ODD</filter><expression>return variant.getStart()%2==0;</expression></vcffilterjdk>"+
        		"<vcftail><count>5</count></vcftail>"+
    		"</config>"
    		);
    	pw.flush();
    	pw.close();
    	
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfXmlAmalgamation().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--xml",tmpXml.getPath(),
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	tmpXml.delete();
    	}
    @Test
    public void testVcfBurdenFilterGenes() throws IOException{   
    	final File tmp = File.createTempFile("_tmp", ".txt",TEST_RESULTS_DIR);
    	PrintWriter pw = new PrintWriter(tmp);
    	pw.println("ISG15");
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfBurdenFilterGenes().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--genes",tmp.getPath(),
        		"--filter","ABCD",
        		"src/test/resources/ExAC.r1.sites.vep.vcf.gz"
        		}));
    	Assert.assertTrue( streamVcf(output).count()>0L);
    	Assert.assertTrue( output.delete());
    	tmp.delete();
    	}
    @Test
    public void testVCFStripAnnotations() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VCFStripAnnotations().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--exclude","INFO/CSQ,ID,FILTER/*",
        		"src/test/resources/ExAC.r1.sites.vep.vcf.gz"
        		}));
    	Assert.assertTrue( streamVcf(output).count()>0L);
    	Assert.assertTrue( streamVcf(output).filter(V->V.hasAttribute("CSQ")).count()==0L);
    	Assert.assertTrue( streamVcf(output).filter(V->V.hasID()).count()==0L);
    	Assert.assertTrue( streamVcf(output).filter(V->V.isFiltered()).count()==0L);
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testGff2kg() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.kg");
    	Assert.assertEquals(0,new Gff2KnownGene().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"src/test/resources/gencode.v19.annotation.gff3"
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testFixVcfMissingGenotypes() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new FixVcfMissingGenotypes().instanceMain(new String[]{
        		"-B",TOY_BAM,
        		"-o",output.getPath(),
        		TOY_VCF_GZ
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testFastqShuffle() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.fq.gz");
    	Assert.assertEquals(0,new FastqShuffle().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"./src/test/resources/SAMPLE1_GATGAATC_L002_R1_001.fastq.gz",
        		"./src/test/resources/SAMPLE1_GATGAATC_L002_R2_001.fastq.gz"
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testPad() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.fq.gz");
    	Assert.assertEquals(0,new PadEmptyFastq().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"./src/test/resources/SAMPLE1_GATGAATC_L002_R1_001.fastq.gz"
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testVcfBed() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
		for(int i=0;i< 2;i++) {
	    	Assert.assertEquals(0,new VCFBed().instanceMain(new String[]{
	        		"-o",output.getPath(),
	        		"-T","VCFBED",
	        		"-e","\"chr\"+bed.get(0)+\":\"+bed.getStart()",
	        		(i==0?"--map":"--bed"),"./src/test/resources/toy.bed.gz",
	        		TOY_VCF_GZ
	        		}));
	    	assertIsVcf(output);
	    	Assert.assertTrue( output.delete());
			}
    	}
    @Test
    public void testCppAlgorithms() throws IOException {
    	final Random rand = new Random(0L);
		final List<Integer> data = new ArrayList<>();
		for(int i=0;i< 50;i++)
			{
			int c = 5+rand.nextInt(20);
			while(c>0) { data.add(i);--c;}
			}
		final int x0 = data.indexOf(40);
		final int x1 = data.indexOf(41);
		Assert.assertTrue(x0>0);
		Assert.assertTrue(x1>x0);
		Assert.assertEquals(x0,Algorithms.lower_bound(data,40));
		Assert.assertEquals(x1,Algorithms.upper_bound(data,40));
		Assert.assertEquals(x0,Algorithms.equal_range(data,40)[0]);
		Assert.assertEquals(x1,Algorithms.equal_range(data,40)[1]);
		Assert.assertEquals(0,Algorithms.lower_bound(data,-1));
		Assert.assertEquals(data.size(),Algorithms.upper_bound(data,60));
    	}
    @Test
    public void testVcfTrap() throws IOException {
    	final File dbFile =new File(TEST_RESULTS_DIR,"chr1.TraPv2.txt");
    	PrintWriter pw =new PrintWriter(dbFile);
    	pw.println("906010\tA\tG\tENSG00000186092\t0.029");
    	pw.println("906010\tA\tG\tENSG00000186092\t1");
    	pw.println("906010\tA\tG\tENSG00000186094\t0");
    	pw.println("906010\tA\tG\tENSG00000186095\t0.0");
    	pw.println("906011\tT\tC\tENSG00000186093\t0.99");
    	pw.flush();
    	pw.close();
    	final File indexFile =new File(TEST_RESULTS_DIR,"chr1.TraPv2.dat");
    	Assert.assertEquals(0,new TrapIndexer().instanceMain(new String[]{
        		"-o",indexFile.getPath(),
        		dbFile.getPath()
        		}));
    	
    	final File manifestFile =new File(TEST_RESULTS_DIR,"jeter.manifest");
    	pw =new PrintWriter(manifestFile);
    	pw.println("1\t"+indexFile.getPath());
    	pw.flush();
    	pw.close();
    	final File outvcf =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfTrap().instanceMain(new String[]{
        		"-o",outvcf.getPath(),
        		"-A","TRAP",
        		"-m",manifestFile.getPath(),
        		"src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz"
        		}));
    	Assert.assertTrue(streamVcf(outvcf).anyMatch(V->V.hasAttribute("TRAP")));
    	Assert.assertTrue(manifestFile.delete());
    	Assert.assertTrue(indexFile.delete());
    	Assert.assertTrue(dbFile.delete());
    	Assert.assertTrue(outvcf.delete());
		}
    
    
    @Test(dataProvider="all_vcfs")
    public void testVcfFileOffsets(final String path) throws IOException {
    	final File vcfFile = new File(path);
    	VcfOffsetsIndexFactory factory= new VcfOffsetsIndexFactory();
    	final File indexFile=factory.indexVcfFile(vcfFile);
    	Assert.assertTrue(indexFile.exists());
    	final VcfList list = VcfList.fromFile(vcfFile);
    	final List<Integer> positions = streamVcf(vcfFile).map(CTX->CTX.getStart()).collect(Collectors.toList());
    	Assert.assertEquals(positions.size(),list.size());
    	Assert.assertNotNull(list.getHeader());
    	for(int i= list.size()-1;i>=0;--i)
    		{
        	Assert.assertNotNull(list.get(i));
        	Assert.assertEquals(list.get(i).getStart(),positions.get(i).intValue());
    		}
    	Assert.assertEquals((long)list.size(),streamVcf(vcfFile).count());
    	list.close();
    	Assert.assertTrue(indexFile.delete());
    	}
    @Test(dataProvider="all_vcfs")
    public void testVcfToRdf(final String path) throws IOException {
    	final File outvcf =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfToRdf().instanceMain(new String[]{
        		"-o",outvcf.getPath(),
        		path
        		}));
    	Assert.assertTrue(outvcf.delete());
		}
    @Test(dataProvider="all_vcfs")
    public void testVcfToSql(final String path) throws IOException {
    	final File outvcf =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfToSql().instanceMain(new String[]{
        		"-o",outvcf.getPath(),
        		path
        		}));
    	Assert.assertTrue(outvcf.delete());
		}

    
    @Test
    public void testBamTile() throws IOException {
    	final File output =new File(TEST_RESULTS_DIR,"jeter.bam");
    	Assert.assertEquals(0,new BamTile().instanceMain(new String[]{
        		"-o",output.getPath(),
        		TOY_BAM
        		}));
    	Assert.assertTrue(output.delete());
		}

    @Test
    public void testRefGenomeFactoryForDAS() throws IOException {
    	final ReferenceGenomeFactory rgf=new ReferenceGenomeFactory();
    	rgf.setBufferSize(10);
    	final ReferenceGenome ref=rgf.openDAS(new URL("http://genome.cse.ucsc.edu/cgi-bin/das/hg19"));
    	Assert.assertTrue(ref.size()>23);
    	ReferenceContig contig = ref.getContig("1");
    	Assert.assertNotNull(contig);
    	Assert.assertEquals(contig.getContig(),"1");
    	Assert.assertEquals(contig.length(),249250621);
    	contig = ref.getContig("M");
    	Assert.assertNotNull(contig);
    	Assert.assertEquals(contig.getContig(),"M");
    	Assert.assertEquals(contig.length(),16571);
    	
    	String dna="gatcacaggtctatcacc";
    	for(int i=0;i< dna.length();++i)
    		{
        	Assert.assertEquals(dna.charAt(i),contig.charAt(i));
    		}
    	dna="cttaaataagacatcacgatg";
    	for(int i=0;i< dna.length();++i)
			{
	    	Assert.assertEquals(dna.charAt(i),contig.charAt(contig.length()-dna.length()+i));
			}
    	ref.close();
    }
    
    @Test
    public void testRefGenomeFactoryForFile() throws IOException {
    	final ReferenceGenomeFactory rgf=new ReferenceGenomeFactory();
    	rgf.setBufferSize(2);
    	final ReferenceGenome ref=rgf.open(TOY_FA);
    	final FastaSequenceReader fsr = new FastaSequenceReader();
    	final List<FastaSequence> seqs = fsr.readAll(new File(TOY_FA));
    	Assert.assertEquals(seqs.size(),ref.size());
    	for(int i=0;i< seqs.size();i++)
    		{
        	Assert.assertEquals(ref.getContig(i).length(),seqs.get(i).length());
        	Assert.assertEquals(ref.getContig(i).getContig(),seqs.get(i).getName());
        	for(int x=0;x<ref.getContig(i).length();++x) {
            	Assert.assertEquals(ref.getContig(i).charAt(x),seqs.get(i).charAt(x));
        		}
    		}

    	ref.close();
		}
    
    @Test
    public void testIlluminaReadName() throws IOException {
    	IlluminaReadName.Parser parser= new IlluminaReadName.Parser();
    	java.util.Optional<IlluminaReadName> n=parser.apply(
    			"@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG");
    	Assert.assertTrue(n.isPresent());
    	Assert.assertEquals(n.get().getInstrument(),"EAS139");
    	Assert.assertEquals(n.get().getRunId(),136);
    	Assert.assertEquals(n.get().getFlowCell(),"FC706VJ");
    	Assert.assertEquals(n.get().getLane(),2);
    	Assert.assertEquals(n.get().getTile(),2104);
    	Assert.assertEquals(n.get().getX(),15343);
    	Assert.assertEquals(n.get().getY(),197393);
    	Assert.assertEquals(n.get().getIndex(),"ATCACG");

    	 n=parser.apply(
     			"@HWUSI-EAS100R:6:73:941:1973#0/1");
     	Assert.assertTrue(n.isPresent());
     	Assert.assertEquals(n.get().getInstrument(),"HWUSI-EAS100R");
     	Assert.assertEquals(n.get().getLane(),6);
     	Assert.assertEquals(n.get().getTile(),73);
     	Assert.assertEquals(n.get().getX(),941);
     	Assert.assertEquals(n.get().getY(),1973);
     	
	    n=parser.apply("@SEQ_ID");
	  	Assert.assertFalse(n.isPresent());
    	}
    
    @Test
    public void testBackLocate() throws IOException {
	    	final File tmp1 = new File(TEST_RESULTS_DIR,"jeter0.txt");
		final PrintWriter pw = new PrintWriter(tmp1);
		pw.println("NOTCH2\tP1090M");
		pw.flush();
		pw.close();
	    	
	    	final File output =new File(TEST_RESULTS_DIR,"jeter1.txt");
	    	Assert.assertEquals(0,new BackLocate().instanceMain(new String[]{
	        		"-o",output.getPath(),
	        		"-R","http://genome.cse.ucsc.edu/cgi-bin/das/hg19/",
	        		tmp1.getPath()
	        		}));
	    	Assert.assertTrue(Files.lines(output.toPath()).anyMatch(S->S.contains("120480548")));
	    	Assert.assertTrue(tmp1.delete());
	    	Assert.assertTrue(output.delete());
		}
    
    @Test
    public void testGoUtils() throws IOException {
	    	final File output = new File(TEST_RESULTS_DIR,"jeter0.txt");
		
	    	Assert.assertEquals(0,new GoUtils().instanceMain(new String[]{
	        		"-A","GO:0005216",
	        		"-go-relations","is_a",
	    			"-o",output.getPath()
	        		}));
	    	assertTableIsConsitent(output, null);
	    	Assert.assertTrue(output.delete());
		}
    
}
