package com.github.lindenb.jvarkit.tools.tests;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Vector;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import java.util.zip.ZipInputStream;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.testng.Assert;
import org.xml.sax.helpers.DefaultHandler;

import com.github.lindenb.jvarkit.dict.ReferenceRegistry;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class TestSupport {
	protected final String SRC_TEST_RESOURCE="./src/test/resources";
	private final List<Path> deletePathsAtExit = new Vector<>();
	
	private ReferenceRegistry refCatalog =new ReferenceRegistry() {
			@Override
			public Optional<Path> getReferenceByName(String name) {
				if(name.equals("rf")) return Optional.of(Paths.get(resource("rotavirus_rf.fa")));
				if(name.equals("toy")) return Optional.of(Paths.get(resource("toy.fa")));
				return ReferenceRegistry.getDefault().getReferenceByName(name);
				}
			
			@Override
			public Optional<Path> getReferenceByDictionary(final SAMSequenceDictionary dict) {
				if(dict.size()==11 && dict.getSequence(0).getSequenceLength()==3_302 && dict.getSequence(0).getSequenceName().equals("RF01")) return getReferenceByName("rf");
				if(dict.size()==2 && 
						dict.getSequence(0).getSequenceLength()==45 &&
						dict.getSequence(0).getSequenceName().equals("ref") &&
						dict.getSequence(1).getSequenceLength()==40 &&
						dict.getSequence(1).getSequenceName().equals("ref2")) 
					return getReferenceByName("toy");
				return ReferenceRegistry.getDefault().getReferenceByDictionary(dict);
			}
		};
	
	public ReferenceRegistry getReferenceRegistry() {
		return this.refCatalog;
	}
		
	public String resource( String fname) {
		if(!fname.startsWith(File.separator)) fname=File.separator+fname;
		return SRC_TEST_RESOURCE+fname;
	}
	
	public Path createTmpPath(final String suffix) throws IOException {
		return deleteOnExit(Files.createTempFile("tmp.", suffix));
	}

	
	protected synchronized Path deleteOnExit(final Path f) {
		if(f!=null) this.deletePathsAtExit.add(f);
		return f;
		}
	
	public synchronized void removeTmpFiles() {
		for(final Path f:this.deletePathsAtExit) try {
			Files.delete(f);
			} catch(IOException err) {}
		System.err.println("[TEST] end "+this.getClass().getSimpleName());
		}
	
	public Stream<VariantContext> variantStream(final Path vcfFile ) {
		final VCFFileReader r = new VCFFileReader(vcfFile,false);
		final CloseableIterator<VariantContext> iter = r.iterator();
		return iter.stream().onClose(()->{CloserUtil.close(iter);CloserUtil.close(r);});
		}

	public void assertIsNotEmpty(final Path f) throws IOException {
		Assert.assertNotNull(f,"File is null");
		Assert.assertFalse(Files.size(f)==0L,"file "+f+" is empty");
		}
	
	public void assertIsValidBam(final Path bamFile) throws IOException {
		final SamReader sr= SamReaderFactory.makeDefault().
				validationStringency(ValidationStringency.SILENT).
				open(bamFile);
		sr.iterator().stream().count();
		sr.close();
		}

	
	public void assertIsXml(final Path f) {
		Exception err=null;
		;
		try {
			Assert.assertTrue(Files.exists(f), "file "+f+" should exist");
			final SAXParser sax=SAXParserFactory.newInstance().newSAXParser();
			try(InputStream in = Files.newInputStream(f))
				{
				sax.parse(in, new DefaultHandler());
				}
			
		}catch (final Exception e) {
			err= e;
			}
		Assert.assertNull(err,"file "+f+" should be xml : "+err);
		}
	public long assertIsBed(final Path f) {
		Exception err=null;
		BufferedReader r = null;
		try {
			r = IOUtils.openPathForBufferedReading(f);
			final BedLineCodec codec = new BedLineCodec();
			return r.lines().map(L->codec.decode(L)).filter(B->B!=null).count();
		}catch (final Exception e) {
			err= e;
		} finally {
			CloserUtil.close(r);
		}
		Assert.assertNull(err,"file "+f+" should be BED : "+err);
		return 0L;
	}

	
	
	public boolean assertZip(final Path f) {
		InputStream is = null;
		try {
			is= Files.newInputStream(f);
			ZipInputStream zin = new ZipInputStream(is);
			zin.close();
			return true;
		} catch(IOException err) {
			return false;
		}
		finally {
			CloserUtil.close(is);
		}
	}
	
	public void assertTsvTableIsConsitent(final Path f,Predicate<String> ignoreLine) {
		final CharSplitter tab=CharSplitter.TAB;
		BufferedReader r=null;
		try {
			Assert.assertTrue(Files.exists(f), "file "+f+" should exist");
			r= IOUtils.openPathForBufferedReading(f);
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

	
	public void assertIsVcf(final Path f) {
		Exception err=null; 
		try {
			variantStream(f).
				flatMap(V->V.getGenotypes().stream()).
				flatMap(G->G.getAlleles().stream());
			
		}catch (final Exception e) {
			err= e;
		}
		Assert.assertNull(err,"file "+f+" should be vcf : "+err);
		}
	
	public Stream<String> allSamOrBams() {
		return Arrays.asList(
		"FAB23716.nanopore.bam",
		"HG02260.transloc.chr9.14.bam",
		"S1.bam", "S2.bam", "S3.bam", "S4.bam", "S5.bam",
		"toy.bam", "toy.sam"
		).stream().map(S->resource(S))
		;
	}
	
	public Stream<String> allVcfOrBcf() {
		return Arrays.asList(
				"ExAC.r1.sites.vep.vcf.gz","gnomad.exomes.r2.0.1.sites.vcf.gz",
				"gnomad.genomes.r2.0.1.sites.1.vcf.gz",
				"rotavirus_rf.ann.vcf.gz",
				"rotavirus_rf.freebayes.vcf.gz",
				"rotavirus_rf.unifiedgenotyper.vcf.gz",
				"rotavirus_rf.vcf.gz",
				"S1.vcf.gz","S2.vcf.gz","S3.vcf.gz",
				"S4.vcf.gz","S5.vcf.gz",
				"test_vcf01.vcf","toy.vcf.gz").
			stream().
			map(S->resource(S))
		;
	}
	
	public Stream<String> allFasta() {
		return Arrays.asList("rotavirus_rf.fa","toy.fa").
				stream().
				map(S->resource(S))
			;
		}
	
	public boolean isBamIndexed(final Path b) {
		try {
			SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			try(SamReader r=srf.open(b))
				{
				return r.hasIndex();
				}
		} catch(Exception err) {
			return false;
		} 
	}
	
	public Stream<String> allIndexedBams() {
		return allSamOrBams().
				filter(S->S.endsWith(".bam")).
				filter(S->isBamIndexed(Paths.get(S)))
				;
	}
	public long wc(final Path f) throws IOException
		{
		if(f.getFileName().toString().endsWith(".bam") || f.getFileName().toString().endsWith(".sam")) {
			final SamReader sr = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.SILENT).
					open(f);
			final SAMRecordIterator iter =sr.iterator();
			final long n= iter.stream().count();
			iter.close();
			sr.close();
			return n;
			}
		
		final BufferedReader br = IOUtils.openPathForBufferedReading(f);
		final long n=br.lines().count();
		br.close();
		return n;
		}
	public Path sortBamOnQueryName(final Path bamFile,final Predicate<SAMRecord> pred) throws IOException {
		final Path sortedBam = this.createTmpPath(".bam");
		SamReader sr = SamReaderFactory.makeDefault().open(bamFile);
		SAMFileHeader outHeader= sr.getFileHeader().clone();
		outHeader.setSortOrder(SortOrder.queryname);
		SAMFileWriter w=new SAMFileWriterFactory().makeBAMWriter(outHeader, false, sortedBam);
		sr.iterator().stream().filter(R->pred==null?true:pred.test(R)).forEach(R->w.addAlignment(R));
		w.close();
		sr.close();
		return sortedBam;
		}

}
