package com.github.lindenb.jvarkit.tools.tests;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.Vector;
import java.util.function.Predicate;
import java.util.stream.Collectors;
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
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class TestSupport {
	public final Random random = new Random();
	protected final String SRC_TEST_RESOURCE="./src/test/resources";
	private final List<Path> deletePathsAtExit = new Vector<>();
	
	public Predicate<String> bamHasIndex = (F)->{
		final Path p = Paths.get(F);
		if(!Files.exists(p)) return false;
		try (SamReader sr = SamReaderFactory.makeDefault().open(p)) {
			return sr.hasIndex();
			}
		catch(Exception err) {
			return false;
			}
		};
	public Predicate<String> vcfhasIndex = (F)->{
		final Path p = Paths.get(F);
		if(!Files.exists(p)) return false;
		try (VCFFileReader sr = new VCFFileReader(p,true)) {
			return sr.isQueryable();
			}
		catch(Throwable err) {
			return false;
			}
		};
	
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

	
	
	
	public synchronized Path deleteOnExit(final Path f) {
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
	
	public void assertIsFastq(final Path f) {
		Exception err=null; 
		FastqReader fqr =null;
		try {
			fqr = new FastqReader(f.toFile());
			Iterator<FastqRecord> iter =fqr.iterator();
			while(iter.hasNext()) iter.next();
		}catch (final Exception e) {
			err= e;
		} finally {
			CloserUtil.close(fqr);
		}
		Assert.assertNull(err,"file "+f+" should be fastq : "+err);
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
	public Path createRandomPedigreeFromFile(final String vcfFile) throws IOException {
		class TmpSample {
			String fam="F1";
			String name;
			String father="0";
			String mother="0";
			int sex=0;
			int status=0;
			}
		final Path vcfIn = Paths.get(vcfFile);
		final VCFFileReader r= new VCFFileReader(vcfIn,false);
		final VCFHeader header=r.getFileHeader();
		r.close();
		final List<String> samples= new ArrayList<>(header.getSampleNamesInOrder());
		final List<TmpSample> ped= new ArrayList<>();
		if(samples.isEmpty()) return null;
		while(!samples.isEmpty())
			{
			final TmpSample indi = new TmpSample();
			indi.name = samples.remove(0);
			indi.fam = ped.stream().
					filter(P->P.father.equals(indi.name) || P.mother.equals(indi.name)).
					map(P->P.fam).findFirst().orElse("F"+this.random.nextInt()
					);
			if(ped.stream().anyMatch(P->P.father.equals(indi.name))) {
				indi.sex=1;
				}
			else if(ped.stream().anyMatch(P->P.mother.equals(indi.name))) {
				indi.sex=2;
				}
			else
				{
				indi.sex = this.random.nextBoolean()?1:2;
				}
			if(random.nextBoolean())
				{
				final List<String> remain=new ArrayList<>(samples);
				if(!remain.isEmpty())
					{
					indi.father = remain.remove(0);
					}
				else
					{
					indi.father = "0";
					}
				if(!remain.isEmpty())
					{
					indi.mother = remain.remove(0);
					}
				else
					{
					indi.mother = "0";
					}
				}
			else
				{
				indi.father = "0";
				indi.mother = "0";
				}
			
			indi.status = this.random.nextInt(2);
			ped.add(indi);
			}
		final Path pedFile = createTmpPath(".ped");
		final PrintWriter pw  = new PrintWriter(Files.newBufferedWriter(pedFile));
		for(TmpSample ts:ped)
			{
			pw.print(ts.fam);
			pw.print('\t');
			pw.print(ts.name);
			pw.print('\t');
			pw.print(ts.father);
			pw.print('\t');
			pw.print(ts.mother);
			pw.print('\t');
			pw.print(ts.sex);
			pw.print('\t');
			pw.print(ts.status);
			pw.println();
			}
		pw.flush();
		pw.close();
		return pedFile;
		}

	public Object[][] toArrayArray(Stream<Object[]> st) {
		final List<Object[]> L=  st.collect(Collectors.toList());
		return L.toArray(new Object[L.size()][]);
		}

	public Path addClippingToBam(final Path bamFile) throws IOException {
		final String bases = "ATGC";
		Path clippedBam = this.createTmpPath(".bam");
		SamReader sr = SamReaderFactory.makeDefault().open(bamFile);
		SAMFileHeader inHeader = sr.getFileHeader();
		boolean createIndex = sr.hasIndex() && inHeader.getSortOrder().equals(SortOrder.coordinate);
		if(createIndex) {
				final Path parent = bamFile.getParent();
				final String baiFname = IOUtil.basename(clippedBam.toFile()) + ".bai";// deprecated? BAMIndex.BAI_INDEX_SUFFIX;
				final Path bai = parent.resolve(baiFname);
			this.deletePathsAtExit.add(bai);
			}
		final SAMFileWriter w=new SAMFileWriterFactory().
				setCreateIndex(createIndex).
				makeBAMWriter(inHeader, true, clippedBam);
		sr.iterator().stream().map(R->{
			if(R.getReadUnmappedFlag() || R.getCigar()==null) return R;
			if(R.getCigar().isClipped()) return R;
			if(R.getBaseQualities().equals(SAMRecord.NULL_QUALS)) return R;
			if(R.getBaseQualityString().equals(SAMRecord.NULL_QUALS_STRING)) return R;
			for(int side=0;side<2;side++) {
				
				final String cigar;
				boolean hard = this.random.nextBoolean();
				final int clipLen = 1+this.random.nextInt(100);
				final StringBuilder seq = new StringBuilder();
				final StringBuilder qual = new StringBuilder();
				
				if(hard) {
					cigar = String.valueOf(clipLen)+"H";
					} 
				else
					{
					cigar = String.valueOf(clipLen)+"S";
					
					
					for(int x=0;x<clipLen;++x) {
						seq.append(bases.charAt(this.random.nextInt(bases.length())));
						qual.append("#");
						}
					}
				if(side==0) {
					R.setReadString(seq.toString()+R.getReadString());
					R.setBaseQualityString(qual.toString()+R.getBaseQualityString());
					R.setCigarString(cigar+R.getCigarString());
					}
				else
					{
					R.setCigarString(R.getCigarString()+cigar);
					R.setReadString(R.getReadString()+seq.toString());
					R.setBaseQualityString(R.getBaseQualityString()+qual.toString());
					}
				}
			return R;
			}).forEach(R->w.addAlignment(R));
		w.close();
		sr.close();
		assertIsValidBam(bamFile);
		return bamFile;
		}
	
	public Object[][] combine2(Stream<?> col1,Stream<?> col2) {
		List<Object[]> L = new ArrayList<>();
		List<Object> L1 = col1.map(O->Object.class.cast(O)).collect(Collectors.toList());
		List<Object> L2 = col2.map(O->Object.class.cast(O)).collect(Collectors.toList());
		for(int x=0;x< L1.size();++x) {
			for(int y=0;y< L2.size();++y) {
				Object row[] = new Object[] {L1.get(x),L2.get(y)};
				L.add(row);
				}
			}
		
		return toArrayArray(L.stream());
	}
	public List<Interval> randomIntervalsFromDict(final Path dictFile,int n,int maxLen) throws IOException{
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(dictFile);
		final List<Interval> rgns = new ArrayList<>(n);
		if(dict==null) return rgns;
		while(n>0)
			{
			final SAMSequenceRecord ssr = dict.getSequence(random.nextInt(dict.size()));
			final int L = 1 + random.nextInt(Math.min(maxLen,ssr.getSequenceLength()-1));
			final int start = 1 + random.nextInt(ssr.getSequenceLength() -L);
			rgns.add(new Interval(ssr.getSequenceName(), start, start+L));
			n--;
			}
		return rgns;
		}
}
