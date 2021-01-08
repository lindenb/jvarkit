/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.onesamplevcf;

import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NoCloseInputStream;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/*
BEGIN_DOC

## Input

if there is only one input with the '.list' suffix, it is interpreted as a file containing the path to the vcf files

A file with the suffixes '.zip' or '.tar' or '.tar.gz' is interpreted as an archive and all the entries looking like a vcf are extracted.

24 fev 2020: refactored, the input is not anymore sorted. Use bcftools sort

## Example

with zip  and tar

```
$ tar tvfz ~/jeter.tar.gz && unzip -l ~/jeter.zip && java -jar dist/vcfmulti2one.jar ~/jeter.tar.gz ~/jeter.zip | bcftools view - | wc -l
-rw-r--r-- lindenb/lindenb 5805 2019-01-11 18:29 src/test/resources/rotavirus_rf.ann.vcf.gz
-rw-r--r-- lindenb/lindenb 27450 2019-01-11 18:29 src/test/resources/rotavirus_rf.freebayes.vcf.gz
-rw-r--r-- lindenb/lindenb  7366 2019-01-11 18:29 src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz
Archive:  /home/lindenb/jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     7366  2019-01-11 18:29   src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz
     5805  2019-01-11 18:29   src/test/resources/rotavirus_rf.ann.vcf.gz
     3661  2019-01-11 18:29   src/test/resources/rotavirus_rf.vcf.gz
    27450  2019-01-11 18:29   src/test/resources/rotavirus_rf.freebayes.vcf.gz
---------                     -------
    44282                     4 files
4883

```


```bash
$ curl -s "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcfmulti2one.jar  -c -r -a  |\
grep -v '##' |\
grep -E '(CHROM|SAMPLENAME)' | head | verticalize 


>>> 2
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00096;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 2

>>> 3
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00097;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 3

>>> 4
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00099;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 4

>>> 5
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .import java.util.Comparator;

$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00100;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 5

>>> 6
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00102;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 6

>>> 7
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00103;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 7

>>> 8
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00105;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 8

>>> 9
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00106;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 9

>>> 10
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00114;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 10
```



END_DOC
 */
@Program(name="vcfmulti2one",
	biostars=130456,
	description="Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column",
	keywords={"vcf","sample"},
	creationDate="20150312",
	modificationDate="20200224"
	)
public class VcfMultiToOne extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfMultiToOne.class).make();

	@Parameter(names={"-c","-nc","--discard_no_call"},description="discard if variant is no-call")
	private boolean discard_no_call = false;
	@Parameter(names={"-r","-hr","--discard_hom_ref"},description="discard if variant is hom-ref")
	private boolean discard_hom_ref = false;
	@Parameter(names={"-a","--discard_non_available"},description="discard if variant is not available")
	private boolean discard_non_available = false;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--regions"},description="Optional. "+IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class)
	private IntervalListProvider userRegions = null;
	@ParametersDelegate
	private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

	public static final String DEFAULT_VCF_SAMPLE_NAME="SAMPLE";
	public static final String DEFAULT_SAMPLE_TAGID="SAMPLENAME";
	public static final String DEFAULT_SAMPLE_FILETAGID="SAMPLESOURCE";
	public static final String SAMPLE_HEADER_DECLARATION="VcfMultiToOne.Sample";
	
	
	private static class NamedVcfIterator {
		private final VCFIterator delegate;
		private final String name;
		NamedVcfIterator(final VCFIterator delegate,final String name) {
			this.name = name;
			this.delegate = delegate;
		}
	}
	
	private interface VCFIteratorSource extends Closeable /* NOT Iterator. Was a mess with the closed stream */{
	public void open() throws IOException;
	/* get next or return null */
	public NamedVcfIterator next() throws IOException;
	}

	
	
	private class PathsIterator implements VCFIteratorSource {
		final String arg;
		boolean first=false;
		PathsIterator(final String arg) {
			this.arg = arg;
			}
		@Override
		public void open() {
			first=true;
			}
		@Override
		public NamedVcfIterator next() throws IOException{
			if(!first) return null;
			first=false;
			return new NamedVcfIterator(new VCFIteratorBuilder().open(arg),arg.toString());
			}
		@Override
		public void close() throws IOException {}
		}
	
	/* reads all vcf-like entries in a zip archive */
	private class ZipIterator implements VCFIteratorSource {
		private final Path path;
		private InputStream in = null;
		private ZipInputStream zin = null;
		ZipIterator(final Path path){
			this.path= path;
			}
		
		@Override
		public void open()  throws IOException {
			this.in = new BufferedInputStream(Files.newInputStream(path));
			this.zin = new ZipInputStream(this.in);		
			}
		
		@Override
		public NamedVcfIterator next() throws IOException {
				for(;;) {
					final ZipEntry entry = zin.getNextEntry();
					if(entry==null) {
						return null;
						}
					if(entry.isDirectory()) {
						continue;
						}
					if(!FileExtensions.VCF_LIST.stream().anyMatch(X->entry.getName().endsWith(X))) {
						continue;
						}
					/* prevent zip from being closed */
					final InputStream do_not_close_in = new NoCloseInputStream(this.zin);

					return new NamedVcfIterator(
							new VCFIteratorBuilder().open(do_not_close_in),
							path.toString()+"!"+entry.getName()
							);
					}
				}
				
			
		@Override
		public void close() throws IOException {
			if(zin!=null) zin.close();
			if(in!=null) in.close();
			zin=null;
			in=null;
			}
		@Override
		public String toString() {
			return this.path.toString();
			}
		}
	private class TarIterator implements VCFIteratorSource {
		final Path path;
		private TarArchiveInputStream tarin=null;
		private InputStream in=null;
		private GZIPInputStream gzin=null;
		TarIterator(final Path path) {
			this.path=path;
		}
		@Override
		public void open() throws IOException {
			this.in = new BufferedInputStream(Files.newInputStream(this.path));
			if( this.path.getFileName().toString().endsWith(".tar")) {
				this.tarin = new TarArchiveInputStream(this.in);
				this.gzin=null;
				}
			else
				{
				this.gzin=new GZIPInputStream(this.in);
				this.tarin = new TarArchiveInputStream(this.gzin);
				}
			}
		
		@Override
		public NamedVcfIterator next() throws IOException {
				for(;;) {
					final TarArchiveEntry entry = tarin.getNextTarEntry();
					if(entry==null) return null;
					
					if(!tarin.canReadEntryData(entry)) continue;
					
					if(entry.isDirectory()) {
						continue;
						}
					if(!FileExtensions.VCF_LIST.stream().anyMatch(X->entry.getName().endsWith(X))) {
						continue;
						}
					/* prevent zip from being closed */
					final InputStream do_not_close_in = new NoCloseInputStream(tarin);
					return new NamedVcfIterator(
							new VCFIteratorBuilder().open(do_not_close_in),
							path.toString()+"!"+entry.getName()
							);
					}
				}
		@Override
		public void close() throws IOException {
			this.tarin.close();
			if(this.gzin!=null) this.gzin.close();
			this.in.close();
			}
		}
	

	
	public VcfMultiToOne()
		{
		}
	
	
	
	/** general utility for program using VCFMulti2One:
	 *  Extract SampleNames
	 */
	static Set<String> extractSampleNames(final VCFHeader header)
		{
		final List<String> sample_list =header.getSampleNamesInOrder();
		if(sample_list.size()!=1 || !sample_list.get(0).equals(DEFAULT_VCF_SAMPLE_NAME))
			{
			throw new IllegalArgumentException("Not a VCF produced by VcfMultiToOne");
			}
		final Set<String> samples = new TreeSet<String>();
		for(final VCFHeaderLine h:header.getMetaDataInInputOrder())
			{
			if(h.getKey().equals(SAMPLE_HEADER_DECLARATION))
				{
				sample_list.add(h.getValue());
				}
			}
		return samples;
		}
	
	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter  out=null;
		
		
		
		try
			{
			final List<String> paths = IOUtils.unrollStrings2018(args);
			if(paths.isEmpty())
				{
				LOG.error("No vcf provided");
				return -1;
				}
			
			final List<VCFIteratorSource> inputFiles = paths.stream().map(fname->{
				if(fname.endsWith(".tar") || fname.endsWith(".tar.gz")) {
					return new TarIterator(Paths.get(fname));
					}
				else if(fname.endsWith(".zip")) {
					return new ZipIterator(Paths.get(fname));
					}
				else
					{
					return new PathsIterator(fname);
					}
				}).collect(Collectors.toList());

			
			SAMSequenceDictionary dict=null;
			final Set<String> sampleNames=new HashSet<String>();

			final Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
			
			// collect dictionary and meta data
			for(final VCFIteratorSource source:inputFiles) {
				source.open();
				for(;;)
					{
					final	NamedVcfIterator nvi = source.next();
					if(nvi==null) break;
					final VCFIterator in = nvi.delegate;
					final VCFHeader header = in.getHeader();
					final SAMSequenceDictionary dict2= header.getSequenceDictionary();
					if(dict2==null) {
						//nothing
						}
					else if(dict==null)
						{
						dict = dict2;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict2))
						{
						LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict, dict2));
						return -1;
						}
					metaData.addAll(in.getHeader().getMetaDataInInputOrder());
					sampleNames.addAll(in.getHeader().getSampleNamesInOrder());
					in.close();
					}
				source.close();
				}
			
			//addMetaData(metaData);
			metaData.add(new VCFInfoHeaderLine(
					DEFAULT_SAMPLE_TAGID,1,VCFHeaderLineType.String,
					"Sample Name from multi-sample vcf"
					));
			metaData.add(new VCFInfoHeaderLine(
					DEFAULT_SAMPLE_FILETAGID,1,VCFHeaderLineType.String,
					"Origin of sample"
					));
			
			for(final String sample:sampleNames)
				{
				metaData.add(
					new VCFHeaderLine(
					SAMPLE_HEADER_DECLARATION,
					sample));
				}
			
			final VCFHeader h2 = new VCFHeader(
					metaData,
					Collections.singleton(DEFAULT_VCF_SAMPLE_NAME)
					);
			this.recalculator.setHeader(h2);
			JVarkitVersion.getInstance().addMetaData(this, h2);
			
			
			final Predicate<VariantContext> pedicateVariantOverlapUserInterval;
			if(this.userRegions!=null) {
				final  IntervalTreeMap<Boolean> userIntervalTreeMap= new IntervalTreeMap<>();
				this.userRegions.dictionary(dict)
					.skipUnknownContigs()
					.stream()
					.forEach(L->userIntervalTreeMap.put(new Interval(L),Boolean.TRUE));
				pedicateVariantOverlapUserInterval = VC->userIntervalTreeMap.containsOverlapping(VC);
			} else
			{
				pedicateVariantOverlapUserInterval = VC->true;
			}
			
			out= this.writingVariantsDelegate.dictionary(dict).open(this.outputFile);
			out.writeHeader(h2);
			
			for(final VCFIteratorSource source:inputFiles) {
				source.open();
				 for(;;)
					{
					if(out.checkError()) break;
					final	NamedVcfIterator nvi = source.next();
					if(nvi==null) break;
					final VCFIterator in = nvi.delegate;
					while(in.hasNext()) {
						final VariantContext ctx = in.next();
						if(!pedicateVariantOverlapUserInterval.test(ctx)) continue;
						// no genotype
						if(ctx.getNSamples()==0)
							{
							if(!this.discard_no_call)
								{
								final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
								vcb.attribute(DEFAULT_SAMPLE_FILETAGID,nvi.name);
								vcb.genotypes(GenotypeBuilder.createMissing(DEFAULT_VCF_SAMPLE_NAME,2));
								out.add(this.recalculator.apply(vcb.make()));
								}
							continue;
							}
						//loop over samples
						for(int i=0;i< ctx.getNSamples();++i)
							{
							final Genotype g= ctx.getGenotype(i);
							final String sample = g.getSampleName();
							
							if(!g.isCalled() && this.discard_no_call) continue;
							if(!g.isAvailable() && this.discard_non_available) continue;
							if(g.isHomRef() && this.discard_hom_ref) continue;
							
							
							final GenotypeBuilder gb=new GenotypeBuilder(g);
							gb.name(DEFAULT_VCF_SAMPLE_NAME);
							
							
							final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
							vcb.attribute(DEFAULT_SAMPLE_TAGID, sample);
							vcb.attribute(DEFAULT_SAMPLE_FILETAGID,nvi.name);
							
							
							vcb.genotypes(gb.make());
							out.add(this.recalculator.apply(vcb.make()));
							}
						} //end while vcfiterator
					//in.close();
					}
				source.close();
				}
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	
	
	public static void main(final String[] args)
		{
		new VcfMultiToOne().instanceMainWithExit(args);
		}
	}
