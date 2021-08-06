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
*/
package com.github.lindenb.jvarkit.tools.gvcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Motivation

find regions for running GATK CombineGVCFs in parallel.

## Input

input is a set of path to the indexed g.vcf files
or it's a file with the '.list' suffix containing the path to the g.vcf files

g.vcf files must be indexed if option `-c` is used.

## Output

output is a BED file containing the calleable GVCFs blocks.

## Example

```
$ java -jar dist/findgvcfsblocks.jar --chrom RF11 S1.g.vcf.gz S2.g.vcf.gz S3.g.vcf.gz 
RF11	0	5
RF11	5	12
RF11	12	15
RF11	15	18
RF11	18	20
RF11	20	21
RF11	21	27
RF11	27	28
RF11	28	30
RF11	30	46
(...)
```
END_DOC
*/
@Program(name="findgvcfsblocks",
	description="Find common blocks of calleable regions from a set of gvcfs",
	keywords={"gvcf","gatk","vcf"},
	creationDate="20210806",
	modificationDate="20210806"
	)
public class FindGVCFsBlocks extends Launcher {
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-c","--chrom","--chromosome","--contig"},description="limit to that contig")
	private String the_contig = null;
	@Parameter(names={"-T"},description="temporary directory")
	private Path tmpDir = null;
	@Parameter(names={"--min-size","--block-size"},description="min block size. "+DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_block_size=0;


	private static final Logger LOG = Logger.build(FindGVCFsBlocks.class).make();
	private static final Allele NON_REF = Allele.create("<NON_REF>",false);
	
	/** closeable iterator over a BED file */
	private class Bedterator extends AbstractCloseableIterator<Locatable> {
		private final BufferedReader br;
		private final BedLineCodec codec = new BedLineCodec();
		Bedterator(final Path bedFile) throws IOException {
			this.br = IOUtils.openPathForBufferedReading(bedFile);
			}
		@Override
		protected Locatable advance()
			{
			try {
				String line;
				while((line= br.readLine())!=null) {
					final BedLine bed = codec.decode(line);
					if(bed==null) continue;
					return bed;
					}
				return null;
			} catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
		@Override
		public void close() {
			try {this.br.close();} catch(Throwable err) {}
			}
		}
	/** closeable iterator over a VCF file, returns start-end of blocks */
	private class GVCFVariantIterator extends AbstractCloseableIterator<Locatable>	{
		private final VCFReader vcfFileReader;
		private final VCFHeader header;
		private SAMSequenceDictionary dict;
		private final CloseableIterator<VariantContext> iter;
		private VariantContext first = null;
		private VariantContext prev = null;
		private final VariantContextComparator comparator;
		GVCFVariantIterator(final Path gvcfFile) {
			this.vcfFileReader = VCFReaderFactory.makeDefault().open(gvcfFile,!StringUtils.isBlank(the_contig));
			this.header = this.vcfFileReader.getHeader();
			this.dict = SequenceDictionaryUtils.extractRequired(this.header);
			if(StringUtils.isBlank(the_contig)) {
				this.iter = this.vcfFileReader.iterator();
				}
			else
				{
				final SAMSequenceRecord ssr = this.dict.getSequence(the_contig);
				if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(the_contig, this.dict);
				this.iter = this.vcfFileReader.query(ssr);
				}
			this.comparator = this.header.getVCFRecordComparator();
			}
		@Override
		protected Locatable advance()
			{
			while(this.iter.hasNext()) {
				final VariantContext ctx = iter.next();
				if(this.prev!=null) {
					if(this.comparator.compare(ctx, this.prev) < 0) {
						throw new RuntimeException("Bad order. Got "+ctx+" after "+this.prev);
						}
					}
				if(this.first==null || !this.first.contigsMatch(ctx)) {
					this.first= ctx;
					}
				this.prev = ctx;
				if(ctx.getAlleles().size()!=2) continue;
				if(!ctx.getAlleles().get(1).equals(NON_REF)) continue;
				if(!ctx.hasAttribute(VCFConstants.END_KEY)) continue;
				final SimpleInterval r = new SimpleInterval(ctx.getContig(),this.first.getStart(),ctx.getAttributeAsInt(VCFConstants.END_KEY, -1));
				this.first=null;
				return r;
				}
			return null;
			}
		
		@Override
		public void close() {
			try {this.iter.close();} catch(Throwable err) {}
			try {this.vcfFileReader.close();} catch(Throwable err) {}
			}
		}
		
	
	@SuppressWarnings("resource")
	@Override
	public int doWork(final List<String> args) {
		Path tmpBedFile0 = null;
		Path tmpBedFile1 = null;
		try {
			
			final List<Path> inputs = IOUtils.unrollPaths(args);
			if(inputs.isEmpty())
				{
				LOG.error("input missing");
				return -1;
				}
			if(this.tmpDir==null) {
				this.tmpDir = (this.outputFile==null?IOUtils.getDefaultTempDir():this.outputFile.getParent());				
				}
			
			IOUtil.assertDirectoryIsWritable(this.tmpDir);
			tmpBedFile0 = Files.createTempFile(this.tmpDir, "tmp.", ".bed");
			tmpBedFile1 = Files.createTempFile(this.tmpDir, "tmp.", ".bed");
			SAMSequenceDictionary dict = null;
			final long initMilliSec = System.currentTimeMillis();
			for(int i=0;i< inputs.size();i++) {
				final long startMilliSec = System.currentTimeMillis();
				LOG.info(inputs.get(i)+" "+(i+1)+"/"+inputs.size());
				try(GVCFVariantIterator r0 = new GVCFVariantIterator(inputs.get(i))) {
					if(dict!=null) {
						SequenceUtil.assertSequenceDictionariesEqual(dict, r0.dict);
						}
					else
						{
						dict = r0.dict;
						}
					long count_variants = 0L;
					try(PrintWriter pw= IOUtils.openPathForPrintWriter(tmpBedFile0)) {
						/* first VCF , just convert to bed */
						if(i==0) {
							while(r0.hasNext()) {
								final Locatable loc = r0.next();
								pw.println(loc.getContig()+"\t"+(loc.getStart()-1)+"\t"+loc.getEnd());
								count_variants++;
								}
							}
						/* merge previous bed with current VCF using INFO/END */
						else
							{
							Locatable start0 = null;
							Locatable start1 = null;
							try(Bedterator r1 = new Bedterator(tmpBedFile1)) {
								PeekableIterator<Locatable> peek0 = new PeekableIterator<>(r0);
								PeekableIterator<Locatable> peek1 = new PeekableIterator<>(r1);
								while(peek0.hasNext() && peek1.hasNext()) {
									final Locatable loc0 = peek0.peek();
									final Locatable loc1 = peek1.peek();
									if(!loc0.contigsMatch(loc1)) {
										throw new IllegalStateException("unexpected: not the same contigs "+loc0+" "+loc1);
										}
									if(start0==null) start0 = loc0;
									if(start1==null) start1 = loc1;
									
									
									final int end0 =  loc0.getEnd();
									final int end1 =  loc1.getEnd();
									if(end0 < end1) {
										peek0.next();
										continue;
										}
									else if(end0 > end1) {
										peek1.next();
										continue;
										}
									else { /* end0==end1 */
										pw.println(loc0.getContig()+"\t"+(Math.min(start0.getStart(),start1.getStart())-1)+"\t"+loc0.getEnd());
										count_variants++;
										peek0.next();//consumme
										peek1.next();//consumme
										start0=null;
										start1=null;
										}
									}
								if(peek0.hasNext()) throw new IllegalStateException("peek0 has Next ?");
								if(peek1.hasNext()) throw new IllegalStateException("peek1 has Next ?");
								peek0.close();
								peek1.close();
								}
							}
						pw.flush();
						final long millisecPerVcf  = (System.currentTimeMillis() - initMilliSec)/(i+1L);
												
						LOG.info("N="+count_variants+". That took: "+StringUtils.niceDuration(System.currentTimeMillis() - startMilliSec)+" Remains: "+ StringUtils.niceDuration((inputs.size()-(i+1))*millisecPerVcf));
						}//end writer
					Files.deleteIfExists(tmpBedFile1);
					Files.move(tmpBedFile0,tmpBedFile1);
					}
			
				}
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				try(Bedterator r1 = new Bedterator(tmpBedFile1)) {
					final PeekableIterator<Locatable> peek1 = new PeekableIterator<>(r1);
					while(peek1.hasNext()) {
						Locatable loc = peek1.next();
						while(this.min_block_size>0 && peek1.hasNext()) {
							final Locatable loc2 = peek1.peek();
							if(!loc2.contigsMatch(loc)) break;
							if(CoordMath.getLength(loc.getStart(), loc2.getEnd()) > this.min_block_size) break;
							loc = new SimpleInterval(loc.getContig(),loc.getStart(),loc2.getEnd());
							//consumme loc2
							peek1.next();
							}
						pw.println(loc.getContig()+"\t"+(loc.getStart()-1)+"\t"+loc.getEnd());
						}
					peek1.close();
					}
				pw.flush();
				Files.deleteIfExists(tmpBedFile1);
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}	
		finally {
			if(tmpBedFile0!=null) try { Files.deleteIfExists(tmpBedFile0);} catch(Throwable err) {}
			if(tmpBedFile1!=null) try { Files.deleteIfExists(tmpBedFile1);} catch(Throwable err) {}
			}
		}
	public static void main(final String[] args) {
		new FindGVCFsBlocks().instanceMainWithExit(args);

	}

}
