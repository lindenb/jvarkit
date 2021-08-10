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
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.IntervalList.IntervalListCodec;
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

input is a set of path to the indexed g.vcf files or a picard-style interval file generated with a previous invocation of findgvcfsblocks.jar with one sample.
or it's a file with the '.list' suffix containing the path to the g.vcf files/interval files

g.vcf files must be indexed if option `-c` is used.

## Output

output is a picard-style **Interval** file containing the calleable GVCFs blocks.

## Example

```
$ java -jar dist/findgvcfsblocks.jar --min-size 100 --chrom RF11 S1.g.vcf.gz S2.g.vcf.gz S3.g.vcf.gz 
@HD	VN:1.6
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@CO	findgvcfsblocks. compilation:20210807160340 githash:b442941 htsjdk:2.24.1 date:20210807160354. cmd:--min-size 100 --chrom RF11 S1.g.vcf.gz S2.g.vcf.gz S3.g.vcf.gz
RF11	1	95	+	.
RF11	96	182	+	.
RF11	183	237	+	.
RF11	238	428	+	.
RF11	429	528	+	.
RF11	529	628	+	.
RF11	629	666	+	.
(...)
```

## Nextflow

```
(...)
Channel.fromPath(params.reference+".fai").
	splitCsv(header: false,sep:'\t',strip:true).
	filter{T->T[0].matches("(chr)?[0-9XY]+")}.
	map{T->[T[0]]}.
	set{each_contig}

process gvcflists {
executor "local"
output:
	path("gvcfs.list") into (gvcfs_list1,gvcfs_list2)
script:
"""
test -s "${params.gvcfs}"

SQRT=`awk 'END{z=sqrt(NR); print (z==int(z)?z:int(z)+1);}' "${params.gvcfs}"`

split -a 9 --additional-suffix=.list --lines=\${SQRT} "${params.gvcfs}" chunck.

find \${PWD} -type f -name "chunck.*.list"  > gvcfs.list

"""
}


gvcfs_list2.splitCsv(header: false,sep:'\t',strip:true).
	map{T->T[0]}.
	set{one_gvcf_list}


process findBlocks {
tag "${file(vcfs).name} ${contig}"
memory "5g"
input:
	tuple vcfs,contig from one_gvcf_list.combine(each_contig)
output:
	tuple contig,path("block.interval_list.gz") into blocks_interval
script:
"""
module load jvarkit/0.0.0

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/findgvcfsblocks.jar \
	--contig "${contig}" \
	-T . \
	-o "block.interval_list.gz" \
	"${vcfs}"

test -s block.interval_list.gz
"""
}


process findCommonBlocks {
tag "${contig} N=${L.size()}"
memory "10g"
input:
	tuple contig,L from blocks_interval.groupTuple()
output:
	path("block.interval_list") into blocks_common
script:
"""
module load jvarkit/0.0.0

cat << EOF > jeter.list
${L.join("\n")}
EOF

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/findgvcfsblocks.jar \
	--contig "${contig}" \
	-T . \
	--block-size "${params.blocksize}" \
	-o "block.interval_list" \
	jeter.list

rm jeter.list

test -s block.interval_list
"""
}


blocks_common.splitCsv(header: false,sep:'\t',strip:true).
	filter{T->!T[0].startsWith("@")}.
	map{T->T[0]+":"+T[1]+"-"+T[2]}.
	set{intervals}

process genotypeRegion {
tag "${region}"
cache "lenient"
errorStrategy "finish"
afterScript "rm -rf TMP"
memory "15g"
input:
	val region from intervals
	path gvcfs from gvcfs_list1
output:
	path("genotyped.vcf.gz")  into genotyped
	path("genotyped.vcf.gz.tbi")
script:
(...)

```



END_DOC
*/
@Program(name="findgvcfsblocks",
	description="Find common blocks of calleable regions from a set of gvcfs",
	keywords={"gvcf","gatk","vcf"},
	creationDate="20210806",
	modificationDate="20210809"
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
	
	private class IntervalListWriter implements Closeable {
		private final BufferedWriter w;
		IntervalListWriter(final Path p,final SAMSequenceDictionary dict) throws IOException {
			w = p==null?
				new BufferedWriter(new OutputStreamWriter(System.out, "UTF-8")):
				IOUtil.openFileForBufferedUtf8Writing(p)
				;
			final SAMFileHeader header = new SAMFileHeader(dict);
			header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
			JVarkitVersion.getInstance().addMetaData(FindGVCFsBlocks.this, header);
			final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
	        codec.encode(w, header);
			}
		public void add(final Locatable loc) throws IOException{
			w.write(loc.getContig()+"\t"+loc.getStart()+"\t"+loc.getEnd()+"\t+\t.");//yes, not bed, interval list so no -1
			w.newLine();
			}
		@Override
		public void close() throws IOException {
			w.flush();
			w.close();
			}
		
		}
	
	
	private abstract class GVCFOrIntervalIterator extends AbstractCloseableIterator<Locatable>	{
		public abstract SAMSequenceDictionary getSAMSequenceDictionary();
		}
	
	/** closeable iterator over a VCF file, returns start-end of blocks */
	private class IntervalListIterator extends GVCFOrIntervalIterator	{
		private BufferedReader br;
		private final SAMSequenceDictionary dict;
		private String line;
	 final IntervalListCodec intervalListCodec;
		IntervalListIterator(final Path intervalFile) throws IOException {
			this.br =  IOUtil.openFileForBufferedReading(intervalFile);
			final StringBuilder builder = new StringBuilder();
			try {
				 while ((this.line = br.readLine()) != null) {
		             if (this.line.startsWith("@")) {
		                 builder.append(this.line).append('\n');
		             	} else {
		                 break;
		             	}
					 	}
					}
			catch(final Throwable err) {
					throw new IOException(err);
					}
			if (builder.length() == 0) {
				throw new IllegalStateException("Interval list file must contain header. ");
				}

			final BufferedLineReader headerReader = BufferedLineReader.fromString(builder.toString());
			final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
			final SAMFileHeader samfileHeader = codec.decode(headerReader,intervalFile.toString());
			this.dict = SequenceDictionaryUtils.extractRequired(samfileHeader);
			this.intervalListCodec = new IntervalListCodec(this.dict);
			}
		@Override
		protected Locatable advance()
			{
			try {
				/* first line was read in the header */
				while(this.line!=null) {
					final Locatable loc = this.intervalListCodec.decode(this.line);
					// now prepare line for next iteration
					this.line = this.br.readLine();
					if(loc==null || (!StringUtils.isBlank(the_contig) && !loc.getContig().equals(the_contig))) continue;
					return loc;
					} 
				} catch(final IOException err) {
					throw new RuntimeIOException(err);
				}
			return null;
			}
		@Override
		public SAMSequenceDictionary getSAMSequenceDictionary() {
			return this.dict;
			}
		@Override
		public void close() {
			this.line = null;
			try{if(this.br!=null) br.close();} catch(IOException err){}
			}
		}
	
	/** closeable iterator over a VCF file, returns start-end of blocks */
	private class GVCFVariantIterator extends GVCFOrIntervalIterator	{
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
		public SAMSequenceDictionary getSAMSequenceDictionary()
			{
			return this.dict;
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
	
	private GVCFOrIntervalIterator openInput(final Path path) throws IOException {
		final String fname = path.getFileName().toString();
		if(FileExtensions.VCF_LIST.stream().anyMatch(S->fname.endsWith(S))) {
			return new GVCFVariantIterator(path);
			}
		else if(fname.endsWith(FileExtensions.INTERVAL_LIST) || fname.endsWith(FileExtensions.COMPRESSED_INTERVAL_LIST)) {
			return new IntervalListIterator(path);
			}
		else
			{
			throw new IOException("unknown file extension : "+path+" not "+FileExtensions.INTERVAL_LIST+"/"+FileExtensions.COMPRESSED_INTERVAL_LIST+"/"+String.join("/", FileExtensions.VCF_LIST));
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
			if(this.outputFile!=null) {
				final String fname = this.outputFile.getFileName().toString();
				if(!fname.endsWith(FileExtensions.INTERVAL_LIST) && !fname.endsWith(FileExtensions.COMPRESSED_INTERVAL_LIST)) {
					LOG.error("Output should end with "+ FileExtensions.INTERVAL_LIST+" or "+FileExtensions.COMPRESSED_INTERVAL_LIST);
					return -1;
					}
				}
			
			if(this.tmpDir==null && this.outputFile!=null) {
				this.tmpDir = this.outputFile.getParent();				
				}
			
			if(this.tmpDir==null) {
				this.tmpDir = IOUtils.getDefaultTempDir();
				}
			IOUtil.assertDirectoryIsWritable(this.tmpDir);
			tmpBedFile0 = Files.createTempFile(this.tmpDir, "tmp.", ".bed");
			tmpBedFile1 = Files.createTempFile(this.tmpDir, "tmp.", ".bed");
			SAMSequenceDictionary dict = null;
			final long initMilliSec = System.currentTimeMillis();
			for(int i=0;i< inputs.size();i++) {
				final long startMilliSec = System.currentTimeMillis();
				LOG.info(inputs.get(i)+" "+(i+1)+"/"+inputs.size());
				try(GVCFOrIntervalIterator r0 = openInput(inputs.get(i))) {
					if(dict!=null) {
						SequenceUtil.assertSequenceDictionariesEqual(dict, r0.getSAMSequenceDictionary());
						}
					else
						{
						dict = r0.getSAMSequenceDictionary();
						}
					long count_variants = 0L;
					try(IntervalListWriter pw= new IntervalListWriter(tmpBedFile0,dict)) {
						/* first VCF , just convert to bed */
						if(i==0) {
							while(r0.hasNext()) {
								final Locatable loc = r0.next();
								pw.add(loc);
								count_variants++;
								}
							}
						/* merge previous bed with current VCF using INFO/END */
						else
							{
							Locatable start0 = null;
							Locatable start1 = null;
							try(IntervalListIterator r1 = new IntervalListIterator(tmpBedFile1)) {
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
										pw.add(new SimpleInterval(loc0.getContig(),(Math.min(start0.getStart(),start1.getStart())),loc0.getEnd()));
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
						final long millisecPerVcf  = (System.currentTimeMillis() - initMilliSec)/(i+1L);
												
						LOG.info("N="+count_variants+". That took: "+StringUtils.niceDuration(System.currentTimeMillis() - startMilliSec)+" Remains: "+ StringUtils.niceDuration((inputs.size()-(i+1))*millisecPerVcf));
						}//end writer
					Files.deleteIfExists(tmpBedFile1);
					Files.move(tmpBedFile0,tmpBedFile1);
					}
			
				}
			try(IntervalListWriter w = new IntervalListWriter(this.outputFile,dict)) {
				try(IntervalListIterator r1 = new IntervalListIterator(tmpBedFile1)) {
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
						w.add(loc);
						}
					peek1.close();
					}
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
