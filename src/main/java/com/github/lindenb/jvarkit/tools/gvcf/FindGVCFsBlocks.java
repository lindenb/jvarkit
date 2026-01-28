/*
The MIT License (MIT)
Copyright (c) 2026 Pierre Lindenbaum

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

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
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
	modificationDate="20220401",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class FindGVCFsBlocks extends Launcher {
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-c","--chrom","--chromosome","--contig"},description="limit to that contig")
	private String the_contig = null;
	@Parameter(names={"-T"},description="option was removed",hidden=true)
	private Path _removed_ignore = null;
	@Parameter(names={"--min-size","--block-size"},description="min block size. "+DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_block_size=0;
	@Parameter(names={"--merge-size","-M"},description="merge adjacent blocks distance. "+DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int merge_blocks_distance = 1;

	@Parameter(names={"--lenient"},description="allow strange GVCF blocks that don't end at the same chromosome end.",hidden=true)
	private boolean _deprecated_lenient_processing = false;
	@Parameter(names={"--bed"},description="Restrict output blocks to those overlapping this bed")
	private Path bedPath=null;


	private static final Logger LOG = Logger.of(FindGVCFsBlocks.class);
	private static final Allele NON_REF = Allele.create("<NON_REF>",false);
	
	
	private static class ContigBlocks {
		Integer minPos = null;
		Integer maxPos = null;
		final Set<Integer> end_positions = new TreeSet<>();
		}
	
	
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
	
	
	
	private ContigBlocks scanVcfFile(final ContigBlocks prevBloc ,final SAMSequenceDictionary mainDict,final SAMSequenceRecord ssr, final Path gvcfFile) throws IOException {
			final ContigBlocks blocks = new ContigBlocks();
			if(prevBloc!=null) {
				blocks.minPos = prevBloc.minPos;
				blocks.maxPos = prevBloc.maxPos;
				}
			try(VCFReader vcfFileReader = VCFReaderFactory.makeDefault().open(gvcfFile,true)) {
				final VCFHeader header = vcfFileReader.getHeader();
				final VariantContextComparator comparator = header.getVCFRecordComparator();
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				SequenceUtil.assertSequenceDictionariesEqual(dict, mainDict);
				try(final CloseableIterator<VariantContext> iter = vcfFileReader.query(ssr)) {
					VariantContext prev = null;
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						if(blocks.minPos==null || blocks.minPos.intValue()> ctx.getStart()) {
							blocks.minPos = ctx.getStart();
							}
						if(blocks.maxPos==null || blocks.maxPos.intValue() < ctx.getEnd()) {
							blocks.maxPos = ctx.getEnd();
							}
						if(prev!=null && comparator.compare(ctx,prev) < 0) {
							throw new RuntimeException("Bad order. Got "+ctx+" after "+prev);
							}
						prev = ctx;
						if(ctx.getAlleles().size()!=2) continue;
						if(!ctx.getAlleles().get(1).equals(NON_REF)) continue;
						if(!ctx.hasAttribute(VCFConstants.END_KEY)) continue;
						final int end_pos = ctx.getAttributeAsInt(VCFConstants.END_KEY, -1);
						if(prevBloc==null || prevBloc.end_positions.contains(end_pos)) {
							blocks.end_positions.add(end_pos);
							}
						}
					}
				}
			return blocks;
			}
		
	
	@Override
	public int doWork(final List<String> args) {
		try {
			if (this.bedPath!=null) {
				IOUtil.assertFileIsReadable(this.bedPath);
			}
		
		
			final List<Path> inputs = IOUtils.unrollPaths(args);
			if(inputs.isEmpty())
				{
				LOG.error("input missing");
				return -1;
				}
			if(merge_blocks_distance < 1) {
				LOG.error("bad merge size : "+merge_blocks_distance);
				return -1;
				}
			
			
			final Predicate<Locatable> inCapture;
			if (this.bedPath!=null) {
				final IntervalTreeMap<Boolean> intervalTreeMap;
				try(BedLineReader blr = new BedLineReader(this.bedPath)) {
					intervalTreeMap = blr.toIntervalTreeMap(BED->Boolean.TRUE);
					}
				inCapture = (L)->intervalTreeMap.containsOverlapping(L);
			} else {
				inCapture = (L)->true;
			}
			
			if(this.outputFile!=null) {
				final String fname = this.outputFile.getFileName().toString();
				if(!fname.endsWith(FileExtensions.INTERVAL_LIST) && !fname.endsWith(FileExtensions.COMPRESSED_INTERVAL_LIST)) {
					LOG.error("Output should end with "+ FileExtensions.INTERVAL_LIST+" or "+FileExtensions.COMPRESSED_INTERVAL_LIST);
					return -1;
					}
				}
			
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(inputs.get(0));
			if(!StringUtils.isBlank(the_contig) && dict.getSequence(the_contig)==null) throw new JvarkitException.ContigNotFoundInDictionary(the_contig, dict);
			
			
			try(IntervalListWriter w = new IntervalListWriter(this.outputFile,dict)) {
				/** loop over chromosomes */
				for(final SAMSequenceRecord ssr :dict.getSequences()){
					ContigBlocks mainBlock = null;
					if(!StringUtils.isBlank(the_contig) && !ssr.getContig().equals(the_contig)) continue;
					final long initMilliSec = System.currentTimeMillis();
					for(int i=0;i< inputs.size();i++) {
						final long startMilliSec = System.currentTimeMillis();
						LOG.info(inputs.get(i)+" "+(i+1)+"/"+inputs.size());
						mainBlock =  scanVcfFile(mainBlock,dict,ssr, inputs.get(i));
						final long count_variants = mainBlock.end_positions.size();
						final long millisecPerVcf  = (System.currentTimeMillis() - initMilliSec)/(i+1L);
						LOG.info("N="+count_variants+". That took: "+StringUtils.niceDuration(System.currentTimeMillis() - startMilliSec)+" Remains: "+ StringUtils.niceDuration((inputs.size()-(i+1))*millisecPerVcf));
						}
					/* nothing was seen */
					if(mainBlock.minPos==null) continue;
					if(mainBlock.end_positions.isEmpty()) {
						final Locatable loc = new SimpleInterval(ssr.getContig(),mainBlock.minPos,mainBlock.maxPos);
						if(inCapture.test(loc)) {
							w.add(loc);
							}
						continue;
						}
					final List<Locatable> intervals = new ArrayList<>(mainBlock.end_positions.size()+1);
					
					final int pos1 = mainBlock.end_positions.stream().mapToInt(P->P.intValue()).min().getAsInt();
					if(mainBlock.minPos< pos1 ) {
						final Locatable loc = new SimpleInterval(ssr.getContig(),mainBlock.minPos,pos1-1);
						if(inCapture.test(loc)) {
							intervals.add(loc);
							}
						}
					Integer prev=null;
					for(Integer p: mainBlock.end_positions) {
						if(prev!=null) {
							final Locatable loc = new SimpleInterval(ssr.getContig(),prev+1,p);
							if(inCapture.test(loc)) {
								intervals.add(loc);
								}
							}
						prev=p;
						}
					
					final int pos2 = mainBlock.end_positions.stream().mapToInt(P->P.intValue()).max().getAsInt();
					if(mainBlock.maxPos> pos2 ) {
						final Locatable loc = new SimpleInterval(ssr.getContig(),pos2+1,mainBlock.maxPos);
						if(inCapture.test(loc)) {
							intervals.add(loc);
							}
						}
					Collections.sort(intervals,(A,B)->Integer.compare(A.getStart(), B.getStart()));
					while(!intervals.isEmpty()) {
						Locatable loc = intervals.remove(0);
						while(this.min_block_size>0 && !intervals.isEmpty() ) {
							final Locatable loc2 = intervals.get(0);
							if(!loc2.withinDistanceOf(loc, merge_blocks_distance)) break;
							if(CoordMath.getLength(loc.getStart(), loc2.getEnd()) > this.min_block_size) break;
							//consumme loc2
							intervals.remove(0);
							loc = new SimpleInterval(loc.getContig(),loc.getStart(),loc2.getEnd());
							}
						if(inCapture.test(loc)) {
							w.add(loc);
							}
						}
					}
				}
				
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			}
		}
	public static void main(final String[] args) {
		new FindGVCFsBlocks().instanceMainWithExit(args);
		}

}
