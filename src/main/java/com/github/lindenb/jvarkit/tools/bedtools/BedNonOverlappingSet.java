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
package com.github.lindenb.jvarkit.tools.bedtools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
/**
BEGIN_DOC

## Motivation

GATK DepthOfCoverage merge overlapping segments (see https://gatkforums.broadinstitute.org/gatk/discussion/1865/ ). I wan't to get the coverage for a set of overlapping windows.

## EXAMPLES

### Example

```bash
$ awk '{printf("%s\t0\t%s\n",$1,$2);}' src/test/resources/rotavirus_rf.fa.fai |\
bedtools makewindows  -w 500 -s 100 -b - |\
java -jar dist/bednonoverlappingset.jar  -o out.zip -m tmp.manifest

[INFO][BedNonOverlappingSet]saving tmp.00001.bed 44
[INFO][BedNonOverlappingSet]saving tmp.00002.bed 39
[INFO][BedNonOverlappingSet]saving tmp.00003.bed 37
[INFO][BedNonOverlappingSet]saving tmp.00004.bed 36
[INFO][BedNonOverlappingSet]saving tmp.00005.bed 33

$ head -n 2 tmp.000*.bed
==> tmp.00001.bed <==
RF01	0	500
RF01	500	1000

==> tmp.00002.bed <==
RF01	100	600
RF01	600	1100

==> tmp.00003.bed <==
RF01	200	700
RF01	700	1200

==> tmp.00004.bed <==
RF01	300	800
RF01	800	1300

==> tmp.00005.bed <==
RF01	400	900
RF01	900	1400

$ cat tmp.manifest 
tmp.00001.bed	44
tmp.00002.bed	39
tmp.00003.bed	37
tmp.00004.bed	36
tmp.00005.bed	33
```

### Example

```bash
(...)
java -jar dist/bednonoverlappingset.jar -x 1 -R ref.fa -o "tmp.__SETID__.bed" -m tmp.manifest input.bed

cut -f 1 tmp.manifest | while read B
do
	${java_exe}   -Djava.io.tmpdir=.  -jar GenomeAnalysisTK.jar \
	   -T DepthOfCoverage -R "ref.fa" \
	   -o "SAMPLE"  -I input.bam  -L "${B}" --omitDepthOutputAtEachBase --omitLocusTable --omitPerSampleStats
 	   grep -v '^Target' "${sample}.sample_interval_summary" | awk -F '	' '{printf("%s\\t%s\\n",\$1,\$3);}' >> tmp.tsv	
	
	   rm "SAMPLE.sample_interval_summary"  "SAMPLE.sample_interval_statistics"
done
	
LC_ALL=C sort -t '	' -k1,1 tmp.tsv >> "SAMPLE.win.cov.tsv"

(...)
```

END_DOC
*/
@Program(
		name="bednonoverlappingset",
		description="Split a Bed file into non-overlapping data set.",
		keywords={"bed"},
		creationDate="20180607",
		modificationDate="20200408"
		)
public class BedNonOverlappingSet extends Launcher {
	private static final Logger LOG = Logger.build(BedNonOverlappingSet.class).make();

	@Parameter(names = { "-o", "--out" }, description = ArchiveFactory.OPT_DESC,required=true)
	private Path archivePath = null;
	@Parameter(names = { "-m", "--manifset" }, description = "Manifest file file containing the generated filenames/number of item.")
	private Path manifestFile = null;
	@Parameter(names = { "-R", "-r","--reference" }, description = INDEXED_FASTA_REFERENCE_DESCRIPTION +
						" If defined, will be used to sort the bed record on chrom/pos before writing the bed records.")
	private Path refFile = null;
	@Parameter(names = { "-x", "--extend" }, description = "Extend intervals by 'x' bases")
	private int extend = 0;
	@Parameter(names={"--compress"},description="Bgzip outut bed files")
	private boolean compressed = false;
	
	private class BedRecord
		implements Locatable
		{
		final BedLine bedline;
		BedRecord(final BedLine bedline) {
			this.bedline = bedline;
			}
		@Override
		public String getContig() {
			return this.bedline.getContig();
			}
		@Override
		public int getStart() {
			return this.bedline.getStart() - extend;
			}
		@Override
		public int getEnd() {
			return  this.bedline.getEnd() + extend;
			}
		}
	
	private final List<IntervalTreeMap<BedRecord>> bedsets= new ArrayList<>();
	
	private void scan(final BufferedReader br) throws IOException
		{
		final BedLineCodec codec = new BedLineCodec();
		br.lines().
			map(L->codec.decode(L)).
			filter(L->L!=null).
			map(B->new BedRecord(B)).
			forEach(B->{
				final Interval key = new Interval(B.getContig(), B.getStart(),B.getEnd());
				int y=0;
				for(y=0;y<  bedsets.size();++y)
					{
					final IntervalTreeMap<BedRecord> itm = bedsets.get(y);
					if(itm.containsOverlapping(key)) continue;
					itm.put(key, B);
					return;
					}
				final IntervalTreeMap<BedRecord> itm = new IntervalTreeMap<>();
				bedsets.add(itm);
				itm.put(key, B);
				});
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.extend<0) {
			LOG.info("extend cannot be negative.");
			return -1;
			}
		ArchiveFactory archiveFactory = null;
		PrintWriter manifest = null;
		BufferedReader br;
		try {
			final Comparator<String> cmpContig;
			if(this.refFile!=null)
				{
				final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(this.refFile);
				cmpContig = new ContigDictComparator(dict);
				}
			else
				{
				cmpContig = (A,B)->A.compareTo(B);
				}
			final Comparator<Locatable> cmpbed=(A,B) ->{
				int i= cmpContig.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(),B.getStart());
				if(i!=0) return i;
				i = Integer.compare(A.getEnd(),B.getEnd());
				return i;
				};
			
			if(this.manifestFile!=null) {
				manifest = IOUtils.openPathForPrintWriter(this.manifestFile);
				}
			else
				{
				manifest = new PrintWriter(new NullOuputStream());
				}
			
			if(args.isEmpty())
				{
				br = super.openBufferedReader(null);
				scan(br);
				br.close();
				}
			else
				{
				for(final String fname:args) {
					br = super.openBufferedReader(fname);
					scan(br);
					br.close();
					}
				}
			if(this.bedsets.isEmpty()) {
				LOG.warn("No set was created");
				}
			
			archiveFactory = ArchiveFactory.open(this.archivePath);
			if(compressed) archiveFactory.setCompressionLevel(0);
			
			for(int i=0;i< this.bedsets.size();i++){
				final BlockCompressedOutputStream bcos;
				final PrintWriter pw;
				String filename =  String.format("%05d", (i+1));
				filename+= FileExtensions.BED;
					
				if(this.compressed) {
					filename += ".gz";
					bcos = new BlockCompressedOutputStream(archiveFactory.openOuputStream(filename),(Path)null);
					pw = new PrintWriter(bcos);
					}
				else
					{
					bcos=null;
					pw = archiveFactory.openWriter(filename);
					}
				
				
				final IntervalTreeMap<BedRecord> itm = this.bedsets.get(i);
				final List<BedLine> list = 
						itm.values().
						stream().
						map(R->R.bedline).
						sorted(cmpbed).
						collect(Collectors.toList());
				LOG.info("saving "+filename+" "+list.size());
				
				
				for(final BedLine bl:list)
					{
					pw.println(bl.join());
					}
				if(bcos!=null) bcos.flush();
				pw.flush();
				pw.close();
				
				manifest.println(filename+"\t"+list.size());
				}
			archiveFactory.close();
			archiveFactory=null;
			manifest.flush();
			manifest.close();
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(archiveFactory);
			CloserUtil.close(manifest);
			}
		}
	
	public static void main(final String[] args) {
		new BedNonOverlappingSet().instanceMainWithExit(args);

	}

}
