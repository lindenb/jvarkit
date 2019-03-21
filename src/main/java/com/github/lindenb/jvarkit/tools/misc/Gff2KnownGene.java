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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.LineIterator;

/**
BEGIN_DOC

## Example

```
$  curl -s "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz" |\
	gunzip -c |\
	java -jar dist/gff2kg.jar
(...)
1826	ENST00000367917.3	chr1	+	162760522	162782607	162760590	162782210	8	162760522,162762448,162766374,162767591,162769532,162774056,162775183,162782087	162760625,162762652,162766467,162767706,162769727,162774113,162775282,162782607	gene_id=ENSG00000132196.9;transcript_id=ENST00000367917.3;gene_type=protein_coding;gene_status=KNOWN;gene_name=HSD17B7;transcript_type=protein_coding;transcript_name=HSD17B7-201;protein_id=ENSP00000356894.3;havana_gene=OTTHUMG00000034420.6;	ENST00000367917.3
(...)
```

In the UCSC (not the structure of konwGene, but we can validate intervals):

```
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -D hg19 -e 'select * from wgEncodeGencodeBasicV19 where name="ENST00000367917.3"' | cat
bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
1826	ENST00000367917.3	chr1	+	162760522	162782607	162760590	162782210	8	162760522,162762448,162766374,162767591,162769532,162774056,162775183,162782087,	162760625,162762652,162766467,162767706,162769727,162774113,162775282,162782607,	0	HSD17B7	cmpl	cmpl	0,2,2,2,0,0,0,0,
```

### From ensembl 

```
$	wget -O -  "ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz" |\
	gunzip -c |\
	java -jar dist/gff2kg.jar
```

## see also

  * Ensembl vs UCSC  [https://twitter.com/yokofakun/status/743751004785545218](https://twitter.com/yokofakun/status/743751004785545218)



END_DOC
 */
@Program(name="gff2kg",
		description="Convert GFF3/GTF format to UCSC knownGene format.",
		keywords={"gff",",gtf","knownGene","ucsc","convert"},
		biostars=276099
		)
public class Gff2KnownGene extends Launcher {
	private static final Logger LOG = Logger.build(Gff2KnownGene.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-bin","--bin"},description="Insert  UCSC 'bin' column as the first column.")
	private boolean writeBin = false;
	@Parameter(names={"-verbose","--verbose"},description="Be verbose, log messages")
	private boolean verbose = false;
	@Parameter(names={"-trid","--trid"},description="Transcript identifiers in the GTF/GFF (column N°3) used to identify a transcript."
			+ "Multiple separated by a semicolon ")
	private String transcriptIdentifiersStr = "transcript;mRNA;snRNA;tRNA;snoRNA";
	@ParametersDelegate
	private GTFCodec.FormatChooser formatChooser = new  GTFCodec.FormatChooser();
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	
	private static final String NO_TRANSCRIPT_NAME="\0\0NOTRANSCRIPT";
	private final static Pattern semicolon=Pattern.compile(";");
	private GTFCodec gtfCodec = null;

	/** a line in the gtf associated to a transcript */
	private class GffLine {
		final GTFLine delegate;
		private final Interval interval;
		private final String transcript;
		GffLine(final GTFLine delegate) {
			this.delegate=delegate;
			final Map<String, String> map = delegate.getAttributes();
			this.interval = new Interval(
					delegate.getContig(),
					delegate.getStart(),
					delegate.getEnd()
					);
			this.transcript = map.getOrDefault("transcript_id",
					 map.getOrDefault("Parent",NO_TRANSCRIPT_NAME)
					 );			
			}
		
		boolean hasTranscript() {
			return !(getTranscript()==null || NO_TRANSCRIPT_NAME.equals(getTranscript()));
			}
		
		
		public String getContig() {
			return delegate.getContig();
		}
		
		public String getTranscript() {
			return transcript;
		}
		
		@Override
		public String toString() {
			return delegate.getLine();
		}
	}
	
	private static class GffLineComparator implements Comparator<GffLine>{
		@Override
		public int compare(final GffLine o1,final GffLine o2)
			{
			int i= o1.getContig().compareTo(o2.getContig());
			if(i!=0) return i;
			i= o1.getTranscript().compareTo(o2.getTranscript());
			if(i!=0) return i;
			return o1.delegate.getLine().compareTo(o2.delegate.getLine());
			}
		}
	
	private class GffLineCodec extends AbstractDataCodec<GffLine> {
		@Override
		public GffLine decode(final DataInputStream dis) throws IOException {
			try {
				return new GffLine(Gff2KnownGene.this.gtfCodec.decode(readString(dis)));
			} catch (final EOFException e) {
			return null;
			}
		}
		@Override
		public void encode(final DataOutputStream dos, final GffLine line) throws IOException {
			writeString(dos, line.delegate.getLine());
		}
		@Override
		public AbstractDataCodec<GffLine> clone() {
			return new GffLineCodec();
			}
	}
	
	
    static int reg2bin(final int beg, int end)
    {
    	return GenomicIndexUtil.regionToBin(beg, end);
    }


	@Override
	public int doWork(final List<String> args) {
		
		this.gtfCodec = this.formatChooser.makeCodec();
		
		LineIterator in =null;
		EqualRangeIterator<GffLine> eq = null;
		CloseableIterator<GffLine> iter = null;
		SortingCollection<GffLine> sorting = null;
		PrintWriter pw = null;
		final Set<String> transcriptIdentifiersSet = 
				new HashSet<>(Arrays.asList(semicolon.split(this.transcriptIdentifiersStr)));
		
		final GffLineComparator comparator = new GffLineComparator();
		try {
			final String input = oneFileOrNull(args);
			in = (input==null?
					IOUtils.openStdinForLineIterator():
					IOUtils.openURIForLineIterator(input)
					);
						
			sorting = SortingCollection.newInstance(
					GffLine.class,
					new GffLineCodec(),
					comparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorting.setDestructiveIteration(true);
			int nRead=0;
			this.gtfCodec.readActualHeader(in);
			
			while(in.hasNext())
			{
				++nRead;
				final String line = in.next();
				if(line.isEmpty() || line.startsWith("#")) continue;
				final GTFLine delegate = this.gtfCodec.decode(line);
				if(delegate.getType().equals("gene")) {
					if(verbose) LOG.info("skipping "+line);
					continue;
				}
				final GffLine gffLine = new GffLine(delegate);
				if(!gffLine.hasTranscript()) {
					if(verbose) LOG.info("skipping "+line);
					continue;
				}
				sorting.add(gffLine);
				if(nRead %50000==0) LOG.info("Read "+nRead+" lines. Last: "+line);
			}
			
			sorting.doneAdding();
			LOG.info("sorting...."+nRead);
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			iter = sorting.iterator();
			eq = new EqualRangeIterator<>(iter,(o1,o2)->{
					final int i= o1.getContig().compareTo(o2.getContig());
					if(i!=0) return i;
					return o1.getTranscript().compareTo(o2.getTranscript());
					}
				);
			while(eq.hasNext()) {
				final List<GffLine> L = eq.next();
				final GffLine first = L.get(0);
				final String firstContig = first.getContig();
				final String firstTranscriptName  =  first.getTranscript();
				if(verbose) LOG.info("processing "+firstTranscriptName);
				final char strand  =  first.delegate.getStrand();
				if(!(strand=='+' || strand=='-')) {
					LOG.error("Bad strand in "+first.delegate.getLine());
					return -1;
				}
				final List<Interval> exons= new ArrayList<>();
				Interval mainTranscriptInterval = null;
				final List<Interval> cds= new ArrayList<>();
				
				for(final GffLine item: L) {
					if(!firstContig.equals(item.getContig())) {
						LOG.error("Conflict in contig!!");
						return -1;
					}
					if(!firstTranscriptName.equals(item.getTranscript())) {
						LOG.error("Conflict in name!! "+firstTranscriptName+":"+item.getTranscript());
						return -1;
					}
					
					if(item.delegate.getType().equals("gene")) {
						if(verbose) LOG.info("ignore line "+item);
						continue;
					}
					
					else if((transcriptIdentifiersSet.contains(item.delegate.getType()))) {
						if(mainTranscriptInterval!=null && !mainTranscriptInterval.equals(item.interval))
							{
							LOG.error("Transcript found twice for "+firstTranscriptName);
							return -1;
							}
						mainTranscriptInterval = item.interval;
						continue;
					} 
					
					
					else if(item.delegate.getType().equals("exon")) {
						exons.add( item.interval);
						continue;
					}
					
					else if(item.delegate.getType().equals("CDS")) {
						cds.add( item.interval);
						continue;
					}
					else //UTR , stop_codon, etc...
					{
						if(verbose) LOG.info("ignore line "+firstTranscriptName+":"+item.delegate.getType());
						continue;
					}
				}
			
			exons.sort((o1,o2)->o1.getStart()-o2.getStart());
			
			if(mainTranscriptInterval==null) {
				LOG.warn(
						"main transcript not found for "+
						firstTranscriptName+" "+first+
						" available feature type where:"
						+ L.stream().map(T->T.delegate.getType()).
							collect(Collectors.toSet()));
				continue;
				}
			
			if(this.writeBin) {
				pw.print(reg2bin(mainTranscriptInterval.getStart()-1, mainTranscriptInterval.getEnd()));
				pw.print("\t");
				}
			pw.print(firstTranscriptName);
			pw.print("\t");
			pw.print(firstContig);
			pw.print("\t");
			pw.print(strand);
			pw.print("\t");
			pw.print(mainTranscriptInterval.getStart()-1);
			pw.print("\t");
			pw.print(mainTranscriptInterval.getEnd());
			pw.print("\t");
			
			if(cds.isEmpty()) {
				pw.print(mainTranscriptInterval.getStart()-1);
				pw.print("\t");
				pw.print(mainTranscriptInterval.getStart()-1);
				pw.print("\t");
			}
			else {
				int minCds=cds.get(0).getStart();
				int maxCds=cds.get(0).getEnd();
				for(int i=1;i< cds.size();++i) {
					minCds = Math.min(cds.get(i).getStart(), minCds);
					maxCds = Math.max(cds.get(i).getEnd(), maxCds);
				}
			pw.print(minCds-1);
			pw.print("\t");
			pw.print(maxCds);
			pw.print("\t");
			}
			pw.print(exons.size());
			pw.print("\t");
			for(int i=0;i< exons.size();++i) {
				if(i>0) pw.print(",");
				pw.print(exons.get(i).getStart()-1);
			}
			pw.print("\t");
			for(int i=0;i< exons.size();++i) {
				if(i>0) pw.print(",");
				pw.print(exons.get(i).getEnd());
			}
			pw.print("\t");
			
			for(final Iterator<Map.Entry<String,String>> metainfoiter = first.delegate.getAttributeIterator();
					metainfoiter.hasNext();)  {
				final Map.Entry<String, String> entry = metainfoiter.next();
				if(entry==null) continue;
				final String s=entry.getKey();
				if(s.equals("gene_id") ||
				   s.equals("transcript_type") ||
				   s.equals("gene_name") ||
				   s.equals("gene_status") ||
				   s.equals("gene_type") ||
				   s.equals("transcript_id") ||
				   s.equals("havana_gene") ||
				   s.equals("havana_transcript") ||
				   s.equals("transcript_name") ||
				   s.equals("protein_id") ||
				   s.equals("ccdsid") ||
				   s.equals("Parent")
				   ) {
					pw.print(entry.getValue());
					pw.print(";");
				}
			}
			pw.print("\t");
			pw.print(firstTranscriptName);
			pw.println();
			}
			
			eq.close();
			iter.close();iter=null;
			sorting=null;
			pw.flush();pw.close();pw=null;
			
			LOG.info("done");
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
			}
		finally {
			CloserUtil.close(eq);
			CloserUtil.close(pw);
			CloserUtil.close(in);
			CloserUtil.close(iter);
			CloserUtil.close(sorting);
			}
		}
		
	public static void main(final String[] args) throws IOException
		{
		new Gff2KnownGene().instanceMainWithExit(args);
		}

	}
