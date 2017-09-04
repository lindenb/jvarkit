/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

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
@Program(name="gff2knowngene",
		description="Convert GFF3 format to UCSC knownGene format.",
		keywords={"gff","knownGene","ucsc","convert"}
		)
public class Gff2KnownGene extends Launcher {
	private static final Logger LOG = Logger.build(Gff2KnownGene.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-bin","--bin"},description="Preppend with 'bin' column")
	private boolean writeBin = false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	
	private static final String NO_TRANSCRIPT_NAME="\0\0NOTRANSCRIPT";
	
	private static class GffLine {
		final String line;
		private Interval interval = null;
		private String transcript = null;
		GffLine(final String line) {
			this.line=line;
		}
		
		boolean hasTranscript()
		{
			return !(getTranscript()==null || NO_TRANSCRIPT_NAME.equals(getTranscript()));
		}
		
		
		private boolean decode() {
			final Pattern tab=Pattern.compile("[\t]");
			
			final String tokens[]=tab.split(line);
			if(tokens.length<9) throw new RuntimeIOException("not enough columns in gff3 line "+line);
			this.interval = new Interval(tokens[0],
					Integer.parseInt(tokens[3]),
					Integer.parseInt(tokens[4])
					);
			final Pattern semicolon=Pattern.compile(";");
			final String metainfo[] = semicolon.split(tokens[8]);
			for(final String meta:metainfo) {
				 Map.Entry<String,String> kv = split(meta);
				
				if(kv.getKey().equals("transcript_id")) {
					this.transcript = kv.getValue();
					break;
				}
			}
			if(transcript==null) {
				//LOG.warn("no transcript_id found in gff3 line "+line+" meta:\n\t"+String.join("\n\t",metainfo));
				this.transcript=NO_TRANSCRIPT_NAME;
				return false;
			}
			return true;
		}
		
		public String getContig() {
			if(interval==null) decode();
			return interval.getContig();
		}
		
		public String getTranscript() {
			if(transcript==null) decode();
			return transcript;
		}
		
		@Override
		public String toString() {
			return line;
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
			return o1.line.compareTo(o2.line);
			}
	}
	
	private static class GffLineCodec extends AbstractDataCodec<GffLine> {
		@Override
		public GffLine decode(final DataInputStream dis) throws IOException {
			try {
				return new GffLine(readString(dis));
			} catch (Exception e) {
			return null;
			}
		}
		@Override
		public void encode(DataOutputStream dos, final GffLine line) throws IOException {
			writeString(dos, line.line);
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

	private  static Map.Entry<String,String> split(String s){
		s=s.trim();
		int x=0;
		final StringBuilder key = new StringBuilder();
		while(x<s.length() && (s.charAt(x)=='_' || Character.isLetter(s.charAt(x)) || Character.isDigit(s.charAt(x))))
			{
			key.append(s.charAt(x));
			x++;
			}
		while(x<s.length() && Character.isWhitespace(s.charAt(x))) {
			if(s.charAt(x)=='=') {
				++x;
				break;
			}
			++x;
			}
		
		String value = s.substring(x).trim();
		if(value.startsWith("\"") && value.endsWith("\"")) {
			final StringBuilder sb=new StringBuilder(value.length());
			
			for(int i=1;i<value.length()-1;++i){
				if(value.charAt(i)=='\\') {
					if(i+1<value.length()-1)
						{
						i++;
						switch(value.charAt(i)) {
						case '"': sb.append("\"");break;
						case '\'': sb.append("\'");break;
						case 't': sb.append("\t");break;
						case 'n': sb.append("\n");break;
						default:break;
						}
						}
					continue;
				} else
				{
					 sb.append(value.charAt(i));
				}
			}
			
			value=sb.toString();
		}
		return new AbstractMap.SimpleEntry<String,String>(key.toString(),value);
	}

	/*
	private abstract class GffCodec
		{
		Map<String,Long> seqdict=new LinkedHashMap<>();
		Set<String> att_keys=new HashSet<>();
		Set<String> sources=new HashSet<>();
		Set<String> types=new HashSet<>();
		protected Pattern tab=Pattern.compile("[\t]");
		protected Pattern semicolon=Pattern.compile("[;]");
		
	
	private class DefaultCodec extends GffCodec
		{
		@Override
		void writeAttributes(XMLStreamWriter w,String attrString) throws XMLStreamException,IOException
			{
			StreamTokenizer st=new StreamTokenizer(new StringReader(attrString));
			st.wordChars('_', '_');
			String key=null;
			while(st.nextToken() != StreamTokenizer.TT_EOF)
				{
				String s=null;
				switch(st.ttype)
					{
					case StreamTokenizer.TT_NUMBER: s=String.valueOf(st.nval);break;
					case '"': case '\'' : case StreamTokenizer.TT_WORD: s=st.sval;break;
					case ';':break;
					default:break;
					}
				if(s==null) continue;
				if(key==null)
					{
					key=s;
					this.att_keys.add(key);
					}
				else 
					{
					w.writeStartElement(key);
					w.writeCharacters(s);
					w.writeEndElement();
					key=null;
					}
				}
			
			}
		}
		*/

	@Override
	public int doWork(List<String> args) {
		BufferedReader in =null;
		EqualRangeIterator<GffLine> eq = null;
		CloseableIterator<GffLine> iter = null;
		SortingCollection<GffLine> sorting = null;
		PrintWriter pw = null;
		final GffLineComparator comparator = new GffLineComparator();
		try {
			in = super.openBufferedReader(oneFileOrNull(args));
			
			
			sorting = SortingCollection.newInstance(
					GffLine.class,
					new GffLineCodec(),
					comparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpDirectories()
					);
			sorting.setDestructiveIteration(true);
			String line;
			int nRead=0;
			while((line=in.readLine())!=null)
			{
				++nRead;
				if(line.isEmpty() || line.startsWith("#")) continue;
				final GffLine gffLine = new GffLine(line);
				if(!gffLine.hasTranscript()) continue;
				sorting.add(gffLine);
				if(nRead %50000==0) LOG.info("Read "+nRead+" lines. Last: "+line);

			}
			sorting.doneAdding();
			LOG.info("sorting....");
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			iter = sorting.iterator();
			eq = new EqualRangeIterator<>(iter,new Comparator<GffLine>() {
				@Override
				public int compare(final GffLine o1,final GffLine o2)
					{
					int i= o1.getContig().compareTo(o2.getContig());
					if(i!=0) return i;
					i= o1.getTranscript().compareTo(o2.getTranscript());
					return i;
					}
				});
			final Pattern tab = Pattern.compile("[\t]");
			final Pattern semicolon=Pattern.compile(";");
			while(eq.hasNext()) {
				final List<GffLine> L = eq.next();
				final GffLine first = L.get(0);
				final String contig = first.getContig();
				final String name  =  first.getTranscript();
				String tokens[] = tab.split(first.line);
				final String strand  =  tokens[6];
				final String meta  =  tokens[8];
				if(!(strand.equals("+") || strand.equals("-"))) {
					LOG.error("Bad strand in "+first.line);
					return -1;
				}
				final List<Interval> exons= new ArrayList<>();
				Interval tx= null;
				final List<Interval> cds= new ArrayList<>();
				
				for(final GffLine item: L) {
					item.decode();//force decode components
					tokens = tab.split(item.line);
					if(!contig.equals(item.getContig())) {
						LOG.error("Conflict in contig!!");
						return -1;
					}
					if(!name.equals(item.getTranscript())) {
						LOG.error("Conflict in name!!");
						return -1;
					}
					if(tokens[2].equals("gene")) {
						continue;
					}
					
					else if(tokens[2].equals("transcript")) {
						if(tx!=null)
							{
							LOG.error("Transcript found twice for "+name);
							return -1;
							}
						tx = item.interval;
						continue;
					}
					
					else if(tokens[2].equals("exon")) {
						exons.add( item.interval);
						continue;
					}
					
					else if(tokens[2].equals("CDS")) {
						cds.add( item.interval);
						continue;
					}
				}
			
			exons.sort(new Comparator<Interval>() {
				@Override
				public int compare(Interval o1, Interval o2) {
					return o1.getStart()-o2.getStart();
				}
			});
			
			if(tx==null) {
				LOG.warn("tx is null for "+name+" "+first);
				continue;
				}
			
			if(this.writeBin) {
				pw.print(reg2bin(tx.getStart()-1, tx.getEnd()));
				pw.print("\t");
				}
			pw.print(name);
			pw.print("\t");
			pw.print(contig);
			pw.print("\t");
			pw.print(strand);
			pw.print("\t");
			pw.print(tx.getStart()-1);
			pw.print("\t");
			pw.print(tx.getEnd());
			pw.print("\t");
			
			if(cds.isEmpty()) {
				pw.print(tx.getStart()-1);
				pw.print("\t");
				pw.print(tx.getStart()-1);
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
			
			for(final String metainfo:semicolon.split(meta)) {
				final Map.Entry<String, String> entry = split(metainfo);
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
				   s.equals("ccdsid")
				   ) {
					pw.print(entry.getValue());
					pw.print(";");
				}
			}
			pw.print("\t");
			pw.print(name);
			pw.println();
			}
			
			eq.close();
			iter.close();iter=null;
			sorting.cleanup();sorting=null;
			pw.flush();pw.close();pw=null;
			
			LOG.info("done");
			return RETURN_OK;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally {
			CloserUtil.close(eq);
			CloserUtil.close(pw);
			CloserUtil.close(in);
			CloserUtil.close(iter);
			if(sorting!=null) sorting.cleanup();
			CloserUtil.close(sorting);
			}
		}
		
	public static void main(String[] args) throws IOException
		{
		new Gff2KnownGene().instanceMainWithExit(args);
		}


	

	}
