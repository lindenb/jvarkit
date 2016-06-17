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

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

public class Gff2KnownGene extends AbstractGff2KnownGene {
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Gff2KnownGene.class);
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
        --end;

        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return  ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return  ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return  ((1<<3)-1)/7 + (beg>>26);
        return 0;
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

	
	@SuppressWarnings("resource")
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		BufferedReader in =null;
		EqualRangeIterator<GffLine> eq = null;
		CloseableIterator<GffLine> iter = null;
		SortingCollection<GffLine> sorting = null;
		PrintWriter pw = null;
		final GffLineComparator comparator = new GffLineComparator();
		try {
			if(inputName==null)  {
				in = IOUtils.openStreamForBufferedReader(stdin());
			} else
				{
				in = IOUtils.openURIForBufferedReading(inputName);
				}
			
			sorting = SortingCollection.newInstance(
					GffLine.class,
					new GffLineCodec(),
					comparator,
					super.maxRecordsInRam
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
			pw = super.openFileOrStdoutAsPrintWriter();
			
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
					throw new IOException("Bad strand in "+first.line);
				}
				final List<Interval> exons= new ArrayList<>();
				Interval tx= null;
				final List<Interval> cds= new ArrayList<>();
				
				for(final GffLine item: L) {
					item.decode();//force decode components
					tokens = tab.split(item.line);
					if(!contig.equals(item.getContig())) {
						throw new IOException("Conflict in contig!!");
					}
					if(!name.equals(item.getTranscript())) {
						throw new IOException("Conflict in name!!");
					}
					if(tokens[2].equals("gene")) {
						continue;
					}
					
					else if(tokens[2].equals("transcript")) {
						if(tx!=null)throw new IOException("Transcript found twice for "+name);
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
			
			if(super.writeBin) {
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
			return wrapException(e);
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
