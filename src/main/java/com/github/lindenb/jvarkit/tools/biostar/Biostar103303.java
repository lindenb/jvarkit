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
package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

/** 

BEGIN_DOC

##Example

```bash
$   curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqA549CellLongnonpolyaAlnRep1.bam" |\
  java -jar dist/biostar103303.jar -g "http://atgu.mgh.harvard.edu/plinkseq/dist/aux/gencodeBasicV11-hg19.gtf.gz"  > result.tsv
```

END_DOC

*/

@Program(name="biostar103303",
description="Calculate Percent Spliced In (PSI).",
keywords= {"sam","bam","psi"},
biostars=103303
)
public class Biostar103303 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar103303.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-g","--gtf"},description="GTF file",required=true)
	private String gtfuri = null;

	private static class GTFGene
		{
		 class Exon
			{
			String exon_id;
			int start;
			int end;
			int index;
			long count_prev_and_next=0L;
			long count_prev_and_curr=0L;
			long count_curr_and_next=0L;
			long count_curr_only=0L;
			long count_others=0L;
			
			
			GTFGene getGene()
				{
				return GTFGene.this;
				}
			
			List<Exon> getPrev()
				{
				if(index==0) return Collections.emptyList();
				return getGene().exons.subList(0,index);
				}
			List<Exon> getNext()
				{
				if(index+1>=getGene().exons.size()) return Collections.emptyList();
				return getGene().exons.subList(index+1,getGene().exons.size());
				}
			boolean contains(int pos)
				{
				return start<=pos && pos<=end;
				}
			@Override
			public String toString()
				{
				return ""+start+"-"+end;
				}
			}
		
		String chrom;
		String gene_name;
		String gene_id;
		String transcript_id;
		List<Exon> exons=new ArrayList<Exon>();
		
		Exon createExon(int start,int end)
			{
			final Exon exon=new Exon();
			exon.start=start;
			exon.end=end;
			this.exons.add(exon);
			return exon;
			}
		@Override
		public String toString()
			{
			return transcript_id+" "+exons;
			}
		
		}

		
	private final IntervalTreeMap<List<GTFGene.Exon>> exonMap=new IntervalTreeMap<List<GTFGene.Exon>>();
	
	private void readGTF(final String uri ,final SAMSequenceDictionary samDict) throws IOException
		{
		final ContigNameConverter ctgNameConverter = ContigNameConverter.fromOneDictionary(samDict);
		int count_exons=0;
		final Set<String> unknown=new HashSet<String>();
		LOG.info("Reading "+uri);
		final Pattern tab=Pattern.compile("[\t]");
		final Map<String,GTFGene> transcript2gene=new HashMap<String, GTFGene>();
		LineIterator iter=IOUtils.openURIForLineIterator(uri);
		while(iter.hasNext())
			{
			final String line=iter.next();
			if(line.startsWith("#")) continue;
			final String tokens[]=tab.split(line);
			if(tokens.length<9) continue;
			if(!tokens[2].equals("exon")) continue;
			final String normContig = ctgNameConverter.apply(tokens[0]);
			if(StringUtil.isBlank(normContig))
				{
				if(!unknown.contains(tokens[0]))
					{
					LOG.warn(JvarkitException.ContigNotFoundInDictionary.getMessage(tokens[0],samDict));
					unknown.add(tokens[0]);
					}
				continue;
				}
			
			String transcript_id=null,gene_id=null,gene_name=null,exon_id=null;
			StreamTokenizer st=new StreamTokenizer(new StringReader(tokens[8]));
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
					}
				else 
					{
					if(key.equals("transcript_id"))
						{
						transcript_id=s;
						}
					else if(key.equals("gene_id"))
						{
						gene_id=s;
						}
					else if(key.equals("gene_name"))
						{
						gene_name=s;
						}
					else if(key.equals("exon_id"))
						{
						exon_id=s;
						}
					key=null;
					}
				}
			if(transcript_id==null || transcript_id.isEmpty()) continue;
			GTFGene gene=transcript2gene.get(normContig+" "+transcript_id);
			if(gene==null)
				{
				gene=new GTFGene();
				gene.transcript_id = transcript_id;
				gene.gene_id = gene_id;
				gene.gene_name = gene_name;
				gene.chrom = normContig;
				transcript2gene.put(normContig+" "+transcript_id, gene);
				}
			GTFGene.Exon exon=gene.createExon(
					Integer.parseInt(tokens[3]),
					Integer.parseInt(tokens[4])
					);
			exon.exon_id=exon_id;
			}
		CloserUtil.close(iter);
		
		for(GTFGene g: transcript2gene.values())
			{
			Collections.sort(g.exons,new Comparator<GTFGene.Exon>()
				{
				@Override
				public int compare(GTFGene.Exon o1, GTFGene.Exon o2)
					{
					return o1.start-o2.start;
					}
				});
			for(int i=0;i< g.exons.size();++i)
				{
				
				final GTFGene.Exon exon=g.exons.get(i);
				exon.index=i;
				
				if(i>0)
					{
					final GTFGene.Exon prev=g.exons.get(i-1);
					if(prev.end>=exon.start)
						{
						throw new IOException("exons "+(i)+" and "+(i+1)+" overlap in "+g);
						}
					}
				
				final Interval interval=new Interval(
						g.chrom,
						exon.start,
						exon.end
						);
				List<GTFGene.Exon> L=exonMap.get(interval);
				if(L==null)
					{
					L=new ArrayList<GTFGene.Exon>(1);
					exonMap.put(interval,L);
					}
				L.add(exon);
				++count_exons;
				}
			}
		LOG.info("End Reading "+uri+ " N="+count_exons);
		}
	
	private static Object notnull(final Object o)
		{
		if(o==null || "".equals(o)) return ".";
		return o;
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.gtfuri==null)
			{
			LOG.error("GTF input misssing");
			return -1;
			}
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		PrintWriter out=null;
		
		try
			{
			out = super.openFileOrStdoutAsPrintWriter(outputFile);
			final SamReaderFactory srf= super.createSamReaderFactory();
			if(args.isEmpty())
				{
				LOG.info("Reading sfomr stdin");
				samReader=srf.open(SamInputResource.of(stdin()));
				}
			else if(args.size()==1)
				{
				final File filename=new File(args.get(0));
				LOG.info("Reading from "+filename);
				samReader=srf.open(filename);
				}
			else
				{
				LOG.error("Illegal number of arguments.");
				return -1;
				}
			
			this.readGTF(gtfuri,SequenceDictionaryUtils.extractRequired(samReader.getFileHeader()));
			
			if(this.exonMap.isEmpty())
				{
				LOG.error("no exon found");
				return -1;
				}
			iter=samReader.iterator();
			
			
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.newInstance().dictionary(samReader.getFileHeader()).logger(LOG).build();
			
			while(iter.hasNext())
				{
				final SAMRecord rec= progress.apply(iter.next());
				
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				
				final Cigar cigar=rec.getCigar();
				if(cigar==null) continue;
				
				
				
				
				for(final List<GTFGene.Exon> L:this.exonMap.getOverlapping(new Interval(
						rec.getReferenceName(),
						rec.getAlignmentStart(),
						rec.getAlignmentEnd()
						)))
					{
					for(final GTFGene.Exon exon:L)
						{
						boolean found_in_prev=false;
						boolean found_in_next=false;
						boolean found_in_curr=false;
						
						final List<GTFGene.Exon> prev=exon.getPrev();
						final List<GTFGene.Exon> next=exon.getNext();
						int refPos=rec.getAlignmentStart();
						for(final CigarElement ce:cigar.getCigarElements())
							{
							switch(ce.getOperator())
								{
								case M:case X:case EQ:
									{
									for(int i=0;i< ce.getLength();++i)
										{
										for(GTFGene.Exon ex2:prev)
											{
											if(ex2.contains(refPos))
												{
												found_in_prev=true;
												}
											}
										for(GTFGene.Exon ex2:next)
											{
											if(ex2.contains(refPos))
												{
												found_in_next=true;
												}
											}
										if(exon.contains(refPos))
											{
											found_in_curr=true;
											}
										refPos++;
										}
									break;
									}
								default:
									{
									if(ce.getOperator().consumesReferenceBases())
										{
										refPos+=ce.getLength();
										}
									break;
									}
								}
							}
						if(found_in_prev && found_in_next && !found_in_curr)
							{
							exon.count_prev_and_next++;
							}
						else if(found_in_prev && !found_in_next && found_in_curr)
							{
							exon.count_prev_and_curr++;
							}
						else if(!found_in_prev && found_in_next && found_in_curr)
							{
							exon.count_curr_and_next++;
							}
						else if(!found_in_prev && !found_in_next && found_in_curr)
							{
							exon.count_curr_only++;
							}
						else if(!found_in_curr && !found_in_next &&!found_in_prev)
							{
							//??
							}
						else 
							{
							exon.count_others++;
							}
						}
					}
				}
			progress.close();
			out.print("#chrom");
			out.print("\t");
			out.print("exon.start");
			out.print("\t");
			out.print("exon.end");
			out.print("\t");
			out.print("exon.exon_id");
			out.print("\t");
			out.print("exon.index5_3");
			out.print("\t");
			out.print("transcript_id");
			out.print("\t");
			out.print("gene_name");
			out.print("\t");
			out.print("gene_id");
			out.print("\t");
			out.print("exon.count_prev_and_next");
			out.print("\t");
			out.print("exon.count_prev_and_curr");
			out.print("\t");
			out.print("exon.count_curr_and_next");
			out.print("\t");
			out.print("exon.count_curr_only");
			out.print("\t");
			out.print("exon.count_others");
			out.println();

			for(List<GTFGene.Exon> L:exonMap.values()) {
			for(GTFGene.Exon exon:L)
				{
				out.print(exon.getGene().chrom);
				out.print("\t");
				out.print(exon.start);
				out.print("\t");
				out.print(exon.end);
				out.print("\t");
				out.print(notnull(exon.exon_id));
				out.print("\t");
				out.print(""+(exon.index+1)+"/"+exon.getGene().exons.size());
				out.print("\t");
				out.print(exon.getGene().transcript_id);
				out.print("\t");
				out.print(notnull(exon.getGene().gene_name));
				out.print("\t");
				out.print(notnull(exon.getGene().gene_id));
				out.print("\t");
				out.print(exon.count_prev_and_next);
				out.print("\t");
				out.print(exon.count_prev_and_curr);
				out.print("\t");
				out.print(exon.count_curr_and_next);
				out.print("\t");
				out.print(exon.count_curr_only);
				out.print("\t");
				out.print(exon.count_others);
				out.println();
				}
				}
			out.flush();
			return 0;
			}
		catch(final Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(out);
			}
		}
		
	public static void main(final String[] args) throws IOException
		{
		new Biostar103303().instanceMainWithExit(args);
		}
		

	}
