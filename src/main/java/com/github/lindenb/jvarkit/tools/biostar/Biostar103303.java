package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class Biostar103303 extends AbstractCommandLineProgram
	{
	private IntervalTreeMap<List<GTFGene.Exon>> exonMap=new IntervalTreeMap<List<GTFGene.Exon>>();
	
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
			Exon exon=new Exon();
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
	private void readGTF(String uri ,SAMSequenceDictionary dict) throws IOException
		{
		int count_exons=0;
		Set<String> unknown=new HashSet<String>();
		info("Reading "+uri);
		Pattern tab=Pattern.compile("[\t]");
		Map<String,GTFGene> transcript2gene=new HashMap<String, GTFGene>();
		LineIterator iter=IOUtils.openURIForLineIterator(uri);
		while(iter.hasNext())
			{
			String line=iter.next();
			if(line.startsWith("#")) continue;
			String tokens[]=tab.split(line);
			if(tokens.length<9) continue;
			if(!tokens[2].equals("exon")) continue;
			if(dict.getSequence(tokens[0])==null)
				{
				if(!unknown.contains(tokens[0]))
					{
					warning("chromosome in "+line+" not in SAMSequenceDictionary ");
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
			GTFGene gene=transcript2gene.get(tokens[0]+" "+transcript_id);
			if(gene==null)
				{
				gene=new GTFGene();
				gene.transcript_id=transcript_id;
				gene.gene_id=gene_id;
				gene.gene_name=gene_name;
				gene.chrom=tokens[0];
				transcript2gene.put(tokens[0]+" "+transcript_id, gene);
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
				
				GTFGene.Exon exon=g.exons.get(i);
				exon.index=i;
				
				if(i>0)
					{
					GTFGene.Exon prev=g.exons.get(i-1);
					if(prev.end>=exon.start)
						{
						throw new IOException("exons "+(i)+" and "+(i+1)+" overlap in "+g);
						}
					}
				
				Interval interval=new Interval(
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
		info("End Reading "+uri+ " N="+count_exons);
		}
	private static Object notnull(Object o)
		{
		if(o==null || "".equals(o)) return ".";
		return o;
		}
	
	
	@Override
	public String getProgramDescription()
		{
		return "Calculate Percent Spliced In (PSI). See also https://www.biostars.org/p/103303/";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar103303";
		}

	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -g (uri) GTF file. required.");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		String gtfuri=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "g:"))!=-1)
			{
			switch(c)
				{
				case 'g': gtfuri= opt.getOptArg();break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		if(gtfuri==null)
			{
			error("GTF input misssing");
			return -1;
			}
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		PrintWriter out=new PrintWriter(System.out);
		try
			{
			SamFileReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			if(opt.getOptInd()==args.length)
				{
				info("Reading sfomr stdin");
				samReader=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				info("Reading from "+filename);
				samReader=SamFileReaderFactory.mewInstance().open(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			
			this.readGTF(gtfuri,samReader.getFileHeader().getSequenceDictionary());
			
			if(this.exonMap.isEmpty())
				{
				error("no exon found");
				return -1;
				}
			iter=samReader.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samReader.getFileHeader().getSequenceDictionary());
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				progress.watch(rec);
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				
				Cigar cigar=rec.getCigar();
				if(cigar==null) continue;
				
				
				
				
				for(List<GTFGene.Exon> L:this.exonMap.getOverlapping(new Interval(
						rec.getReferenceName(),
						rec.getAlignmentStart(),
						rec.getAlignmentEnd()
						)))
					{
					for(GTFGene.Exon exon:L)
					{
					
					boolean found_in_prev=false;
					boolean found_in_next=false;
					boolean found_in_curr=false;
					
					List<GTFGene.Exon> prev=exon.getPrev();
					List<GTFGene.Exon> next=exon.getNext();
					int refPos=rec.getAlignmentStart();
					for(CigarElement ce:cigar.getCigarElements())
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
			progress.finish();
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
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(out);
			}
		return 0;
		}

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException
		{
		new Biostar103303().instanceMainWithExit(args);
		}
		

	}
