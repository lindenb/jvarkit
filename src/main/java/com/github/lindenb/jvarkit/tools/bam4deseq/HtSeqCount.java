package com.github.lindenb.jvarkit.tools.bam4deseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.bin.SamSequenceRecordBinMap;

public class HtSeqCount extends AbstractCommandLineProgram
	{
	/** features we want in the GTF */
	private Set<String> features=null;
	/** first Dict found, to print chrom names and compare with others */
	private SAMSequenceDictionary firstDict=null;
	/** mapping position to transcript */
	private SamSequenceRecordBinMap<Transcript> pos2transcript=null;
	/** all filenames */
	private List<String> filenames=new ArrayList<String>();
	/** user GTF file */
	private File gtfFile=null;
	
	/* current file scanner in filenames */
	private int file_index=0;
	/* user feature  in the GTF */
	private String transcript_id="transcript_id";
	/* all transcripts to print at the end */
	private HashMap<String, Transcript> name2transcript=new HashMap<String, Transcript>();
	/* first ouput line is header */
	private boolean print_header=false;
	/* remove lines with 0 coverage */
	private boolean removeZero=false;
	/* print chrom/start-end  on the left*/
	private boolean print_chromstartend=false;

	
	private static class Transcript
		{
		String name;
		int count[];
		boolean bad_flag=false;//multiple chroms
		int tid;
		int start=Integer.MAX_VALUE;
		int end=0;
		}
	
	@Override
	public String getProgramDescription() {
		return "java version of htseqcount";
		}
	
	private void parseGTF(BufferedReader in,SAMSequenceDictionary dict) throws IOException
		{
		this.pos2transcript=new SamSequenceRecordBinMap<Transcript>(dict);
		final Pattern tab=Pattern.compile("[\t]");
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line);
			if(tokens.length<=8)
				{
				throw new IOException("Bad GTF line " +line);
				}
			if(!this.features.contains(tokens[2])) continue;
			String chrom=tokens[0];
			int tid=dict.getSequenceIndex(chrom);
			if(tid<0)
				{
				warning("Unknown chromosome "+chrom+" Ignoring "+line);
				continue;
				}
			String transcript_name=null;
			String info=tokens[8];
			StreamTokenizer scanner=new StreamTokenizer(new StringReader(info));
			scanner.wordChars('_', '_');
			scanner.wordChars('0', '9');
			while(scanner.nextToken()!=StreamTokenizer.TT_EOF)
				{
				String key=scanner.sval;
				//info("K "+key);
				if(scanner.nextToken()==StreamTokenizer.TT_EOF)
					{
					throw new IOException("Bad gtf in "+line+" after "+key);
					}
			
				String value=scanner.sval;
				//info("v "+value);
				if(key.equals(this.transcript_id))
					{
					transcript_name=value;
					break;
					}
				if(scanner.nextToken()==StreamTokenizer.TT_EOF) break;
				if(scanner.ttype!=';')
					{
					throw new IOException("Bad gtf in "+line+" expect ; ");
					}
				}
			if(transcript_name==null)
				{
				info("No "+transcript_id+" in "+line);
				continue;
				}
			
			
			Transcript transcript=name2transcript.get(transcript_name);
			if(transcript==null)
				{
				transcript=new Transcript();
				transcript.tid=tid;
				transcript.name=transcript_name;
				transcript.count=new int[filenames.size()];
				Arrays.fill(transcript.count,0);
				name2transcript.put(transcript_name,transcript);
				if(name2transcript.size()%1000==0)
					{
					info("Transcripts: "+name2transcript.size()+" "+transcript_name);
					}
				}
			else
				{
				if(transcript.tid!=tid)
					{
					warning("Multiple chromosomes for "+transcript);
					transcript.bad_flag=true;
					continue;
					}
				}
			
			int start1=Integer.parseInt(tokens[3]);
			int end1=Integer.parseInt(tokens[4]);
			transcript.start=Math.min(transcript.start,start1);
			transcript.end=Math.max(transcript.end,end1);
			
			
			this.pos2transcript.put(tid, start1-1, end1,transcript);
			
			/**
			boolean ok=false;
			Iterator<Transcript>x=this.pos2transcript.overlapping(tid, start1-1, end1);
			while(x.hasNext()) if(transcript==x.next()) {ok=true;}
			if(!ok) throw new IllegalStateException("boum");*/
			}
		info("Done Reading transcripts:"+name2transcript.size());
		}
	
	private void touch(int tid,int start1,int end1)
		{
		Set<String> seen=null;
		Iterator<Transcript> iter=this.pos2transcript.overlapping(tid, start1-1, end1);
		while(iter.hasNext())
			{
			Transcript transcript=iter.next();
			if(seen==null) seen=new HashSet<String>();
			if(!seen.add(transcript.name)) continue;
			transcript.count[this.file_index]++;
			}
		}
	
	private void run(SAMFileReader sfr) throws IOException
		{
		sfr.setValidationStringency(ValidationStringency.LENIENT);
		SAMFileHeader header=sfr.getFileHeader();
		
		if(this.file_index==0)
			{
			firstDict=header.getSequenceDictionary();
			
			info("Reading "+this.gtfFile);
			BufferedReader gtfIn=IOUtils.openFileForBufferedReading(this.gtfFile);
			parseGTF(gtfIn,header.getSequenceDictionary());
			gtfIn.close();
			}
		else
			{
			if(!SequenceUtil.areSequenceDictionariesEqual(firstDict, header.getSequenceDictionary()))
				{
				throw new IOException("not same sequence dictionaires between "+
					filenames.get(0)+" && "+
					filenames.get(this.file_index)
					);
				}
			}
		long nReads=0L;
		SAMRecordIterator iter=sfr.iterator();
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(rec.getNotPrimaryAlignmentFlag()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(nReads++%1E6==0)
				{
				info("Read "+nReads+" in "+this.filenames.get(this.file_index));
				}
			touch(rec.getReferenceIndex(),
					rec.getAlignmentStart(),
					rec.getAlignmentEnd()
					);
			}
		iter.close();
		}
	

	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -g (file) path to GTF file. REQUIRED");
		out.println(" -z remove 0 covered file");
		out.println(" -c print chrom/start/end");
		out.println(" -H print header");
		out.println(" -F <string> add this feature. default are : CDS and exon");
		out.println(" -T <type> key for clustering info. Default is "+this.transcript_id);
		
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, "hvL:g:zcHF:T:"))!=-1)
			{
			switch(c)
				{
				
				case 'F':
					{
					if(features==null) features=new HashSet<String>();
					features.add(opt.getOptArg());
					break;
					}
				case 'T':
					{
					this.transcript_id=opt.getOptArg();
					break;
					}
				case 'g': this.gtfFile=new File(opt.getOptArg());break;
				case 'z': this.removeZero=true; break;
				case 'c': this.print_chromstartend=true; break;
				case 'H': this.print_header=true; break;
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(opt.getOptArg()));break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+opt.getOptOpt());return -1;
				}
			}
		if(features==null)
			{
			features=new HashSet<String>();
			for(String F:new String[]{"CDS","exon"})
				{
				info("Adding "+F+" as default feature.");
				features.add(F);
				}
			}
		if( this.gtfFile==null)
			{
			error("undefined GTF file");
			return -1;
			}
		
		
		SAMFileReader sfr=null;
		try
			{
			this.file_index=0;
			
			if(opt.getOptInd()==args.length)
				{
				info("Opening stdin");
				filenames.add("stdin");
				sfr=new SAMFileReader(System.in);
				run(sfr);
				sfr.close();
				}
			else 
				{
				for(int optind=opt.getOptInd();optind< args.length;++optind)
					{
					String filename=args[optind];
					filenames.add(filename);
					}
				
				for(String filename:filenames)
					{
					info("Opening "+filename);
					sfr=new SAMFileReader(new File(filename));
					run(sfr);
					sfr.close();
					this.file_index++;
					}
				}
			
			if(this.print_header)
				{
				System.out.print("#");
				if(this.print_chromstartend)
					{
					System.out.print("chrom\tstart\tend\t");
					}
				System.out.print("Transcript");
				for(String filename:filenames)
					{
					System.out.print("\t"+filename);
					}
				System.out.println();
				}
			
			List<Transcript> ordered=new ArrayList<Transcript>(this.name2transcript.values());
			Collections.sort(ordered,new java.util.Comparator<Transcript>()
					{
					@Override
					public int compare(Transcript a,Transcript b)
						{
						int i= a.tid-b.tid;
						if(i!=0) return i;
						i= a.start-b.start;
						if(i!=0) return i;
						i= a.end-b.end;
						if(i!=0) return i;
						return a.name.compareTo(b.name);
						}
					});
			
			for(Transcript tr:ordered)
				{
				if(tr.bad_flag) continue;
				
				if(this.removeZero)
					{
					//remove unmapped transcripts
					int x=0;
					for(x=0;x< tr.count.length;++x)
						if(tr.count[x]!=0) break;
					if(x==tr.count.length) continue;
					}
				if(this.print_chromstartend)
					{
					System.out.print(firstDict.getSequence(tr.tid).getSequenceName());
					System.out.print("\t");
					System.out.print(tr.start);
					System.out.print("\t");
					System.out.print(tr.end);
					System.out.print("\t");
					}
				System.out.print(tr.name);
				for(int C:tr.count)
					{
					System.out.print("\t");
					System.out.print(C);
					}
				System.out.println();
				}
		
		
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			}
		}
	public static void main(String[] args) {
		new HtSeqCount().instanceMainWithExit(args);
		}
	}
