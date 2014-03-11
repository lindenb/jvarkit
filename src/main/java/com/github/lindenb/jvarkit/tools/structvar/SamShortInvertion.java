package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
//import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;

public class SamShortInvertion extends AbstractCommandLineProgram
	{
	private int max_size_inversion=2000;
	private int min_coverage=10;

	
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -R (ref) "+getMessageBundle("reference.faidx"));
		out.println(" -m (int) max size of inversion . Default: "+max_size_inversion);
		out.println(" -c (int) min coverage . Default: "+min_coverage);
		super.printOptions(out);
		}
	
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:m:c:"))!=-1)
			{
			switch(c)
				{
				case 'm':max_size_inversion= Integer.parseInt(opt.getOptArg());break;
				case 'c':min_coverage= Integer.parseInt(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
	
		SAMFileReader r=null;
		PrintStream out=System.out;
		//SAMFileWriter w=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				r=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				r=new SAMFileReader(new File(filename));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			r.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=r.getFileHeader();
			OtherCanonicalAlignFactory xpalignFactory=new OtherCanonicalAlignFactory(header);
			int prev_tid=-1;
			short coverage[]=null;
			short max_coverage=0;
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			SAMRecordIterator it= r.iterator();
			for(;;)
				{
				SAMRecord rec=null;
				if(it.hasNext())
					{
					rec=it.next();
					progress.watch(rec);
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getDuplicateReadFlag()) continue;
					}
				if(rec==null || prev_tid==-1 || (prev_tid!=-1 && prev_tid!=rec.getReferenceIndex()))
					{
					if(coverage!=null)
						{

						while(max_coverage>=Math.max(1,this.min_coverage))
							{
							info("Scanning "+header.getSequence(prev_tid).getSequenceName()+" for cov:"+max_coverage);

							int chromStart0=0;
							while(chromStart0 < coverage.length)
								{
								if(coverage[chromStart0]==max_coverage)
									{
									coverage[chromStart0]=0;//reset that pos
									int chromEnd0=chromStart0+1;
									while(chromEnd0 < coverage.length && coverage[chromEnd0]==max_coverage)
										{
										coverage[chromEnd0]=0;//reset that pos
										++chromEnd0;
										}
									out.print(header.getSequence(prev_tid).getSequenceName());
									out.print('\t');
									out.print(chromStart0);
									out.print('\t');
									out.print(chromEnd0);
									out.print('\t');
									out.print(max_coverage);
									out.println();
									
									//reset 3'
									for(int x=chromEnd0;x<coverage.length && coverage[x]>0;++x)
										{
										coverage[x]=0;
										}
									//reset 5'
									for(int x=chromStart0-1;x>=0 && coverage[x]>0;--x)
										{
										coverage[x]=0;
										}
									chromStart0=chromEnd0;
									}
								else
									{
									++chromStart0;
									}
								}
							max_coverage--;
							}
						coverage=null;
						}
					
					if(rec==null) break;
					prev_tid=rec.getReferenceIndex();
					info("Alloc sizeof(short)*"+header.getSequence(prev_tid).getSequenceLength());
					coverage=new short[header.getSequence(prev_tid).getSequenceLength()];
					Arrays.fill(coverage,(short)0);
					max_coverage=0;
					}
				
				List<OtherCanonicalAlign> saList=xpalignFactory.getXPAligns(rec);
				if(saList.isEmpty()) continue;
				for(OtherCanonicalAlign xp:saList)
					{
					if(!xp.getChrom().equals(rec.getReferenceName())) continue;
					
					if(!rec.getReadNegativeStrandFlag()) //read is plus
						{
						if(xp.getStrand()=='+')
							{
							//info(xp+" vs strand "+rec.getReadNegativeStrandFlag());
							continue;//ignore both same strand
							}
						}
					else //read.strand=='-'
						{
						if(xp.getStrand()=='-')
							{
							//info(xp+" vs strand "+rec.getReadNegativeStrandFlag());
							continue;//ignore both same strand
							}
						}
					if(Math.abs(rec.getUnclippedStart() - xp.getPos()) > max_size_inversion) 
						{
						//info(xp+" vs pos "+rec.getUnclippedStart());
						continue;
						}
					int chromStart=Math.min(rec.getUnclippedStart(), xp.getPos())-1;
					int chromEnd=Math.max(rec.getUnclippedEnd(), xp.getPos())-1;
					for(int x=chromStart;x<=chromEnd && x< coverage.length;++x )
						{
						if(coverage[x]<Short.MAX_VALUE) coverage[x]++;
						if(max_coverage< coverage[x])
							{
							info("Max coverage "+max_coverage);
							max_coverage=coverage[x];
							}
						}
					}
				}
			
			it.close();
			progress.finish();
			return 0;
			}

		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			//CloserUtil.close(w);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
