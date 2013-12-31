package com.github.lindenb.jvarkit.tools.biostar;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import net.sf.picard.PicardException;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class Biostar84452 extends AbstractCommandLineProgram
	{
	private Biostar84452()
		{
		}
	
	@Override
	public String getProgramDescription()
		{
		return "remove clipped bases from BAM. See: http://www.biostars.org/p/84452/";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar84452";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (filename) output file. default: stdout.");
		super.printOptions(out);
		}
	@Override
	public int doWork(String[] args)
		{
		File fileout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:"))!=-1)
			{
			switch(c)
				{
				case 'o': fileout=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		SAMFileWriter sfw=null;
		SAMFileReader sfr=null;
		try
			{
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading sfomr stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				info("Reading from "+filename);
				sfr=new SAMFileReader(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.LENIENT);
			SAMFileHeader header=sfr.getFileHeader();
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setProgramName(getProgramName());
			prg.setProgramVersion(getVersion());
			prg.setCommandLine(getProgramCommandLine());
			
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setCreateMd5File(false);
			sfwf.setCreateIndex(false);
			sfwf.setTempDirectory(getTmpDirectories().get(0));
			
			final boolean sorted=true;
			if(fileout==null)
				{
				info("writing to stdout");
				sfw=sfwf.makeSAMWriter(header,sorted, System.out);
				}
			else
				{
				info("Writing to "+fileout);
				sfw=sfwf.makeSAMOrBAMWriter(header, sorted, fileout);
				}
			long nChanged=0L;
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				if(rec.getReadUnmappedFlag())
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				Cigar cigar=rec.getCigar();
				if(cigar==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				byte bases[]= rec.getReadBases();
				if(bases==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				ArrayList<CigarElement> L=new ArrayList<CigarElement>();
				ByteArrayOutputStream nseq=new ByteArrayOutputStream();
				ByteArrayOutputStream nqual=new ByteArrayOutputStream();
				
				byte quals[]= rec.getBaseQualities();
				int indexBases=0;
				for(CigarElement ce:cigar.getCigarElements())
					{
					switch(ce.getOperator())
						{
						case S: indexBases+=ce.getLength(); break;
						case H://cont
						case P: //cont
						case N: //cont
						case D:
							{
							L.add(ce);
							break;
							}
							
						case I:
						case EQ:
						case X:
						case M:
							{
							L.add(ce);
							nseq.write(bases,indexBases,ce.getLength());
							if(quals.length!=0) nqual.write(quals,indexBases,ce.getLength());
							indexBases+=ce.getLength(); 
							break;
							}
						default:
							{
							throw new PicardException("Unsupported Cigar opertator:"+ce.getOperator());
							}
						}
					
					}
				if(indexBases!=bases.length)
					{
					throw new PicardException("ERRROR "+rec.getCigarString());
					}
				if(L.size()==cigar.numCigarElements())
					{
					sfw.addAlignment(rec);
					continue;
					}
				++nChanged;
				rec.setCigar(new Cigar(L));
				rec.setReadBases(nseq.toByteArray());
				if(quals.length!=0)  rec.setBaseQualities(nqual.toByteArray());
				sfw.addAlignment(rec);
				}
			info("Num records changed:"+nChanged);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar84452().instanceMainWithExit(args);
		}

	}
