package com.github.lindenb.jvarkit.tools.samfixcigar;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.Deflater;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileReader.ValidationStringency;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

public class SamFixCigar extends AbstractCommandLineProgram
	{
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	@Override
	public String getProgramDescription() {
		return "Fix Cigar String in SAM replacing 'M' by 'X' or '='";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SamFixCigar";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -r (file) reference indexed with samtools faidx . Required.");
		out.println(" -o (file) BAM/SAM fileout. default:stdout.");
		out.println(" -C (int) set zip compression level.");

		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		GenomicSequence genomicSequence=null;
		int maxRecordsInRam=10000;
		File faidx=null;
		File fout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, getGetOptDefault()+"r:o:C:"))!=-1)
			{
			switch(c)
				{
				case 'r': faidx=new File(opt.getOptArg());break;
				case 'o': fout=new File(opt.getOptArg());break;
				case 'C': 
					{
					BlockCompressedOutputStream.setDefaultCompressionLevel(
							Math.max(Deflater.NO_COMPRESSION, Math.min(Deflater.BEST_COMPRESSION, Integer.parseInt(opt.getOptArg())))
							);
					break;
					}
				default: switch(handleOtherOptions(c, opt, null))
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			}
		
		if(faidx==null)
			{
			error("Reference was not specified.");
			return -1;
			}
		long nReads=0L;
		long nX=0L;
		SAMFileReader sfr=null;
		SAMFileWriter sfw=null;
		SAMFileHeader header;
		try
			{
			info("Loading reference");
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File fin=new File(args[opt.getOptInd()]);
				info("Reading "+fin);
				sfr=new SAMFileReader(fin);
				}
			else
				{
				error("Illegal number of arguments");
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.LENIENT);
			header=sfr.getFileHeader();
			
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setCreateIndex(false);
			sfwf.setCreateMd5File(false);
			sfwf.setMaxRecordsInRam(maxRecordsInRam);
			
			if(fout==null)
				{
				sfw=sfwf.makeSAMWriter(header, header.getSortOrder()==SortOrder.coordinate,System.out);
				}
			else
				{
				sfw=sfwf.makeSAMOrBAMWriter(header, header.getSortOrder()==SortOrder.coordinate,fout);
				}
			List<CigarElement> newCigar=new ArrayList<CigarElement>();
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				++nReads;
				if(nReads%1E6==0)
					{
					info("Reads "+nReads+" operator:X="+nX);
					}
				SAMRecord rec=iter.next();
				Cigar cigar=rec.getCigar();
				byte bases[]=rec.getReadBases();
				if( rec.getReadUnmappedFlag() ||
					cigar==null ||
					cigar.getCigarElements().isEmpty() ||
					bases==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				if(genomicSequence==null ||
					genomicSequence.getSAMSequenceRecord().getSequenceIndex()!=rec.getReferenceIndex())
					{
					genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
					}
				
				newCigar.clear();
				int refPos1=rec.getAlignmentStart();
				int readPos0=0;
				
				
				for(CigarElement ce:cigar.getCigarElements())
					{
					switch(ce.getOperator())
						{
    					case H:break;

    					
    					case P://cont
						case N://cont
						case D:
							{
							newCigar.add(ce);
							refPos1+=ce.getLength();
							break;
							}
						case S:
						case I:
							{
							newCigar.add(ce);
							readPos0+=ce.getLength();							
							break;
							}
						case EQ://cont
						case X:
							{
							newCigar.add(ce);
							refPos1+=ce.getLength();
							readPos0+=ce.getLength();	
							break;
							}
						case M:
							{
							boolean X=false;
							for(int i=0;i< ce.getLength();++i)
	    		    			{
								char c1=Character.toUpperCase((char)bases[readPos0]);
								char c2=Character.toUpperCase(refPos1-1< genomicSequence.length()?genomicSequence.charAt(refPos1-1):'*');
								
								if(c2=='N' || c1==c2)
									{
									newCigar.add(new CigarElement(1, CigarOperator.EQ));
									}
								else
									{
									newCigar.add(new CigarElement(1, CigarOperator.X));
									X=true;
									}
								
	    						refPos1++;
	    						readPos0++;
    		    				}
							if(X) nX++;
							break;
							}
						default: throw new PicardException("Cannot parse cigar "+rec.getCigarString()+" in "+rec.getReadName());
						}
					}
				int i=0;
				while(i< newCigar.size())
					{
					CigarOperator op1 = newCigar.get(i).getOperator();
					int length1 = newCigar.get(i).getLength();
					
					if( i+1 <  newCigar.size() &&
						newCigar.get(i+1).getOperator()==op1)
						{
						CigarOperator op2= newCigar.get(i+1).getOperator();
						int length2=newCigar.get(i+1).getLength();

						 newCigar.set(i,new CigarElement(length1+length2, op2));
						 newCigar.remove(i+1);
						}
					else
						{
						++i;
						}
					}
				cigar=new Cigar(newCigar);
				String newCigarStr=TextCigarCodec.getSingleton().encode(cigar);
				//info("changed "+rec.getCigarString()+" to "+newCigarStr+" "+rec.getReadName()+" "+rec.getReadString());
				rec.setCigar(cigar);
				rec.setCigarString(newCigarStr);
				
				sfw.addAlignment(rec);
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
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamFixCigar().instanceMainWithExit(args);

	}

}
