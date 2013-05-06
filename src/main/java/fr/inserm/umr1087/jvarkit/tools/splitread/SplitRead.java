package fr.inserm.umr1087.jvarkit.tools.splitread;

import java.io.File;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.TextCigarCodec;

public class SplitRead {
	private static final Logger LOG=Logger.getLogger(SplitRead.class.getSimpleName());
	private final Pattern semiColon=Pattern.compile("[;]");
	private final Pattern comma=Pattern.compile("[,]");
	private float maxFractionCommon=0.1f;
	
	private class Fragment
		{
		String chrom;
		int pos;
		char strand;
		String cigar;
		
		public int compareTo(Fragment other)
			{
			int i=chrom.compareTo(other.chrom);
			if(i!=0) return i;
			return pos-other.pos;
			}
		
		void print()
			{
			System.out.print(chrom+"\t"+pos+"\t"+strand+"\t"+cigar);
			}
		}
	
	
	private void scanRecord(final SAMRecord record) throws Exception
		{
		if(record.getReadUnmappedFlag()) return;
		String xp=record.getStringAttribute("XP");
		if(xp==null) return;
		LOG.info(xp);
		Cigar cigar1=record.getCigar();
		int readPos=0;
		int readMap1[]=new int[record.getReadLength()];
		for(CigarElement ce:cigar1.getCigarElements())
			{
			switch(ce.getOperator())
				{
				case I: case S:
					{
					readPos+=ce.getLength();
					break;
					}
				case M:case X:case EQ:
					{
					for(int i=0;i< ce.getLength();++i)
						{
						readMap1[readPos]+=1;
						readPos++;
						}
					break;
					}
				case P: case H: case D: case N: break;
				default: throw new RuntimeException("cigar operator not handled:"+ce.getOperator());
				}
			}
		
		for(String s:this.semiColon.split(xp))
			{
			if(s.isEmpty()) continue;
			
			
			String tokens[]=this.comma.split(s);
			Cigar cigar2=TextCigarCodec.getSingleton().decode(tokens[2]);
			
			readPos=0;
			float common=0f;
			for(CigarElement ce:cigar2.getCigarElements())
				{
				switch(ce.getOperator())
					{
					case I: case S:
						{
						readPos+=ce.getLength();
						break;
						}
					case M:case X:case EQ:
						{
						for(int i=0;i< ce.getLength();++i)
							{
							if(readMap1[readPos]==1)
								{
								common++;
								}
							readPos++;
							}
						break;
						}
					case P: case H: case D: case N: break;
					default: throw new RuntimeException("cigar operator not handled:"+ce.getOperator());
					}
				}
			float fraction=common/readMap1.length;
			if(  fraction > this.maxFractionCommon)
				{
				continue;
				}

			
			Fragment f1=new Fragment();
			f1.chrom=record.getReferenceName();
			f1.pos=record.getAlignmentStart();
			f1.strand=(record.getReadNegativeStrandFlag()?'-':'+');
			f1.cigar=record.getCigarString();
			
			Fragment f2=new Fragment();
			f2.chrom=tokens[0];
			f2.pos=Integer.parseInt(tokens[1].substring(1));
			f2.strand=tokens[1].charAt(0);
			f2.cigar=tokens[2];

			System.out.print(
				record.getReadName()+"\t"+
				(record.getFirstOfPairFlag()?'F':'R')+"\t"
				);
			if(f1.compareTo(f2)<0)
				{
				f1.print();
				System.out.print("\t");
				f2.print();
				}
			else
				{
				f2.print();
				System.out.print("\t");
				f1.print();
				}	
			System.out.println("\t"+fraction);
			}
		}
	

	private void scan(SAMFileReader reader) throws Exception
		{
		long nrecords=0L;
		for(Iterator<SAMRecord> iter=reader.iterator();
				iter.hasNext(); )
			{
			SAMRecord record=iter.next();
			++nrecords;
			if(nrecords%1E6==0)
				{
				LOG.info("nRecord:"+nrecords);
				System.out.flush();
				}
			scanRecord(record);
			}
		}
	
	
	private void run(String[] args)
		throws Exception
		{
		
		int optind=0;
		while(optind< args.length)
			{
			if(args[optind].equals("-h") ||
			   args[optind].equals("-help") ||
			   args[optind].equals("--help"))
				{
				System.err.println("Pierre Lindenbaum PhD. 2013");
				System.err.println("Options:");
				System.err.println(" -h help; This screen.");
				System.err.println(" -R (reference file) REQUIRED.");
				System.err.println(" -L|--level (log level): default:"+LOG.getLevel());
;
				}
			
			else if((args[optind].equals("-L") || args[optind].equals("--level")) && optind+1< args.length)
				{
				LOG.setLevel(Level.parse(args[++optind]));
				}

			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unknown option "+args[optind]);
				return;
				}
			else 
				{
				break;
				}
			++optind;
			}
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.LENIENT);
		
		if(optind==args.length)
			{
			SAMFileReader r=new SAMFileReader(System.in);
			scan(r);
			r.close();
			}
		else if(optind+1==args.length)
			{
			File file=new File(args[optind++]); 
			SAMFileReader r=new SAMFileReader(file);
			scan(r);
			r.close();
			}
		else 
			{
			System.err.println("illegal number of arguments.");
			System.exit(-1);
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		LOG.setLevel(Level.OFF);
		new SplitRead().run(args);
		}
	
	}
