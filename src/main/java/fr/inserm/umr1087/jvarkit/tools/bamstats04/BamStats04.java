package fr.inserm.umr1087.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.logging.Logger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamStats04
	{
	private static final Logger LOG=Logger.getLogger("bamstats04");
	private File bedFile=null;
	private boolean skipDuplicates=true;
	private int minQual=0;
	private int basesperbin=10;
	private int num_bin=20;
	private boolean cumulative=true;
	private BamStats04()
		{
		
		
		}
	
	private void scan(File bam) throws Exception
		{
		LOG.info("Scanning "+bam);
		long bases2count[]=new long[num_bin];
		Arrays.fill(bases2count, 0L);
		Pattern tab=Pattern.compile("[\t]");
		String tokens[];
		BufferedReader in=new BufferedReader(new FileReader(this.bedFile));
		SAMFileReader samReader = new SAMFileReader(bam);
		long total=0L;
		String line=null;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			tokens=tab.split(line,5);
			if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.bedFile);
			String chrom=tokens[0];
			int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
			int chromEnd1= Integer.parseInt(tokens[2]);
			/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
			
			int counts[]=new int[chromEnd1-chromStart1+1];
			Arrays.fill(counts, 0);
			
			/**
			 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
			 */
			SAMRecordIterator r=samReader.queryOverlapping(chrom, chromStart1, chromEnd1);
			while(r.hasNext())
				{
				SAMRecord rec=r.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(skipDuplicates && rec.getDuplicateReadFlag() ) continue;
				if(!rec.getReferenceName().equals(chrom)) continue;
				if(rec.getMappingQuality()< this.minQual) continue;
				int end=rec.getAlignmentEnd();
				for(int pos=rec.getAlignmentStart();
						pos<=end;
						++pos)
					{
					if(pos< chromStart1) continue;
					if(pos> chromEnd1) continue;
					counts[pos-chromStart1]++;
					}
				}
			
			r.close();
			
			for(int depth:counts)
				{
				int cat=depth/this.basesperbin;
				if(cat>=bases2count.length) cat=bases2count.length-1;
				if(cumulative)
					{
					for(int j=0;j<=cat;++j)
						{
						bases2count[j]++;
						}
					}
				else
					{
					bases2count[cat]++;
					}
				++total;
				}
			
			}
		in.close();
		samReader.close();
		System.out.print(bam.toString()+"\t"+total);
		for(long bases:bases2count)
			{
			System.out.print("\t");
			System.out.print(bases);
			}
		System.out.println();
		}

	
	private void run(String args[]) throws Exception
		{
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				System.out.println("Pierre Lindenbaum PhD. 2013.");
				System.out.println(" -b bedfile (required).");
				System.out.println(" -m (int) min-qual.");
				System.out.println(" -D do NOT ignore duplicates.");
				System.out.println(" -c NOT cumulative results.");
				return;
				}
			else if(args[optind].equals("-b") && optind+1< args.length)
				{
				this.bedFile=new File(args[++optind]);
				}
			else if(args[optind].equals("-m") && optind+1< args.length)
				{
				this.minQual=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-D"))
				{
				this.skipDuplicates=false;
				}
			else if(args[optind].equals("-c"))
				{
				this.cumulative=false;
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return;
				}
			else
				{
				break;
				}
			++optind;
			}
		if(optind==args.length)
				{
				System.err.println("Bam input missing");
				System.exit(-1);
				}
		
		System.out.print("#filename\ttotal_bases");
		

		for(int i=0;i< this.num_bin;++i)
			{
			System.out.print("\t[" + (i*this.basesperbin));
			if(i+1==this.num_bin)
				{
				System.out.print("-all[");
				}
			else
				{
				System.out.print("-" + ((i+1)*this.basesperbin) + "[");
				}
			}
		System.out.println();

		while(optind< args.length)
			{
			File bam=new File(args[optind++]);
			scan(bam);
			}
				
		
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception
		{
		BamStats04 app=new BamStats04();
		app.run(args);
		}

}

