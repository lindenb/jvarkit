package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.File;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ReadNameSortMethod;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
/**
 BEGIN_DOC
 
## Example

### Example 1 

``` 
 $ java -jar dist/commbams.jar --samtools \
 	 -f but_metadata -delim '\n' \
 	B00GWFP_std.hg19.qname.bam B00GWFP_S1.hg19.qname.bam
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/1	83	chr5	21564864	40	151M	=	21564540	-475	CTCCCAGAGAGAAGCATCAACAGCTTAGGGTGTAGTCTAAACAGAAATCTTGCACTCCTCCTGCAGTAGCGTCTCTATTTTTTATGCTGAACATTATTTGCTAATTCCAACTGGCTCTAAGCTAATGTGTTTCCCAGGTTTTCTCAATGAN	AFAA<,,,<,,,,,,,,7,,,,7,,,,A7KKF<,F,,7,7,A,A7F7,K<A,,,,7,,,7KKKFAFA,7,A7F7,7,,,KFF,,AKKFFFF<<K<<KAFAKA,A,,A7,AAAFKFA,A,FKKAA,AKKKKFFFKF<KKKKKKKFFAAAA<#
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/1	77	*	0	0	*	*	0	0	NTCATTGAGAAAACCTGGGAAACACATTAGCTTAGAGCCAGTTGGAATTAGCAAATAATGTTCAGCATAAAAAATAGAGACGCTACTGCAGGAGGAGTGCAAGATTTCTGTTTAGACTACACCCTAAGCTGTTGATGCTTCTCTCTGGGA	!<AAAAFFKKKKKKK<FKFFFKKKKA,AAKKF,A,AFKFAAA,7A,,A,AKAFAK<<K<<FFFFKKA,,FFK,,,7,7F7A,7,AFAFKKK7,,,7,,,,A<K,7F7A,A,7,7,,F,<FKK7A,,,,7,,,,7,,,,,,,,<,,,<AAF
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/2	163	chr5	21564540	60	8S106M37S	=	21564864	475	NTAAGAATATTTCACACTTAAAACAAAATCTGATTAGACAAACACTTTGATTGTTATTATTCGCGTATATCATCTACCAGAAGCAAATAGACATCTACTACATCTTTCAAGAAAGTTTACCTATCAATATTACTCAACTGGACCCAATAAT	#<A,<,,A,,K<7FKFF,7,,AF,,7,7AAF,,7<<,,7,AF,,7,7A<,7FA,,7,,7F,,A7FKK7,7,,,,,7,,,,,<,<,,,,7,<,,,,,,7FF7AF<7,,<,,,,7,7,,,,,,,<,,,,,,,,,,,,,,,,,,<,,,,,,,,,
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/2	141	*	0	0	*	*	0	0	NTAAGAATATTTCACACTTAAAACAAAATCTGATTAGACAAACACTTTGATTGTTATTATTCGCGTATATCATCTACCAGAAGCAAATAGACATCTACTACATCTTTCAAGAAAGTTTACCTATCAATATTACTCAACTGGACCCAATAA	!<A,<,,A,,K<7FKFF,7,,AF,,7,7AAF,,7<<,,7,AF,,7,7A<,7FA,,7,,7F,,A7FKK7,7,,,,,7,,,,,<,<,,,,7,<,,,,,,7FF7AF<7,,<,,,,7,7,,,,,,,<,,,,,,,,,,,,,,,,,,<,,,,,,,,
```

### Example 2 

```
 $ java -jar dist/commbams.jar --samtools \
 	B00GWFP_std.hg19.qname.bam B00GWFP_S1.hg19.qname.bam

.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/2
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:6513/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:6513/2
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:7181/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:7181/2
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10205/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10205/2
.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10380/1	.
.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10380/2	.
```


 END_DOC
 */
@Program(name="commbams",
	description="Equivalent of unix 'comm' for bams sorted on queryname",
	keywords={"sam","bam","comm","compare"}
	)
public class CommBams extends Launcher {
	
	private enum WhatToPrint { name,but_metadata,all}
	
	private static final Logger LOG=Logger.build(CommBams.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-1","--hide1"},description="suppress read unique to file 1")
	private boolean hide1=false;
	@Parameter(names={"-2","--hide2"},description="suppress read unique to file 2")
	private boolean hide2=false;
	@Parameter(names={"-3","--hide3"},description="suppress reads present in both files")
	private boolean hide3=false;
	@Parameter(names={"-sortmethod","--sortmethod"},description="[20171110]"+ReadNameSortMethod.DESCRIPTION)
	private ReadNameSortMethod sortMethod = ReadNameSortMethod.picard;
	@Parameter(names={"-delim","--delimiter"},description="Output delimiter")
	private String delim = "\t";
	@Parameter(names={"-empty","--empty"},description="Empty content symbol")
	private String emptySymbol = ".";
	@Parameter(names={"-f","--format"},description="What should I print ? (only the read name ? etc...)")
	private WhatToPrint whatPrint=WhatToPrint.name;
	
	private final Pattern tab=Pattern.compile("[\t]");
	
	private static int side(final SAMRecord rec)
		{
		return CompareBams4.side(rec);
		}
	
	private static String readName(final SAMRecord rec) {
		if(rec.getReadPairedFlag()) {
			if(rec.getFirstOfPairFlag()) return rec.getReadName()+"/1";
			if(rec.getSecondOfPairFlag()) return rec.getReadName()+"/2";
			throw new RuntimeException("Side for flag "+rec.getReadName()+":"+rec.getFlags()+"?");
		}
		else {
			return rec.getReadName();
		}
	}
	
	private void dump(final PrintWriter out,
			final String name,
			int col
			)
		{
		boolean first=true;
		for(int i=1;i<=3;++i)
			{
			if(i==1 && hide1) continue;
			if(i==2 && hide2) continue;
			if(i==3 && hide3) continue;
			if(!first) out.print(this.delim);
			if(i==col)
				{
				out.print(name);
				}
			else
				{
				out.print(this.emptySymbol);
				}
			first=false;
			}
		out.println();
		}
	
	
	
	private void dump(final PrintWriter out,
			final SAMRecord rec0,
			final SAMRecord rec1
			)
		{
		switch(this.whatPrint)
			{
			case name:
				if(rec0==null && rec1!=null)
					{
					if(!hide2) dump(out,readName(rec1),2);
					}
				else if(rec0!=null && rec1==null)
					{
					if(!hide1) dump(out,readName(rec0),1);
					}
				else
					{
					if(!hide3) dump(out,readName(rec0),3);
					}
				break;
			case but_metadata:
			case all:
				{
				for(int side=0;side<2;++side) {
					final SAMRecord r= (side==0?rec0:rec1);
					if(r==null)
						{
						for(int j=0;j< 11;++j)
							{
							if(j>0) out.print('\t');
							out.print(this.emptySymbol);
							}
						continue;
						}
					
					String samStr = r.getSAMString();
					if(samStr.charAt(samStr.length()-1)=='\n')
						{
						samStr= samStr.substring(0,samStr.length()-1);
						}
					if(this.whatPrint.equals(WhatToPrint.all))
						{
						out.print(samStr);
						}
					else
						{
						final String tokens[] = this.tab.split(samStr);
						tokens[0]=readName(r);
						for(int j=0;j< 11;++j)
							{
							if(j>0) out.print('\t');
							out.print(j<tokens.length?tokens[j]:this.emptySymbol);
							}
						}					
					if(side==0) out.print(this.delim);
					}
				out.println();
				break;
				}
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(args.size() !=2)
			{
			LOG.info("Expected two and only two bams please, but got "+args.size());
			return -1;
			}
		final Comparator<SAMRecord> comparator = this.sortMethod.get();

		final SamReader samFileReaders[]={null,null};
		@SuppressWarnings("unchecked")
		final PeekableIterator<SAMRecord> iters[]=new PeekableIterator[]{null,null};
		PrintWriter out=null;
		Optional<SAMRecord> rec0  = Optional.empty();
		Optional<SAMRecord> rec1  = Optional.empty();
		
		if(whatPrint.equals(WhatToPrint.name) && hide1 && hide2 && hide3) {
			LOG.error("all flags hide** are on");
			return -1;
		}
		try
			{
			if( this.delim.equals("\\n")) this.delim = "\n";
			for(int i=0;i< args.size() && i< samFileReaders.length;++i)
				{
				final String samFile=args.get(i);
				LOG.info("opening "+samFile);
				samFileReaders[i]=super.openSamReader(samFile);
				final SAMFileHeader header = samFileReaders[i].getFileHeader();
				if(header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
					LOG.error("Expected "+samFile+" to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder());
					return -1;
					}
				
				iters[i] = new PeekableIterator<>(samFileReaders[i].iterator());
				}
			out= super.openFileOrStdoutAsPrintWriter(outputFile);
			
			for(;;) {
				for(int i=0;i< 2;++i) {
					Optional<SAMRecord> optRec = (i==0?rec0:rec1);
					if(!optRec.isPresent())
						{
						while(iters[i].hasNext()) {
							final SAMRecord rec = iters[i].peek();
							
							if(optRec.isPresent() && comparator.compare(optRec.get(),rec)>0)
								{
								LOG.error("Something is wrong in sort order of "+args.get(i)+" : got\n\t"
										+rec+" "+side(rec)+"\nafter\n\t"+ optRec.get()+" "+side(optRec.get())+"\nSee also option (samtools querysort)"
										);
								return -1;
								}
							// equals
							else if( optRec.isPresent()&& comparator.compare(optRec.get(),rec)==0)
								{
								iters[i].next();//consumme
								}
							//it's a new
							else if(!optRec.isPresent()){
								optRec = Optional.of(iters[i].next());//consumme
								if(i==0)
									{
									rec0=optRec;
									}
								else
									{
									rec1=optRec;
									}
								}
							// compare <0
							else
								{	
								break;
								}
							}
						}
				}
				
				
				if(!rec0.isPresent() && !rec1.isPresent()) break;
	
						
				if((!rec0.isPresent() && rec1.isPresent()) ||
				   (rec0.isPresent() && rec1.isPresent() && comparator.compare(rec0.get(),rec1.get())>0)
					)
					{
					dump(out,null,rec1.get());
					rec1 =Optional.empty();
					}
				else if((rec0.isPresent() && !rec1.isPresent()) ||
						(rec0.isPresent() && rec1.isPresent() && comparator.compare(rec0.get(), rec1.get())<0))
					{
					dump(out,rec0.get(),null);
					rec0 =Optional.empty();
					}
				else
					{
					dump(out,rec0.get(),rec1.get());
					rec0 =Optional.empty();
					rec1 =Optional.empty();
					}
				}
			out.flush();
			out.close();out=null;
			return 0;
			}
		catch(Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			for(int i=0;i< samFileReaders.length;++i){
				CloserUtil.close(iters[i]);
				CloserUtil.close(samFileReaders[i]);
				}
			CloserUtil.close(out);
		}
		}
	
	public static void main(final String[] args) {
		new CommBams().instanceMainWithExit(args);

	}

}
