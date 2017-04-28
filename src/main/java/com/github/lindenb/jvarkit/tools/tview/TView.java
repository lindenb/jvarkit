package com.github.lindenb.jvarkit.tools.tview;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
//import htsjdk.samtools.util.Log;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;


public class TView implements Closeable
	{
	private static final Logger LOG = Logger.build(TView.class).make();
	
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	private enum AnsiColor {
    	BLACK (30),
    	RED (31),
    	GREEN (32),
    	YELLOW (33),
    	BLUE (34),
    	MAGENTA (35),
    	CYAN (36),
    	WHITE (37)
		;
    	
    	AnsiColor(final int opcode) {
    		this.opcode=opcode;
    		}
    	final int opcode;
    	int pen() { return (opcode);}
    	int paper() { return (opcode+10);}
    	}
	

	private class Colorizer
		{
		protected Set<Integer> opcodes=new HashSet<>();
		protected PrintStream out;
		
		Colorizer(PrintStream out) {this.out=out;}
		public Colorizer pen(AnsiColor c) {
			if(c==null)
				{
				for(AnsiColor ansi:AnsiColor.values()) opcodes.remove(ansi.pen());
				}
			else
				{
				opcodes.add(c.pen());
				}
			return this;
		}
		public Colorizer paper(AnsiColor c) {
			if(c==null)
			{
			for(AnsiColor ansi:AnsiColor.values()) opcodes.remove(ansi.paper());
			}
		else
			{
			opcodes.add(c.paper());
			}
		return this;
	}
		public Colorizer bold(boolean b) {
			if(b) {opcodes.add(1);}
			else { opcodes.remove(1);}
			return this;
		}
		
		public Colorizer underscore(boolean b) {
			if(b) {opcodes.add(4);}
			else { opcodes.remove(4);}
			return this;
		}
		
		
		public Colorizer print(final Object o)
			{
			out.print(o);
			return this;
			}
		public Colorizer println(final Object o)
			{
			this.print(o);
			out.println();
			return this;
			}
		}
	
	private class AnsiColorizer extends Colorizer {
		AnsiColorizer(PrintStream out) {super(out);}

		@Override
		public Colorizer print(Object o) {
			if(!opcodes.isEmpty())
				{
				final StringBuilder sb=new StringBuilder(ANSI_ESCAPE);
				sb.append(opcodes.stream().map(I->String.valueOf(I)).collect(Collectors.joining(";")));
				sb.append("m");
				out.print(sb.toString());
				}
			
			out.print(o);
			if(!opcodes.isEmpty()) {
				this.opcodes.clear();
				out.print(ANSI_RESET);
				}
			
			return this;
			}
		}
	
	
	private boolean showClip=false;
	private Interval interval=null;
	@Parameter(names={"-R","--reference"},description="Indexed Reference file.",required=true)
	private File referenceFile=null;
	@Parameter(names={"--insert"},description="Show insertions")
	private boolean showInsertions=false;
	@Parameter(names={"--readName"},description="Show read name")
	private boolean showReadName=false;
	@Parameter(names={"--plain"},description="Plain text, disable ANSI colors")
	private boolean disableANSIColors=false;
	@Parameter(names={"--hideBases"},description="Hide bases")
	private boolean hideBases=false;
	@Parameter(names={"-V","--variant","--variants","--vcf"},description="Variant file. "+IOUtils.UNROLL_FILE_MESSAGE)
	private File variantFiles = null;
	@Parameter(names=SamFilterParser.DEFAULT_OPT,description=SamFilterParser.FILTER_DESCRIPTION)
	private SamRecordFilter samRecordFilter = SamFilterParser.buildDefault();

	
	private int distance_between_reads=2;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private final List<SamInputResource> samInputResources=new ArrayList<>();
	private final List<SamReader> samReaders=new ArrayList<>();
	private final List<VCFFileReader> vcfReaders=new ArrayList<>();

	
	public TView() {
		
	}
	
	
	public void setBamFiles(List<SamInputResource> bamFiles)
		{
		this.samInputResources.clear();
		this.samInputResources.addAll(bamFiles);
		}
	
	public void setInterval(final Interval interval)
		{
		this.interval = interval;
		}

	
	
	public int initialize() throws IOException
		{
		if(this.referenceFile!=null) {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.referenceFile);
			}
		if(this.samRecordFilter==null) {
			this.samRecordFilter = SamFilterParser.ACCEPT_ALL;
		}
		final SamReaderFactory srf=SamReaderFactory.makeDefault().
				referenceSequence(this.referenceFile).
				validationStringency(ValidationStringency.LENIENT)
				; 
		for(final SamInputResource sir:this.samInputResources)
			{
			final SamReader samReader= srf.open(sir);
			this.samReaders.add(samReader);
			}
		
		for(final File vcfFile:IOUtils.unrollFile(this.variantFiles))
			{
			final VCFFileReader vcfReader= new VCFFileReader(vcfFile,true);
			this.vcfReaders.add(vcfReader);
			}
		
		return 0;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.indexedFastaSequenceFile);
		this.indexedFastaSequenceFile=null;
		
		for(final SamReader r: this.samReaders)
			{
			CloserUtil.close(r);
			}

		for(final VCFFileReader r: this.vcfReaders)
			{
			CloserUtil.close(r);
			}
		this.samInputResources.clear();
		this.samReaders.clear();
		this.vcfReaders.clear();
		}
	
	public SamRecordFilter getRecordFilter() {
		return this.samRecordFilter;
	}
	
	public Function<SAMRecord,String> getRecordGroup() {
		
		return R->{
			final SAMReadGroupRecord rg= R.getReadGroup();
			if(rg!=null && rg.getSample()!=null) return rg.getSample();
			return "UNDEFINED";
		};
	}
	
	
	public Function<SAMRecord,Integer> left() {
		return R->showClip?
			R.getUnclippedStart():
			R.getAlignmentStart();
	}
	public Function<SAMRecord,Integer> right() {
		return R->showClip?
			R.getUnclippedEnd():
			R.getAlignmentEnd();
	}

	
	void paint(final PrintStream out) {
		final Colorizer colorizer = disableANSIColors?
				new Colorizer(out):
				new  AnsiColorizer(out)
				;
				
		if(interval==null) 
			{
			LOG.warn("No interval defined");
			return;
			}
		final GenomicSequence contigSequence;
		final Function<Integer, Character> refPosToBase;
		if(indexedFastaSequenceFile!=null)
			{
				
			if(indexedFastaSequenceFile.getSequenceDictionary().getSequence(this.interval.getContig())==null)
				{
				LOG.warn("No interval with contig "+interval+" in REF");
				return;
				}
			contigSequence= new GenomicSequence(indexedFastaSequenceFile, interval.getContig());
			refPosToBase = POS->{
				if(POS<0 || POS >= contigSequence.length()) return 'N';
				return contigSequence.charAt(POS);
				};
			}
		else
			{
			contigSequence = null;
			refPosToBase = POS -> 'N';
			}
		
		/** test if genomic position is in interval */
		final Predicate<Integer> testInInterval = new Predicate<Integer>() {
			@Override
			public boolean test(final Integer pos) {
				return interval.getStart() <= pos && pos <= interval.getEnd();
			}
		};
		
		final int pixelWidth= this.interval.length() ;
		final Map<Integer, Integer> genomicpos2insertlen = new TreeMap<>();
		

		final Map<String, List<SAMRecord>> group2record=new TreeMap<>();
		
		for(final SamReader samReader:this.samReaders)
			{
			SAMRecordIterator iter = samReader.query(
					this.interval.getContig(),
					this.interval.getStart(),
					this.interval.getEnd(),
					false
					);
			while(iter.hasNext())
				{
				final SAMRecord rec = iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getCigar()==null) continue;
				if(getRecordFilter().filterOut(rec)) continue;
				if( !rec.getContig().equals(interval.getContig())) continue;
				if(right().apply(rec) < this.interval.getStart()) continue;
				if(this.interval.getEnd() < left().apply(rec) ) continue;
				final String group = getRecordGroup().apply(rec);
				if(group==null || group.isEmpty()) continue;
				List<SAMRecord> records = group2record.get(group);
				if( records == null) {
					records = new ArrayList<>();
					group2record.put(group,records);
					}
				records.add(rec);
				
				//loop over cigar, get the longest insert
				int refpos=rec.getAlignmentStart();
				for(final CigarElement ce:rec.getCigar().getCigarElements()) {
					if(!this.showInsertions) break;
					final CigarOperator op = ce.getOperator();

					if(op.equals(CigarOperator.I) && testInInterval.test(refpos))
						{
						final Integer longestInsert= genomicpos2insertlen.get(refpos);
						if(longestInsert==null|| longestInsert.compareTo(ce.getLength())<0)
							{
							genomicpos2insertlen.put(refpos, ce.getLength());
							}
						}
					if(op.consumesReferenceBases())
						{
						refpos += ce.getLength();
						}
					if(refpos > interval.getEnd()) break;
					}
				}
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			}
		

		
		
		
		/** compute where are the insertions */
		LOG.debug(genomicpos2insertlen);
		final Predicate<Integer> insertIsPresentAtX = SCREENX ->{
			int x=0;
			int ref= interval.getStart();
			while(x < pixelWidth)
				{
				if(x>SCREENX) return false;
				final Integer insertLen = genomicpos2insertlen.get(ref);
				if(insertLen==null)
					{
					++x;
					++ref;
					}
				else
					{
					if(x<=SCREENX && SCREENX <x+insertLen) return true;
					x+=(insertLen+1);//(+1) I DON'T UNDERSTAND WHY, BUT IT WORKS
					++ref;
					}
				}
			return false;
			};

		final Function<Character,AnsiColor> base2ansiColor = BASE->{
			switch(Character.toUpperCase(BASE))
				{
				case 'A': return AnsiColor.BLUE;
				case 'T': return AnsiColor.GREEN;
				case 'G': return AnsiColor.CYAN;
				case 'C': return AnsiColor.YELLOW;
				default: return null;
				}
			};
		
		/** paint base position */
		int ref = this.interval.getStart();
		int x=0;
		
		while(x < pixelWidth)
			{
			
			if(insertIsPresentAtX.test(x))
				{
				colorizer.pen(AnsiColor.RED).print("^");
				++x;
				}
			else if((ref-this.interval.getStart())%10==0)
				{
				final String f=String.format("%d", ref);
				for(int i=0;i< f.length() && x < pixelWidth;++i)
					{
					colorizer.pen(AnsiColor.GREEN).print(f.charAt(i));
					if(!insertIsPresentAtX.test(x)) ++ref;
					++x;
					}
				}
			else
				{
				out.print(".");
				++ref;
				++x;
				}
			}
		out.println();
		
		
		/* paint ref base */
		ref = this.interval.getStart();
		x=0;
		while(x < pixelWidth)
			{
			if(insertIsPresentAtX.test(x))
				{
				colorizer.paper(AnsiColor.YELLOW).print("*");
				++x;
				}
			else
				{
				char refBase= refPosToBase.apply(ref-1);
				colorizer.pen(base2ansiColor.apply(refBase)).print(refBase);
				++ref;
				++x;
				}
			}
		out.println();

		
		for(x=0;x< pixelWidth;++x)
			{
			out.print(insertIsPresentAtX.test(x)?"^":" ");
			}
		out.println();
		
		
		for(final String groupName : group2record.keySet())
			{
			colorizer.pen(AnsiColor.BLUE).println(groupName);
			final List<List<SAMRecord>> rows = new ArrayList<>();
			
			/* pileup reads */
			for(final SAMRecord rec: group2record.get(groupName)) {
				int y=0;
				for(y=0;y< rows.size();++y)
					{
					final List<SAMRecord> row = rows.get(y);
					final SAMRecord last = row.get(row.size()-1);
					if(right().apply(last) + this.distance_between_reads < left().apply(rec)) {
						row.add(rec);
						break;
						}
					}
				if(y==rows.size()) {
					final List<SAMRecord> row = new ArrayList<>();
					row.add(rec);
					rows.add(row);
					}
				}
			/* print each row */
			for(final List<SAMRecord> row : rows)
				{
				ref = interval.getStart();

				x=0;
				
				for(final SAMRecord rec: row) {
					int readRef=  left().apply(rec);
					// pad before record
					while( x < pixelWidth &&
							ref <  readRef &&
							testInInterval.test(ref))
						{
						if(!insertIsPresentAtX.test(x) ) ++ref;
						++x;
						out.print(' ');
						}
										
					int readpos=0;
					
					/* get read base function */
					final Function<Integer,Character> baseAt=new Function<Integer, Character>() {
						@Override
						public Character apply(final Integer readpos) {
							if(readpos<0 || readpos>=rec.getReadLength()) return '?';
							return (char) rec.getReadBases()[readpos];
							}
						};
					
					for(final CigarElement ce:rec.getCigar())
						{
						final CigarOperator op = ce.getOperator();
						if(op.equals(CigarOperator.PADDING)) continue;
						
						if(this.showInsertions && op.equals(CigarOperator.I)) {
							int cigarIdx =0;
							while( x < pixelWidth &&
									cigarIdx < ce.getLength()
									)
								{
								if(testInInterval.test(readRef)) {
									colorizer.paper(AnsiColor.RED).print(baseAt.apply(readpos));
									++x;
									}
								++cigarIdx;
								++readpos;
								}
							continue;
							}

												
						int cigarIdx =0;
						while( x < pixelWidth &&
								cigarIdx < ce.getLength() 
								)
							{
							colorizer.opcodes.clear();
							//pad before base
							while( x < pixelWidth &&
									testInInterval.test(readRef) &&
									(insertIsPresentAtX.test(x))
									)
								{
								++x;
								colorizer.paper(AnsiColor.YELLOW).print("*");
								continue;
								}

							
							switch(op)
								{
								case I: 
									{
									//if visible, processed above
									if(showInsertions) throw new IllegalStateException();
									readpos++;
									break;
									}
								case P: break;
								case H:
									{
									if(showClip)
										{
										if(testInInterval.test(readRef)) {
											colorizer.paper(AnsiColor.YELLOW).print('N');
											++x;
											}
										++readRef;
										}
									break;
									}
								case S:
									{
									if(showClip)
										{
										if(testInInterval.test(readRef)) {
											colorizer.paper(AnsiColor.YELLOW).print((char)rec.getReadBases()[readpos]);
											++x;
											}
										++readpos;
										++readRef;
										}
									else
										{
										readpos++;
										}	
									break;
									}
								case D: case N:
									{
									if(testInInterval.test(readRef)) {
										colorizer.paper(AnsiColor.RED).print('-');
										++x;
										}
									++readRef;
									break;
									}
								case EQ: case M : case X:
									{
									if(testInInterval.test(readRef)) {
										final char refBase =  Character.toUpperCase(refPosToBase.apply(readRef-1));
										char readBase = Character.toUpperCase(baseAt.apply(readpos));
										
										colorizer.pen(base2ansiColor.apply(readBase));
										
										if(op.equals(CigarOperator.X) || (refBase!='N' &&  readBase!='N'  && readBase!=refBase))
											{
											colorizer.pen(AnsiColor.RED);
											}
										else if(hideBases)
											{
											
											if(rec.getReadNegativeStrandFlag())
												{
												readBase=',';
												}
											else
												{
												readBase='.';
												}
											}
										
										if(showReadName)
											{
											final String readName = rec.getReadName();
											if(readpos<0 || readpos>=readName.length()) {
												readBase = '_';
												}
											else
												{
												readBase = readName.charAt(readpos);
												}
											}
										if(rec.getReadNegativeStrandFlag())
											{
											readBase=Character.toLowerCase(readBase);
											}
										else
											{
											readBase=Character.toUpperCase(readBase);
											}
										colorizer.print(readBase);
										++x;
										}
									++readpos;
									++readRef;
									break;
									}
								}
							++cigarIdx;
							}
						}//end of loop cigar
					ref=readRef;
					}//end of loop read
				
				//out.println( " "+ref+" "+row.get(0).getAlignmentStart()+" "+row.get(0).getCigarString()+" "+row.get(0).getReadString());
				while(x<pixelWidth)
					{
					out.print(" ");
					++x;
					}
				out.println();
				if(out.checkError()) break;
				}
			if(out.checkError()) break;
			}
		
		/** variant section*/
		if(this.vcfReaders.isEmpty() && !out.checkError()) {
			out.println();
			for(final VCFFileReader r:vcfReaders)
				{
				if(out.checkError()) break;
				final VCFHeader header = r.getFileHeader();
				final CloseableIterator<VariantContext> iter = r.query(this.interval.getContig(), interval.getStart(), interval.getEnd());
				final List<VariantContext> variants = new ArrayList<>();
				while(iter.hasNext())
					{}
				iter.close();
				if(header.hasGenotypingData())
					{}
				else
					{}
				}
			}
		LOG.debug("done");
		}
	
    
	}
