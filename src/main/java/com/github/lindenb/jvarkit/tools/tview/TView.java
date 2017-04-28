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
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
//import htsjdk.samtools.util.Log;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


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
			opcodes.add(c.pen());
			return this;
		}
		public Colorizer paper(AnsiColor c) {
			opcodes.add(c.paper());
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
		public Colorizer print(Object o)
			{
			out.print(o);
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
				out.print(ANSI_RESET);
				}
			opcodes.clear();
			return this;
			}
		}
	
	
	private final static int INSERTION_OPCODE = -1; 
    private final static int NOT_VISIBLE = -1; 
	private boolean showClip=false;
	private Interval interval=null;
	@Parameter(names={"-R","--reference"},description="Indexed Reference file.",required=true)
	private File referenceFile=null;
	@Parameter(names={"--insert"},description="Show insertions")
	private boolean showInsertions=false;
	@Parameter(names={"--readname"},description="Show read name")
	private boolean showReadName=false;
	
	private int distance_between_reads=2;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private final List<SamInputResource> samInputResources=new ArrayList<>();
	private final List<SamReader> samReaders=new ArrayList<>();

	
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
		final SamReaderFactory srf=SamReaderFactory.makeDefault().
				referenceSequence(this.referenceFile).
				validationStringency(ValidationStringency.LENIENT)
				; 
		for(final SamInputResource sir:this.samInputResources)
			{
			 final SamReader samReader= srf.open(sir);
			this.samReaders.add(samReader);
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
		this.samInputResources.clear();
		this.samReaders.clear();
		}
	
	public Predicate<SAMRecord> getRecordFilter() {
		return R->true;
	}
	
	public Function<SAMRecord,String> getRecordGroup() {
		return R->"ALL";
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
		Colorizer colorizer = new Colorizer(out);//new AnsiColorizer(out);
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
				if(!getRecordFilter().test(rec)) continue;
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
				}
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			}
		
		/** test if genomic position is in interval */
		final Predicate<Integer> testInInterval = new Predicate<Integer>() {
			@Override
			public boolean test(final Integer pos) {
				return interval.getStart() <= pos && pos <= interval.getEnd();
			}
		};
		
		/** convert pixel to genomic position , value can be INSERTION_OPCODE */
		final int pixelToRef[] = new int[ this.interval.length() ];
		
		
		/** compute where are the insertions */
		int x=0;
		int ref = interval.getStart();
		while( x < pixelToRef.length )
			{
			Integer longestInsert = null;
			for(final SAMRecord rec: group2record.values().stream().
					flatMap(G->G.stream()).collect(Collectors.toList())) 
				{	
				if(!showInsertions) break;
				if(rec.getAlignmentEnd() < ref) continue;
				if(ref < rec.getAlignmentStart()) continue;
				int readref = rec.getAlignmentStart();
				
				for(final CigarElement ce:rec.getCigar().getCigarElements()) {
					final CigarOperator op = ce.getOperator();

					if(op.equals(CigarOperator.I) && 
						ref==readref &&
						(longestInsert==null || longestInsert.compareTo(ce.getLength())<0)
						)
						{
						longestInsert = ce.getLength();
						}
					
					if(op.consumesReferenceBases())
						{
						readref += ce.getLength();
						if(readref > ref) break;
						}
					}
				}
			if(!this.showInsertions || longestInsert==null)
				{
				pixelToRef[x]=ref;
				++x;
				++ref;
				}
			else
				{
				LOG.debug(longestInsert+" at "+ref);
				for(int i=0;i<longestInsert && x< pixelToRef.length;++i )
					{
					pixelToRef[x]=INSERTION_OPCODE;
					++x;
					}
				++ref;
				}
			}
		
		
		//final Map<Integer,Integer> posToInsertLength =new HashMap<>(this.interval.length());
		

		
		
		/** paint base position */
		ref = this.interval.getStart();
		x=0;
		while(x < pixelToRef.length)
			{
			
			if(pixelToRef[x]==INSERTION_OPCODE)
				{
				out.print("^");
				++x;
				}
			else if((ref-this.interval.getStart())%10==0)
				{
				final String f=String.format("%d", ref);
				for(int i=0;i< f.length() && x < pixelToRef.length;++i)
					{
					out.print(f.charAt(i));
					if(pixelToRef[x]!=INSERTION_OPCODE) ++ref;
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
		while(x < pixelToRef.length)
			{
			if(pixelToRef[x]==INSERTION_OPCODE)
				{
				out.print("*");
				++x;
				}
			else
				{
				out.print(refPosToBase.apply(ref-1));
				++ref;
				++x;
				}
			}
		out.println();

		
		for(x=0;x< pixelToRef.length;++x)
			{
			out.print(pixelToRef[x]==INSERTION_OPCODE?"^":" ");
			}
		out.println();
		
		
		for(final String groupName : group2record.keySet())
			{
			out.println(groupName);
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
					while( x < pixelToRef.length &&
							ref <  readRef &&
							testInInterval.test(ref))
						{
						if( pixelToRef[x] != NOT_VISIBLE ) ++ref;
						++x;
						out.print(' ');
						}
										
					int readpos=0;
					
					/* get read base function */
					final Function<Integer,Character> baseAt=new Function<Integer, Character>() {
						@Override
						public Character apply(final Integer readpos) {
							if(readpos<0 || readpos>=rec.getReadLength()) return '?';
							final byte 	c= rec.getReadBases()[readpos];
								
							if(rec.getReadNegativeStrandFlag())
								{
								return (char)Character.toLowerCase(c);
								}
							else
								{
								return (char)Character.toUpperCase(c);
								}
							}
						};
					
					for(final CigarElement ce:rec.getCigar())
						{
						final CigarOperator op = ce.getOperator();
						if(op.equals(CigarOperator.PADDING)) continue;
						
						if(this.showInsertions && op.equals(CigarOperator.I)) {
							int cigarIdx =0;
							while( x < pixelToRef.length &&
									cigarIdx < ce.getLength()
									)
								{
								if(testInInterval.test(readRef)) {
									out.print(baseAt.apply(readpos));
									++x;
									}
								++cigarIdx;
								++readpos;
								}
							continue;
							}

												
						int cigarIdx =0;
						while( x < pixelToRef.length &&
								cigarIdx < ce.getLength() 
								)
							{
							
							//pad before base
							while( x < pixelToRef.length &&
									testInInterval.test(readRef) &&
									(pixelToRef[x] == NOT_VISIBLE || pixelToRef[x]<readRef)
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
											out.print('N');
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
											out.print((char)rec.getReadBases()[readpos]);
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
										colorizer.paper(AnsiColor.YELLOW).print('-');
										++x;
										}
									++readRef;
									break;
									}
								case EQ: case M : case X:
									{
									if(testInInterval.test(readRef)) {
										char refBase =  Character.toUpperCase(refPosToBase.apply(readRef));
										char readBase = Character.toUpperCase(baseAt.apply(readpos));
										if(!(refBase!='N' &&  readBase!='N'  && readBase!=refBase))
											{
											colorizer.pen(AnsiColor.RED);
											//readBase=',';
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
										
										colorizer.print(readBase);
										//out.print(readBase);
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
				while(x<pixelToRef.length)
					{
					out.print(" ");
					++x;
					}
				out.println();
				if(out.checkError()) break;
				}
			if(out.checkError()) break;
			}
		
		}
	
    
	}
