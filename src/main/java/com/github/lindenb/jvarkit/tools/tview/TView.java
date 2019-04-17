/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.tview;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.TabixKnownGeneFileReader;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
//import htsjdk.samtools.util.Log;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;


public class TView implements Closeable
	{
	private static final String BLACK_SQUARE="\u25A0";
	private static final Logger LOG = Logger.build(TView.class).make();
	private enum LayoutReads {pileup,name};
	public enum Formatout {tty,plain,html};
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	
	@Parameter(names={"--clip"},description="Show clip")
	private boolean showClip=false;
	private Interval interval=null;
	@Parameter(names={"-R","--reference"},description=Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File referenceFile=null;
	@Parameter(names={"--insert"},description="Show insertions")
	private boolean showInsertions=false;
	@Parameter(names={"--readName"},description="Show read name")
	private boolean showReadName=false;
	@Parameter(names={"--format","--outputformat"},description="Output format")
	private Formatout formatOut=Formatout.tty;
	@Parameter(names={"--hideBases"},description="Hide bases")
	private boolean hideBases=false;
	@Parameter(names={"-V","--variant","--variants","--vcf"},description="Variant file. "+IOUtils.UNROLL_FILE_MESSAGE)
	private File variantFiles = null;
	@Parameter(names="--filter",description=SamRecordFilterFactory.FILTER_DESCRIPTION,splitter=NoSplitter.class)
	private SamRecordFilter samRecordFilter = SamRecordFilterFactory.getDefault();
	@Parameter(names={"-left","--leftmargin"},description="left margin width")
	private int leftMarginWidth=15;
	@Parameter(names={"-maxrows","--maxrowss"},description="maximum number of rows per read group. -1 == all")
	private int maxReadRowPerGroup=-1;
	@Parameter(names={"--hideHomRef"},description="Hide HOM_REF variations")
	private boolean hideHomRef=false;
	@Parameter(names={"--hideNoCall"},description="Hide NO_CALL variations")
	private boolean hideNoCall=false;
	@Parameter(names={"--groupby"},description="Group Reads by. " +SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"--noconsensus"},description="Hide Consensus line")
	private boolean hideConsensus=false;
	@Parameter(names={"--coverage","--depth"},description="Number of rows for coverage (hide:<=0)")
	private int numCoverageRows=10;
	@Parameter(names={"-layout","--layout"},description="Layout reads")
	private LayoutReads layoutReads=LayoutReads.pileup;
	@Parameter(names={"-kg","--knownGenes"},description="Tabix indexed UCSC knownGene File",hidden=true)
	private String knownGeneUri = null;

	
	
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
    	String rgb() {
    		return this.name().toLowerCase();
    		}
    	}
	

	private class Colorizer
		{
		protected AnsiColor _pen = null;
		protected AnsiColor _paper = null;
		protected PrintStream out;
		
		Colorizer(final PrintStream out) {this.out=out;}
		public Colorizer pen(AnsiColor c) {
			this._pen = c;
			return this;
		}
		public Colorizer paper(AnsiColor c) {
			this._paper= c;
			return this;
		}
		public Colorizer clear() {
			this._pen=null;
			this._paper=null;
			return this;
		}
		
		public Colorizer print(final Object o)
			{
			out.print(BLACK_SQUARE.equals(o)?"*":o);
			return this;
			}
		@SuppressWarnings("unused")
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
			final Set<Integer> opcodes = new HashSet<>();
			if(this._paper!=null) opcodes.add(this._paper.paper());
			if(this._pen!=null) opcodes.add(this._pen.pen());
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
			clear();
			return this;
			}
		}

	private class HtmlColorizer extends Colorizer {
		private final XMLStreamWriter wout;
		HtmlColorizer(final PrintStream out) {
			super(out);
			try {
				this.wout=XMLOutputFactory.newFactory().createXMLStreamWriter(out, "UTF-8");
				}
			catch(Throwable err) {
				throw new RuntimeException(err);
				}
			}
	
		@Override
		public Colorizer print(final Object o) {
			try {
				if( _pen!=null ||  _paper!=null) {
					wout.writeStartElement("span");
					final StringBuilder sb=new StringBuilder();
					if(_pen!=null) {
						sb.append("color:").append(_pen.rgb()).append(";");
						}
					if(_paper!=null) {
						sb.append("background-color:").append(_paper.rgb()).append(";");
						}
					wout.writeAttribute("style", sb.toString());
					}
				if(BLACK_SQUARE.equals(o)) {
					wout.writeEntityRef("#x25A0");
					}	
				else
					{
					wout.writeCharacters(String.valueOf(o));
					}
				
				if( _pen!=null ||  _paper!=null) {
					wout.writeEndElement();
					}
				clear();
				}
			catch(XMLStreamException err) 
				{
				throw new RuntimeException(err);
				}
			return this;
			}
		}

	
	
	private static class VcfSource
		{
		File vcfFile;
		VCFFileReader vcfFileReader;
		}
	
	private static class Consensus
		{
		private final Counter<Character> count=new Counter<>();
		void watch(char c) {
			if(Character.isWhitespace(c)) return;
			count.incr(Character.toUpperCase(c));
			}
		
		public char getConsensus() {
			if(count.isEmpty()) return ' ';
			if(count.getCountCategories()!=1) return 'N';
			return  count.getMostFrequent();
			}
		public int getCoverage() {
			return (int)this.count.getTotal();
			}
		}
	
	private static class ConsensusBuilder
		{
		private int x=0;
		private final List<Consensus> bases = new ArrayList<>();
		void add(final char c) {
			while(bases.size() <= this.x) this.bases.add(new Consensus());
			this.bases.get(this.x).watch(c);
			++this.x;
			}
		// end of line builder, reset to begin of line
		public void eol() {
			this.x=0;
			}
		
		}

	

	
	private int distance_between_reads=2;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private final List<SamInputResource> samInputResources=new ArrayList<>();
	private final List<SamReader> samReaders=new ArrayList<>();
	private final List<VcfSource> vcfReaders=new ArrayList<>();
	private TabixKnownGeneFileReader tabixKnownGene =null;
	
	public TView() {
		
	}
	
	public File getReferenceFile()
		{
		return referenceFile;
		}
	
	public void setReferenceFile(final File referenceFile) {
		this.referenceFile = referenceFile;
		}
	
	public void setShowClip(boolean showClip) {
		this.showClip = showClip;
	}
	
	public void setShowReadName(boolean showReadName) {
		this.showReadName = showReadName;
	}
	
	public void setShowInsertions(boolean showInsertions) {
		this.showInsertions = showInsertions;
	}
	
	public void setHideBases(boolean hideBases) {
		this.hideBases = hideBases;
	}
	
	public void setMaxReadRowPerGroup(int maxReadRowPerGroup) {
		this.maxReadRowPerGroup = maxReadRowPerGroup;
	}
	
	private String margin(final Object o)
		{
		if(this.leftMarginWidth<=0) return "";
		StringBuilder str = new StringBuilder(o==null?"":o.toString());
		if(str.length() < this.leftMarginWidth)
			{
			str.append(" ");
			}
		while(str.length() < this.leftMarginWidth)
			{
			str.insert(0, " ");
			}
		while(str.length() > this.leftMarginWidth)
			{
			str.deleteCharAt(0);
			}
		return str.toString();
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

	public void setFormatOut(Formatout formatOut) {
		this.formatOut = formatOut;
		}
	
	public void setSamRecordFilter(SamRecordFilter samRecordFilter) {
		this.samRecordFilter = samRecordFilter;
	}
	
	public int initialize() throws IOException
		{
		if(this.referenceFile!=null) {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.referenceFile);
			}
		if(this.samRecordFilter==null) {
			this.samRecordFilter = SamRecordFilterFactory.ACCEPT_ALL;
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
			final VcfSource vcfSource = new VcfSource();
			LOG.debug("OPEN "+vcfFile);
			vcfSource.vcfFile = vcfFile;
			vcfSource.vcfFileReader = new VCFFileReader(vcfFile,true);
			this.vcfReaders.add(vcfSource);
			}
		if(this.tabixKnownGene!=null) {
			this.tabixKnownGene = new TabixKnownGeneFileReader(this.knownGeneUri);
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

		for(final VcfSource r: this.vcfReaders)
			{
			CloserUtil.close(r.vcfFileReader);
			}
		this.samInputResources.clear();
		this.samReaders.clear();
		this.vcfReaders.clear();
		CloserUtil.close(this.tabixKnownGene);
		this.tabixKnownGene =null;
		}
	
	public SamRecordFilter getRecordFilter() {
		return this.samRecordFilter;
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
		final Colorizer colorizer;
		switch(this.formatOut)
			{
			case html: colorizer = new HtmlColorizer(out);break;
			case tty: colorizer = new AnsiColorizer(out);break;
			case plain: colorizer = new Colorizer(out);break;
			default: throw new IllegalStateException();
			}
	
				
		if(interval==null) 
			{
			LOG.warn("No interval defined");
			return;
			}
		final GenomicSequence contigSequence;
		final Function<Integer, Character> refPosToBase;
		if(indexedFastaSequenceFile!=null)
			{
			final SAMSequenceDictionary dict=SAMSequenceDictionaryExtractor.extractDictionary(referenceFile);	
			if(dict.getSequence(this.interval.getContig())==null)
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
				String group = this.groupBy.getPartion(rec);
				if(group==null || group.isEmpty()) {
					group="undefined_"+this.groupBy.name();
					}
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
		//LOG.debug(genomicpos2insertlen);
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
		/** print interval title */
		out.println(interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
		/** paint base position */
		int ref = this.interval.getStart();
		int x=0;
		out.print(margin("POS:"));
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
		out.print(margin("REF:"));
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

		/*
		for(x=0;x< pixelWidth;++x)
			{
			out.print(insertIsPresentAtX.test(x)?"^":" ");
			}
		out.println();
		*/
		/* loop over samples **/
		for(final String groupName : group2record.keySet())
			{
			if(this.maxReadRowPerGroup==0) continue;
			final ConsensusBuilder consensus = new ConsensusBuilder();
			int y_group=0;
			final List<List<SAMRecord>> rows = new ArrayList<>();
			out.println(margin(""));
			
			switch(this.layoutReads) 
				{
				case name:
						{
						rows.addAll(group2record.get(groupName).
							stream().
							sorted((R1,R2)->R1.getReadName().
							compareTo(R2.getReadName())).
							map(R->Collections.singletonList(R)).
							collect(Collectors.toList())
							);
						break;
						}
				default:
					{
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
					break;
					}
				}
			//each row is only one read, so we need to print the groupName
			if(layoutReads==LayoutReads.name)
					{
					out.print(margin(groupName));
					out.println();
					}
			/* print each row */
			for(final List<SAMRecord> row : rows)
				{
				++y_group;
				boolean print_this_line = (this.maxReadRowPerGroup<0 || y_group <= this.maxReadRowPerGroup);
				
				if(print_this_line) {
					//each row is only one read, print the read name
					if(layoutReads==LayoutReads.name)
						{
						String readName = row.get(0).getReadName();
						if(row.get(0).getReadPairedFlag()) {
							if(row.get(0).getFirstOfPairFlag()) {
								readName+="/1";
								}
							if(row.get(0).getSecondOfPairFlag()) {
								readName+="/2";
								}
							}
						out.print(margin(readName));
						}
					else
						{
						out.print(margin(y_group==1?groupName:""));
						}
					}
				
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
						if(print_this_line) out.print(' ');
						consensus.add(' ');
						}
										
					int readpos=0;
					
					/* get read base function */
					final Function<Integer,Character> baseAt=new Function<Integer, Character>() {
						@Override
						public Character apply(final Integer readpos) {
							final byte readBases[] = rec.getReadBases();
							if( readBases == SAMRecord.NULL_SEQUENCE) return 'N';
							if(readpos<0 || readpos>=rec.getReadLength()) return '?';
							return (char) readBases[readpos];
							}
						};
					
					for(final CigarElement ce:rec.getCigar())
						{
						final CigarOperator op = ce.getOperator();
						if(op.equals(CigarOperator.PADDING)) continue;
						
						/* IN INSERTION, only print if showInsertions is true */
						if(this.showInsertions && op.equals(CigarOperator.I)) {
							int cigarIdx =0;
							while( x < pixelWidth &&
									cigarIdx < ce.getLength()
									)
								{
								if(testInInterval.test(readRef)) {
									final char readbase = baseAt.apply(readpos);
									if(print_this_line) colorizer.paper(AnsiColor.RED).print(readbase);
									consensus.add(readbase);
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
							colorizer.clear();
							//pad before base
							while( x < pixelWidth &&
									testInInterval.test(readRef) &&
									(insertIsPresentAtX.test(x))
									)
								{
								++x;
								if(print_this_line) colorizer.paper(AnsiColor.YELLOW).print("*");
								consensus.add(' ');
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
											if(print_this_line) colorizer.paper(AnsiColor.YELLOW).print('N');
											consensus.add(' ');//CLIPPED base not part of consensus 
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
											final char readBase = baseAt.apply(readpos);
											if(print_this_line) colorizer.paper(AnsiColor.YELLOW).print(readBase);
											consensus.add(' ');//CLIPPED base not part of consensus
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
										if(print_this_line) colorizer.paper(AnsiColor.RED).print('-');
										consensus.add(' ');//deletion not not part of consensus
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
										consensus.add(readBase);
										
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
										if(print_this_line) colorizer.print(readBase);
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
					if(print_this_line) out.print(" ");
					++x;
					}
				if(print_this_line) out.println();
				consensus.eol();
				if(out.checkError()) break;
				}
			if(out.checkError()) break;
			
			if(!this.hideConsensus &&  consensus.bases.stream().anyMatch(C->C.getCoverage()>0))
				{
				out.print(margin(groupName+" CONSENSUS"));
				x=0;
				ref=interval.getStart();
				while(x< consensus.bases.size())
					{
					final char refBase =  Character.toUpperCase(refPosToBase.apply(ref-1));
					final char consensusBase  = consensus.bases.get(x).getConsensus();
					if( Character.isWhitespace(consensusBase) ) {
						//nothing
						}
					else if( refBase!='N' &&  consensusBase!=refBase)
						{
						colorizer.pen(AnsiColor.RED);
						}
					else
						{
						colorizer.pen(base2ansiColor.apply(consensusBase));
						}
					if(!insertIsPresentAtX.test(x)) ++ref;
					colorizer.print(consensusBase);
					++x;
					}
				out.println();
				}
			if(this.numCoverageRows>0)
				{
				int minCov = consensus.bases.stream().mapToInt(C->C.getCoverage()).min().orElse(0);
				final int maxCov = consensus.bases.stream().mapToInt(C->C.getCoverage()).max().orElse(0);
				for(int y=0;maxCov>0 && y<this.numCoverageRows;++y)
					{
					if(minCov==maxCov) minCov--;
					double fract = (maxCov-minCov)/((double)this.numCoverageRows);
					int inverse_y = (this.numCoverageRows-1)-y;
					int d0 = (int)((fract) * inverse_y);
					//int d1 = (int)((fract) * (inverse_y+1));
					
					out.print(margin(y==0 ? groupName+" "+maxCov:(y+1==this.numCoverageRows?String.valueOf(minCov):"")));
					for(x=0;x< consensus.bases.size();++x)
						{
						int depth = consensus.bases.get(x).getCoverage() - minCov;
						colorizer.print(depth>=d0?BLACK_SQUARE:' ');
						}
					out.println();
					}
				}
			}/* end of loop over sample */
		/** known gene section */
		
		if(this.tabixKnownGene!=null && this.indexedFastaSequenceFile!=null)
			{
			final List<KnownGene> genes =  this.tabixKnownGene.getItemsInInterval(this.interval);
			if(!genes.isEmpty()) {
				out.println(this.knownGeneUri);
				for(final KnownGene gene: genes) {
					final KnownGene.CodingRNA codingRna=gene.getCodingRNA(contigSequence);
					final KnownGene.Peptide peptide = codingRna.getPeptide();

					out.print(margin(gene.getName()));
					x = 0;
					int ref0 = this.interval.getStart()-1;
					while(x<pixelWidth)
						{
						if(insertIsPresentAtX.test(x))
							{
							out.print("*");
							++x;
							}
						else
							{
							char pepChar = ' ';
							if(ref0>=gene.getTxStart() && ref0<gene.getTxEnd())
								{
								pepChar=(gene.isPositiveStrand()?'>':'<');
								int pepIdx = peptide.convertGenomicToPeptideCoordinate(ref0);
								if(pepIdx!=-1)
									{
									final String aa3 = GeneticCode.aminoAcidTo3Letters(peptide.charAt(pepIdx));
									final int offset[] = peptide.convertToGenomicCoordinates(pepIdx);
									if(offset!=null && offset.length==3 && aa3!=null && aa3.length()==3) {
										if(offset[0]==ref0 ) pepChar=aa3.charAt(0);
										else if(offset[1]==ref0 ) pepChar=aa3.charAt(1);
										else if(offset[2]==ref0 ) pepChar=aa3.charAt(2);
										else pepChar='?';
										}
									else
										{
										pepChar='?';
										}
									}
								}
							out.print(pepChar);
							++ref0;
							++x;
							}
						}
					
					while(x<pixelWidth)
						{
						out.print(" ");
						++x;
						}
					out.println();
					}
				}
			out.println();
			}
		
		
		/** variant section*/
		if(!this.vcfReaders.isEmpty() && !out.checkError()) {
			final Function<GenotypeType, Character> gTypeToSymbol = new Function<GenotypeType, Character>()
				{
				@Override
				public Character apply(final GenotypeType gt) {
					switch(gt)
						{
						case NO_CALL: return '?';
						case HOM_REF: return '0';
						case HET: return '1';
						case HOM_VAR: return '2';
						case MIXED: return 'm';
						case UNAVAILABLE: return 'u';
						default: return '.';
						}			
					}
				};
			out.println();
			for(final VcfSource r:this.vcfReaders)
				{
				if(out.checkError()) break;
				final VCFHeader header = r.vcfFileReader.getFileHeader();
				final CloseableIterator<VariantContext> iter = r.vcfFileReader.query(this.interval.getContig(), interval.getStart(), interval.getEnd());
				final List<VariantContext> variants = new ArrayList<>();
				while(iter.hasNext())
					{
					variants.add(iter.next());
					}
				iter.close();
				if(variants.isEmpty()) continue;
				
				out.println(r.vcfFile.getPath());
				if(header.hasGenotypingData())
					{
					for(final String sample:header.getSampleNamesInOrder()) {
						if( !variants.stream().
								map(V->V.getGenotype(sample)).
								filter(G->!hideNoCall || (hideNoCall && !G.isNoCall())).
								filter(G->!hideHomRef || (hideHomRef && !G.isHomRef())).
								findAny().isPresent()
								)
							{
							continue;
							}
						out.print(margin(sample));
						ref = this.interval.getStart();
						x=0;
						while(x < pixelWidth)
							{
							if(insertIsPresentAtX.test(x))
								{
								out.print("*");
								++x;
								}
							else
								{
								char refBase= ' ';
								for(final VariantContext ctx:variants)
									{
									if( ctx.getStart() == ref ) {
										final Genotype g = ctx.getGenotype(sample);
										if( g.isNoCall() && this.hideNoCall) continue;
										if( g.isHomRef() && this.hideHomRef) continue;
										refBase = gTypeToSymbol.apply(g.getType());
										break;
										}
									}
								out.print(refBase);
								++ref;
								++x;
								}
							}
						out.println();
						}
					}
				else //no genotype
					{
					for(final VariantContext ctx:variants)
						{
						out.print(
							margin(String.valueOf(ctx.getStart())+":"+
							ctx.getReference().getDisplayString()+"/"+
							ctx.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))
							));
						
						ref = this.interval.getStart();
						x=0;
						while(x < pixelWidth)
							{
							if(insertIsPresentAtX.test(x))
								{
								out.print("*");
								++x;
								}
							else
								{
								out.print(ctx.getStart() == ref?'+':' ');
								++ref;
								++x;
								}
							}
						out.println();
						}
					}
				}
			}
		}
    
	}
