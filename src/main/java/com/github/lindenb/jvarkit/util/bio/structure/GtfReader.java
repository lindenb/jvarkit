/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.bio.structure;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;

/**
 * A GTF/KnownGene Reader
 */
public class GtfReader implements Closeable {
	private static final Logger LOG = Logger.build(GtfReader.class).make();
	public static final String OPT_DESC="A GTF (General Transfer Format) file. See https://www.ensembl.org/info/website/upload/gff.html . "
			+ "Please note that CDS are only detected if a start and stop codons are defined.";
	/** available files extensions for GTF files */
	public static List<String> SUFFIXES = Arrays.asList(".gtf",".gtf.gz");
	
	private final GtfResource resource;
	/** type of  input detected: ucsc knownGene or gtf */
	private enum InputFormat {undefined,knowngene,gtf};
	private InputFormat format = InputFormat.undefined;
	private Function<String,String> contigNameConverter  = S->S;
	private TabixReader tabixReader = null;
	
	public GtfReader(final InputStream in) {
		this.resource = new InputStreamGtfResource(in);
	}
	
	public GtfReader(final Path path) {
		this(path,false);
		}
	
	private GtfReader(final Path path,boolean requireIndex) {
		IOUtil.assertFileIsReadable(path);
		this.resource = new PathGtfResource(path);
		if(requireIndex) {
			try {
				this.tabixReader = new TabixReader(path.toString());
			} catch (final IOException e) {
				throw new RuntimeIOException(e);
				}
			}
		}
	
	public GtfReader(final String uri) {
		if(uri==null) {
			this.resource = new InputStreamGtfResource(System.in);
			}
		else if(IOUtil.isUrl(uri)) {
			this.resource = new RemoteGtfResource(uri);
			}
		else
			{
			final Path path= Paths.get(uri);
			IOUtil.assertFileIsReadable(path);
			this.resource = new PathGtfResource(path);
			}
		}
	
	public void setContigNameConverter(final Function<String, String> contigNameConverter) {
		this.contigNameConverter = contigNameConverter;
		}
	
	/** exon coordinate */
	private static class Coords {
		int start;
		int end;
		}
	
	private class State {
		final Map<String,GeneImpl> id2gene = new HashMap<>();
		final Map<String,TranscriptImpl> id2transcript = new HashMap<>();
		final Map<String,List<Coords>> transcript2exons = new HashMap<>();
		final IntervalTreeMap<Locatable> treemap;
		final GTFCodec codec = new GTFCodec();

		/** @param intervals can be null */
		State(final Collection<Locatable> intervals)
			{
			if(intervals==null)
				{
				this.treemap = null;
				}
			else
				{
				this.treemap = new IntervalTreeMap<>();
				intervals.stream().forEach(R->treemap.put(new Interval(R), R));
				}
			}
		
		private GeneImpl getGene(final String gene_id) {
			GeneImpl g  = this.id2gene.get(gene_id);
			if(g==null) {
				g = new GeneImpl();
				g.gene_id = gene_id;
				this.id2gene.put(gene_id,g);
				}
			return g;
			}
		private TranscriptImpl getTranscript(final String transcript_id) {
			TranscriptImpl g  = this.id2transcript.get(transcript_id);
			if(g==null) {
				g = new TranscriptImpl();
				g.transcript_id = transcript_id;
				this.id2transcript.put(transcript_id,g);
				}
			return g;
			}
		private void visitKg(final String line) {

			final String tokens[] = CharSplitter.TAB.split(line);
			final int binIdx=tokens[2].equals("+") || tokens[2].equals("-")?0:1;
			
			final GeneImpl gene = new GeneImpl();
			gene.gene_id = tokens[binIdx + 0];
			gene.contig = contigNameConverter.apply(tokens[binIdx + 1]);
			if(StringUtils.isBlank(gene.contig)) return;
			final TranscriptImpl transcript =new TranscriptImpl();
			transcript.gene = gene;
			gene.transcripts.add(transcript);
			transcript.txStart =  1+ Integer.parseInt(tokens[binIdx + 3]);
			transcript.txEnd =  Integer.parseInt(tokens[binIdx + 4]);
			transcript.transcript_id = tokens[binIdx + 0];
			transcript.strand = tokens[binIdx + 2].charAt(0);
			transcript.codon_start = transcript.new StartCodonImpl();
			transcript.codon_end = transcript.new StopCodonImpl();

			if(transcript.strand=='+') {
				int p = 1+Integer.parseInt(tokens[binIdx + 5]);
				for(int i=0;i< 3;i++)   transcript.codon_start.pos[i]=p+i;
				p = 1+Integer.parseInt(tokens[binIdx + 6]);
				for(int i=0;i< 3;i++)   transcript.codon_end.pos[i]=p+i;
				}
			else
				{
				int p = 1+Integer.parseInt(tokens[binIdx + 5]);
				for(int i=0;i< 3;i++)   transcript.codon_end.pos[i]=p+i;
				p = 1+Integer.parseInt(tokens[binIdx + 6]);
				for(int i=0;i< 3;i++)   transcript.codon_start.pos[i]=p+i;
				}
			transcript.coding = !tokens[binIdx + 5].equals(tokens[binIdx + 6]);
			gene.start = transcript.txStart;
			gene.end = transcript.txEnd;
	        gene.strand = transcript.strand;
	        
	        final int exonCount=Integer.parseInt(tokens[binIdx + 7]);
	        transcript.exonStarts =  Arrays.stream(CharSplitter.COMMA.split(tokens[binIdx + 8])).limit(exonCount).mapToInt(S->Integer.parseInt(S)).toArray();
	        transcript.exonEnds =  Arrays.stream(CharSplitter.COMMA.split(tokens[binIdx + 9])).limit(exonCount).mapToInt(S->Integer.parseInt(S)).toArray();
	        
	        this.id2gene.put(gene.gene_id, gene);
			}		
		
		private  void visitGtf(final String line) {
			final GTFLine T = this.codec.decode(line);
			if(T==null) return;
			
			final String contig = contigNameConverter.apply(T.getContig());
			if(StringUtils.isBlank(contig)) return;
			
			if(T.getStart()<=0) throw new IllegalArgumentException("Bad start in "+line);
			if(T.getEnd()<=0) throw new IllegalArgumentException("Bad start in "+line);

			
			if(this.treemap!=null && this.treemap.debugGetTree(T.getContig())==null) return;
			
			if(T.getType().equals("gene"))
				{
				final GeneImpl g = this.getGene(getRequiredProperty(T, "gene_id"));
				g.properties.putAll(T.getAttributes());
				g.contig = contig;
				g.start = T.getStart();
				g.end = T.getEnd();
				g.strand = T.getStrand();
				}
			else if(T.getType().equals("transcript"))
				{
				final GeneImpl g = this.getGene(getRequiredProperty(T, "gene_id"));
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				t.properties.putAll(T.getAttributes());
				t.gene = g;
				t.txStart = T.getStart();
				t.txEnd = T.getEnd();
				t.strand = T.getStrand();
				g.transcripts.add(t);
				}
			else if(T.getType().equals("start_codon"))
				{
				// codon start can be spliced of ENST00000429923 : multiple start_codon
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				if(t.codon_start == null) t.codon_start = t.new StartCodonImpl();
				t.codon_start.visit(T);
				t.coding = true;
				}
			else if(T.getType().equals("stop_codon"))
				{
				// codon start can be spliced of ENST00000429923 : multiple start_codon
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				if(t.codon_end == null) t.codon_end = t.new StopCodonImpl();
				t.codon_end.visit(T);
				t.coding = true;
				}
			else if(T.getType().equals("exon"))
				{
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				
				List<Coords> coords = this.transcript2exons.get(t.transcript_id);
				if(coords==null) {
					coords = new ArrayList<>();
					this.transcript2exons.put(t.transcript_id,coords);
					}
				final Coords coord = new Coords();
				coord.start = T.getStart();
				coord.end = T.getEnd();
				coords.add(coord);
				}
			else if(T.getType().equals("five_prime_utr") || 
					T.getType().equals("CDS") || 
					T.getType().equals("three_prime_utr")) {
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				t.saw_cds_flag = true;
				}
			else
				{
				}
			}
		
		List<Gene> finish() {
		     // it happens! ENST00000541351
			// gunzip -c ~/jeter.gtf.gz | grep ENST00000541351 | grep codon
			// 12	havana	stop_codon	10854686	10854688	.	-	0	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "5"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
			
			
			for(final TranscriptImpl tr: this.id2transcript.values())
				{
				/*
				if(tr.codon_start==-1 && tr.codon_end==-1)
					{
					if(tr.saw_cds_flag) {
						// e.g. ENST00000580917.
						//LOG.error("Saw CDS/UTR for "+tr.getId()+". But not codon defined.");
						}
					}
				// happens wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" | gunzip -c | grep ENST00000327956 | grep stop
				else if(tr.codon_start>0 && tr.codon_end==-1)
					{
					
					}
				else if(tr.codon_start>0 && tr.codon_end>0)
					{
					if(tr.isPositiveStrand() && tr.codon_start>tr.codon_end) throw new IllegalStateException("cds bad start/end in "+tr);
					if(tr.isNegativeStrand() && tr.codon_end>tr.codon_start) throw new IllegalStateException("cds bad start/end in "+tr);
					}
				else //codon_end > 0 && codon_start==-1
					{
					//throw new IllegalStateException("only codon_start XOR codon_end defined for "+
					//tr+" codon_start="+tr.codon_start+" codon_end="+tr.codon_end);
					}*/
				}
			
			for(final String transcript_id:this.transcript2exons.keySet())
				{
				final TranscriptImpl tr = this.id2transcript.get(transcript_id);
				if(tr==null) throw new IllegalStateException();
				final List<Coords> coords = this.transcript2exons.get(transcript_id);
				if(coords.isEmpty()) {
					throw new IllegalStateException("No exon defined for "+transcript_id);
					}
				
				Collections.sort(coords,(A,B)->Integer.compare(A.start, B.start));
				tr.exonStarts = coords.stream().mapToInt(E->E.start).toArray();
				tr.exonEnds = coords.stream().mapToInt(E->E.end).toArray();
				}
			
			this.id2gene.values().stream().filter(T->T.start==-1 || T.end==-1 || StringUtils.isBlank(T.contig)).findAny().ifPresent(T->new RuntimeIOException("gene without data : "+T.gene_id));
			this.id2transcript.values().stream().filter(T->T.txStart==-1 || T.txEnd==-1 || !(T.strand=='+' || T.strand=='-')).findAny().ifPresent(T->new RuntimeIOException("transcript without data : "+T.transcript_id));
			
			
			if(this.treemap!=null)
				{
				this.id2gene.values().removeIf(G->!this.treemap.containsOverlapping(G));
				}
			return new ArrayList<>(this.id2gene.values());
			}
		}
	
	private static String getRequiredProperty(final GTFLine line,String prop)  {
		final String s = line.getAttribute(prop);
		if(StringUtils.isBlank(s)) {
			throw new RuntimeIOException("cannot find property \""+prop+"\" in "+line.getLine());
			}
		return s;
		}

	
	public List<Gene> getAllGenes() {
		return fetchGenes(null);
		}
	
	/** not tested */
	private List<Gene> queryGenes(final Locatable interval) {
		if(interval==null) throw new IllegalArgumentException("interval cannot be null");
		if(this.tabixReader==null) throw new IllegalArgumentException("Not a opened as a bgzipped+tabix indexed gtf file.");
		String tabixContig= interval.getContig();
		if(this.contigNameConverter!=null) {
			tabixContig = this.tabixReader.getChromosomes().
				stream().
				map(C->this.contigNameConverter.apply(C)).
				filter(C->C!=null && C.equals(interval.getContig())).
				findFirst().
				orElse(null);
			}
		
		if(StringUtils.isBlank(tabixContig)) return Collections.emptyList();
		
		try {
			 
			final GTFCodec codec = new GTFCodec();
			int start = interval.getStart();
			int end  = interval.getEnd();
			TabixReader.Iterator iter = this.tabixReader.query(tabixContig, start, end);
			for(;;) {
				final String line = iter.next();
				if(line==null) break;
				GTFLine L=codec.decode(line);
				if(L==null) continue;
				if(!L.getType().equals("gene")) continue;
				start = Math.min(start, L.getStart());
				end = Math.max(end, L.getEnd());
				}
			final State state = new State(null);
			iter = this.tabixReader.query(tabixContig, start, end);
			for(;;) {
				final String line = iter.next();
				if(line==null) break;
				state.visitGtf(line);
				}
			return state.finish();
		} catch(final IOException err) {
			throw new RuntimeIOException(err);
		}
	}
	
	
	private List<Gene> fetchGenes(final Collection<Locatable> intervals) {
		final State state = new State(intervals);
		
		try(final BufferedReader br=this.resource.openReader())
			{
			br.lines().
				filter(S->!S.startsWith("#")).
				filter(S->!StringUtils.isBlank(S)).
				forEach(L->{
					if(this.format.equals(InputFormat.undefined)) {
						final String tokens[] = CharSplitter.TAB.split(L);
						if(tokens.length>6 && (tokens[6].equals(".") || tokens[6].equals("+") || tokens[6].equals("-"))) {
							this.format = InputFormat.gtf;
							}
						else
							{
							this.format = InputFormat.knowngene;
							}
						}
					switch(this.format) {
						case gtf: state.visitGtf(L);break;
						case knowngene: state.visitKg(L);break;
						default: throw new IllegalStateException();
					}
				});
			}
		catch (final IOException e) {
			throw new RuntimeIOException(e);
			}
		return  state.finish();
		}
	
	@Override
	public void close()  {
		this.resource.close();
		if(tabixReader!=null) tabixReader.close();
		}
	
	private abstract class GtfResource  implements Closeable {
		public abstract BufferedReader openReader() throws IOException;
		@Override
		public void close() {
			}
		}
	
	private class InputStreamGtfResource  extends GtfResource{
		private final InputStream in;
		short count_open = 0;
		InputStreamGtfResource( final InputStream in) {
			this.in = in;
			}
		@Override
		public BufferedReader openReader() throws IOException {
			if(count_open>0) throw new IOException("Cannot read GTF file from inputsream more than one time.");
			++count_open;
			return IOUtils.openStreamForBufferedReader(this.in);
			}
		@Override
		public void close() {
			CloserUtil.close(in);
			}
		}
	
	private class RemoteGtfResource  extends GtfResource{
		private final String url;
		RemoteGtfResource(final String url) {
			this.url = url;
			}
		@Override
		public BufferedReader openReader() throws IOException {
			return IOUtils.openURIForBufferedReading(this.url);
			}
		}
	
	private class PathGtfResource  extends GtfResource{
		private final Path path;
		PathGtfResource(final Path path) {
			this.path = path;
			}
		@Override
		public BufferedReader openReader() throws IOException {
			return IOUtils.openPathForBufferedReading(this.path);
			}
		}
	

	
	private static class TranscriptImpl implements Transcript
		{
		GeneImpl gene;
		String transcript_id;
		int txStart;
		int txEnd;
		StartCodonImpl codon_start= null;
		StopCodonImpl codon_end = null;
		char strand = '?';
		int exonStarts[]=null;
		int exonEnds[]=null;
		final Map<String,String> properties = new HashMap<String, String>();
		boolean coding = false;
		boolean saw_cds_flag = false;
		
		private abstract class AbstractCodonImpl implements Codon {
			final int pos[]=new int[] {-1,-1,-1};
			
			private void _visit(int loc) {
				int i=0;
				while( i< pos.length) {
					if(pos[i]==-1) break;
					i++;
					}
				if(i>= pos.length) throw new IllegalArgumentException();
				pos[i]=loc;
				}
			
			void visit(Locatable loc) {
				int L=loc.getLengthOnReference();
				if(L>3) throw new IllegalStateException("large codon for "+getName()+" from "+loc);
				for(int i=0;i< L ;i++)
					{
					_visit(loc.getStart()+i);
					}
				Arrays.sort(pos);
				}
			
			private int at(int idx) {
				final int p = pos[idx];
				if(p<1) throw new IllegalStateException("incomplete codon for "+getName()+" "+getTranscript()+ " "+Arrays.toString(this.pos));
				return p;
				}
			
			@Override
			public int getStart() {
				return at(0);
				}
			
			@Override
			public int getMiddle() {
				return at(1);
				}
			
			@Override
			public int getEnd() {
				return at(2);
				}
			
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			}
		
		private class StartCodonImpl extends AbstractCodonImpl {
			@Override
			public String getName() {
				return "Start Codon";
				}
			}
		
		private class StopCodonImpl extends AbstractCodonImpl {
			@Override
			public String getName() {
				return "Stop Codon";
				}
			}
		
		private abstract class AbstractUTR implements UTR {
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			@Override
			public int hashCode() {
				return (getTranscript().getId().hashCode() * 31 +this.getStart())*31 + getEnd();
				}

			
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof AbstractUTR)) return false;
				final AbstractUTR other = AbstractUTR.class.cast(obj);
				return other.getTranscript().equals(this.getTranscript()) &&
						other.getStart() == this.getStart() &&
						other.getEnd() == this.getEnd()
						;
				}
			}
		
		private class UTR5Impl extends AbstractUTR {
			
			
			@Override
			public int getStart() {
				return getTranscript().getTxStart();
				}
			public int getEnd() {
				final Optional<Codon> codon;
				if( isPositiveStrand())
					{
					codon= getTranscript().getCodonStart();
					}
				else 
					{
					codon = getTranscript().getCodonStop();
					}
				
				if(!codon.isPresent()) throw new IllegalStateException("no codon for UTR5' for "+getTranscript().getId());
				return codon.get().getStart() - 1 /* one base before the codon */;
				}
			@Override
			public String getName() {
				return (getTranscript().isNegativeStrand()? "3":"5")+"' UTR of "+getTranscript().getId();
				}
			@Override
			public String toString() {
				return getName();
				}
			}

		private class UTR3Impl extends AbstractUTR {

			@Override
			public int getStart() {
				final Optional<Codon> codon;
				if( isPositiveStrand())
					{
					codon= getTranscript().getCodonStop();
					}
				else 
					{
					codon = getTranscript().getCodonStart();
					}
				if(!codon.isPresent()) throw new IllegalStateException("no codon for UTR3' for "+getTranscript().getId());
				return codon.get().getEnd() + 1 /* one base after the codon */;
				}
			
			public int getEnd() {
				return getTranscript().getTxEnd();
				}
			@Override
			public String getName() {
				return (getTranscript().isNegativeStrand()? "5":"3")+"' UTR of "+getTranscript().getId();
				}
			@Override
			public String toString() {
				return getName();
				}
			}
		private class CdsImpl implements Cds
			{
			private final int exon_index;
			CdsImpl(int exon_index) {
				this.exon_index = exon_index;
				}
			@Override
			public Exon getExon() {
				return getTranscript().getExon(this.exon_index);
				}
			@Override
			public int getStart() {
				return Math.max(
						getTranscript().getExonStart(this.exon_index),
						getTranscript().getLeftmostCodon().get().getStart()
						);
				}
			@Override
			public int getEnd() {
				return Math.min(
						getTranscript().getExonEnd(this.exon_index),
						getTranscript().getRightmostCodon().get().getEnd()
						);
				}
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			
			@Override
			public int hashCode() {
				return TranscriptImpl.this.transcript_id.hashCode() * 31 +this.exon_index;
				}
			
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof CdsImpl)) return false;
				final CdsImpl c = CdsImpl.class.cast(obj);
				return c.getTranscript().equals(this.getTranscript()) &&
						c.exon_index == this.exon_index
						;
				}
			
			
			@Override
			public String getName() {
				return "CDS."+getExon().getName();
				}
			@Override
			public String toString() {
				return getName();
				}
			}
		
		private class ExonImpl implements Exon
			{
			private final int index0;
			ExonImpl(final int index0) {
				this.index0 = index0;
				}
			@Override
			public int getIndex() {
				return this.index0;
				}
			@Override
			public int getStart() {
				return getTranscript().getExonStart(this.index0);
				}
			@Override
			public int getEnd() {
				return getTranscript().getExonEnd(this.index0);
				}
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			@Override
			public int hashCode() {
				return TranscriptImpl.this.transcript_id.hashCode() * 31 +this.index0;
				}

			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof ExonImpl)) return false;
				final ExonImpl ex = ExonImpl.class.cast(obj);
				return ex.getTranscript().equals(this.getTranscript()) &&
						ex.getStart() == this.getStart() &&
						ex.getEnd() == this.getEnd()
						;
				}
			
			@Override
			public String getName() {
				return TranscriptImpl.this.transcript_id+ ".Exon"+
							(isNegativeStrand()?
							getExonCount()-this.index0
							:1+this.index0
							);
				}

			
			@Override
			public String toString() {
				return "Exon "+getContig()+":"+getStart()+"-"+getEnd();
				}
			}
		
		private class IntronImpl implements Intron
			{
			private final int index0;
			IntronImpl(final int index0) {
				this.index0 = index0;
				}
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			@Override
			public int getStart() {
				return TranscriptImpl.this.exonEnds[this.index0] + 1;
				}
			@Override
			public int getEnd() {
				return TranscriptImpl.this.exonStarts[this.index0+1] -1;
				}
			
			@Override
			public int hashCode() {
				return TranscriptImpl.this.transcript_id.hashCode() * 31 +this.index0;
				}
			
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof IntronImpl)) return false;
				final IntronImpl ex = IntronImpl.class.cast(obj);
				return ex.getTranscript().equals(this.getTranscript()) &&
						ex.getStart() == this.getStart() &&
						ex.getEnd() == this.getEnd()
						;
				}
			@Override
			public String toString() {
				return "Intron "+getContig()+":"+getStart()+"-"+getEnd();
				}

			@Override
			public String getName() {
				return TranscriptImpl.this.transcript_id+ ".Intron"+
							(isNegativeStrand()?
							getIntronCount()-this.index0
							:1+this.index0
							);
				}
	
			}

		@Override
		public String getId() {
			return this.transcript_id;
			}
		
		@Override
		public Map<String, String> getProperties() {
			return this.properties;
			}
		
		@Override
		public int hashCode() {
			return this.transcript_id.hashCode();
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof TranscriptImpl)) return false;
			final TranscriptImpl tr = TranscriptImpl.class.cast(obj);
			return tr.transcript_id.equals(this.transcript_id);
			}
		
		@Override
		public String getContig() {
			return getGene().getContig();
			}
		@Override
		public int getTxStart() {
			return txStart;
			}
		@Override
		public int getTxEnd() {
			return txEnd;
			}
		
		
		@Override
		public int getStart() {
			return this.getTxStart();
			}
		@Override
		public int getEnd() {
			return this.getTxEnd();
			}
		@Override
		public Gene getGene() {
			return this.gene;
		}
		@Override
		public Optional<Codon> getCodonStart() {
			return Optional.ofNullable(this.codon_start);
		}
		
		@Override
		public boolean hasCodonStartDefined() {
			return this.codon_start!=null;
			}
		
		@Override
		public boolean hasCodonStopDefined() {
			return this.codon_end!=null;
			}
		
		@Override
		public Optional<Codon> getCodonStop() {
			return Optional.ofNullable(this.codon_end);
		}
		
		@Override
		public int getExonStart(int index0) {
			return this.exonStarts[index0];
			}
		@Override
		public int getExonEnd(int index0) {
			return this.exonEnds[index0];
			}
		
		@Override
		public int getExonCount() {
			return this.exonStarts.length;
			}
		
		@Override
		public Exon getExon(final int index0) {
			return new ExonImpl(index0);
			}
		
		@Override
		public List<Exon> getExons() {
			return IntStream.range(0, this.getExonCount()).
					mapToObj(T->getExon(T)).
					collect(Collectors.toList());
			}
		
		@Override
	    public int getIntronCount() {
	    	return this.getExonCount()-1;
	    	}
		@Override
		public List<Intron> getIntrons() {
			return IntStream.range(0, this.getIntronCount()).
					mapToObj(T->getIntron(T)).
					collect(Collectors.toList());
			}
		
		@Override
		public Intron getIntron(int index0) {
			return new IntronImpl(index0);
			}
		@Override
		public Optional<UTR> getUTR5() {
			if(isNonCoding()) return Optional.empty();
			if(
				(isPositiveStrand() && hasCodonStartDefined()) || 
				(isNegativeStrand() && hasCodonStopDefined())
				) return Optional.of(new UTR5Impl());
			return Optional.empty();
			}
		
		@Override
		public Optional<UTR> getUTR3() {
			if(isNonCoding()) return Optional.empty();
			if(
					(isPositiveStrand() && hasCodonStopDefined()) || 
					(isNegativeStrand() && hasCodonStartDefined())
					) return Optional.of(new UTR3Impl());
			
			return Optional.empty();
			}

		@Override
		public boolean isCoding() {
			return this.coding;
			}
		@Override
		public boolean isNonCoding() {
			return !this.coding;
			}
		
		@Override
		public char getStrand() {
			return this.strand;
			}
		@Override
		public List<Cds> getAllCds() {
			if(!isCoding()) return Collections.emptyList();
			if(!hasCodonStartDefined()) throw new IllegalStateException("transcript has not start codon  defined");
			if(!hasCodonStopDefined()) throw new IllegalStateException("transcript has not end codon defined");
			final List<Cds> L = new ArrayList<>(this.getExonCount());
			for(int i=0;i< this.getExonCount();i++)
				{
				if(this.getExonStart(i) > this.getRightmostCodon().get().getEnd()) break;
				if(this.getExonEnd(i) < this.getLeftmostCodon().get().getStart()) continue;
				L.add(new CdsImpl(i));
				}
			return L;
			}
		
		@Override
		public String toString() {
			return this.transcript_id+" "+getContig()+":"+getStart()+"-"+getEnd();
			}
		}
	
	private static class GeneImpl implements Gene
		{
		String gene_id;
		String contig;
		int start;
		int end;
		char strand;
		final List<Transcript> transcripts= new ArrayList<>();
		final Map<String,String> properties = new HashMap<String, String>();
		
		
		@Override
		public String getId() {
			return this.gene_id;
			}
		
		@Override
		public List<Transcript> getTranscripts() {
			return transcripts;
			}
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		@Override
		public Map<String, String> getProperties() {
			return properties;
			}
		@Override
		public char getStrand() {
			return this.strand;
			}
		@Override
		public int hashCode() {
			return this.gene_id.hashCode();
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof GeneImpl)) return false;
			final GeneImpl tr = GeneImpl.class.cast(obj);
			return tr.gene_id.equals(this.gene_id);
			}

		@Override
		public String toString() {
			return this.gene_id+" "+getContig()+":"+getStart()+"-"+getEnd();
			}
		
		}
	
	}
