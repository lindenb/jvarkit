package com.github.lindenb.jvarkit.util.bio.structure;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

public class GftReader implements Closeable {
	private final GtfResource resource;
	private enum InputFormat {undefined,knowngene,gtf};
	private InputFormat format = InputFormat.undefined;
	GftReader(final String uri) {
		if(IOUtil.isUrl(uri)) {
			this.resource = new RemoteGtfResource(uri);
			}
		else
			{
			final Path path= Paths.get(uri);
			IOUtil.assertFileIsReadable(path);
			this.resource = new PathGtfResource(path);
			}
		}
	
	private static class State {
		final Map<String,GeneImpl> id2gene = new HashMap<>();
		final Map<String,TranscriptImpl> id2transcript = new HashMap<>();
		final IntervalTreeMap<Locatable> treemap;
		final GTFCodec codec = GTFCodec.createGtfCodec();

		
		State(Collection<Locatable> intervals)
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
		
		GeneImpl getGene(final String gene_id) {
			GeneImpl g  = this.id2gene.get(gene_id);
			if(g==null) {
				g = new GeneImpl();
				g.gene_id = gene_id;
				this.id2gene.put(gene_id,g);
				}
			return g;
			}
		TranscriptImpl getTranscript(final String transcript_id) {
			TranscriptImpl g  = this.id2transcript.get(transcript_id);
			if(g==null) {
				g = new TranscriptImpl();
				g.transcript_id = transcript_id;
				this.id2transcript.put(transcript_id,g);
				}
			return g;
			}
		void visitKg(final String line) {
			final String tokens[] = CharSplitter.TAB.split(line);
			final int binIdx=tokens[2].equals("+") || tokens[2].equals("-")?0:1;
			
			final GeneImpl gene = new GeneImpl();
			gene.gene_id = tokens[binIdx + 0];
			gene.contig = tokens[binIdx + 1];
			final TranscriptImpl transcript =new TranscriptImpl();
			transcript.gene = gene;
			gene.transcripts.add(transcript);
			transcript.txStart =  1+ Integer.parseInt(tokens[binIdx + 3]);
			transcript.txEnd =  Integer.parseInt(tokens[binIdx + 4]);
			transcript.transcript_id = tokens[binIdx + 0];
			transcript.strand = tokens[binIdx + 2].charAt(0);
	        transcript.cdsStart= 1+Integer.parseInt(tokens[binIdx + 5]);
	        transcript.cdsEnd= Integer.parseInt(tokens[binIdx + 6]);

			gene.start = transcript.txStart;
			gene.end = transcript.txEnd;
	        gene.strand = transcript.strand;
	        
	        final int exonCount=Integer.parseInt(tokens[binIdx + 7]);
	        int exonStarts[] =  Arrays.stream(CharSplitter.COMMA.split(tokens[binIdx + 8])).mapToInt(S->Integer.parseInt(S)).toArray();
	        int exonEnds[] =  Arrays.stream(CharSplitter.COMMA.split(tokens[binIdx + 9])).mapToInt(S->Integer.parseInt(S)).toArray();
	        for(int i=0;i<exonCount;i++ )
	        	{
	        	ExonImpl ex= new ExonImpl();
	        	ex.transcript = transcript;
	        	transcript.exons.add(ex);
	        	ex.start = 1 + exonStarts[i];
	        	ex.end =  exonEnds[i];
	        	}
	        this.id2gene.put(gene.gene_id, gene);
			}		
		
		void visitGtf(final String line) {
			
			final GTFLine T = this.codec.decode(line);
			if(T==null) return;
			
			if(this.treemap!=null && this.treemap.debugGetTree(T.getContig())==null) return;
			
			if(T.getType().equals("gene"))
				{
				final GeneImpl g = this.getGene(getRequiredProperty(T, "gene_id"));
				g.start = T.getStart();
				g.end = T.getEnd();
				g.strand = T.getStrand();
				}
			else if(T.getType().equals("transcript"))
				{
				final GeneImpl g = this.getGene(getRequiredProperty(T, "gene_id"));
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				t.gene = g;
				t.txStart = T.getStart();
				t.txEnd = T.getEnd();
				t.strand = T.getStrand();
				}
			else if(T.getType().equals("start_codon"))
				{
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				t.cdsStart = T.getStart();
				}
			else if(T.getType().equals("stop_codon"))
				{
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				t.cdsEnd = T.getStart();
				}
			else if(T.getType().equals("exon"))
				{
				final TranscriptImpl t = this.getTranscript(getRequiredProperty(T, "transcript_id"));
				final ExonImpl ex = new ExonImpl();
				ex.transcript = t;
				ex.start = T.getStart();
				ex.end = T.getEnd();
				t.exons.add(ex);
				}
			}
		List<Gene> finish() {
			this.id2gene.values().stream().filter(T->T.start==-1 || T.end==-1 || StringUtils.isBlank(T.contig)).findAny().ifPresent(T->new RuntimeIOException("gene without data : "+T.gene_id));
			this.id2transcript.values().stream().filter(T->T.txStart==-1 || T.txEnd==-1 || !(T.strand=='+' || T.strand=='-')).findAny().ifPresent(T->new RuntimeIOException("transcript without data : "+T.transcript_id));
			
			if(this.treemap!=null)
				{
				this.id2gene.values().removeIf(G->!this.treemap.containsOverlapping(G));
				}
			
			// sort exons for each transcript
			this.id2gene.values().
				stream().
				flatMap(G->G.transcripts.stream()).
				forEach(T->{
					TranscriptImpl tr= TranscriptImpl.class.cast(T);
					Collections.sort(
						tr.exons,
						(A,B)->Integer.compare(A.getStart(), B.getStart()));
					for(int i=0;i< tr.exons.size();i++)
						{
						ExonImpl.class.cast(tr.exons.get(i)).exon_index = i;
						}
					
				});
			
			
			
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
		}
	
	private abstract class GtfResource  implements Closeable {
		public abstract BufferedReader openReader() throws IOException;
		@Override
		public void close() {
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
	
	private static class TranscriptIntervalImpl implements TranscriptInterval
		{
		TranscriptImpl transcript;
		int start;
		int end;

		@Override
		public Transcript getTranscript() {
			return null;
			}
		@Override
		public String getContig() {
			return getTranscript().getContig();
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		}
	
	private static class ExonImpl extends TranscriptIntervalImpl implements Exon
		{
		int exon_index = 0;
		}

	
	private static class TranscriptImpl implements Transcript
		{
		GeneImpl gene;
		String transcript_id;
		List<Exon> exons= new ArrayList<>();
		int txStart;
		int txEnd;
		int cdsStart;
		int cdsEnd;
		char strand;
		@Override
		public String getContig() {
			return getGene().getContig();
			}
		public int getTxStart() {
			return txStart;
			}
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
		public List<Exon> getExons() {
			return exons;
			}
		@Override
		public char getStrand() {
			return this.strand;
			}
		}
	private static class GeneImpl implements Gene
		{
		String gene_id;
		String contig;
		int start;
		int end;
		char strand;
		List<Transcript> transcripts= new ArrayList<>();
		Map<String,String> properties = new HashMap<String, String>();
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
		
		}
	
	}
