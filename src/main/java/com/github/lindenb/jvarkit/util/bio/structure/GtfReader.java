/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.util.function.Function;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.log.Logger;

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
public class GtfReader extends AbstractGxxReader {
	private static final Logger LOG = Logger.build(GtfReader.class).make();
	public static final String OPT_DESC="A GTF (General Transfer Format) file. See https://www.ensembl.org/info/website/upload/gff.html . "
			+ "Please note that CDS are only detected if a start and stop codons are defined.";
	/** available files extensions for GTF files */
	public static List<String> SUFFIXES = Arrays.asList(".gtf",".gtf.gz");
	
	private final GtfResource resource;
	/** type of  input detected: ucsc knownGene or gtf */
	private enum InputFormat {undefined,knowngene,gtf,gff};
	private InputFormat format = InputFormat.undefined;
	private Function<String,String> contigNameConverter  = S->S;
	private TabixReader tabixReader = null;
	private static final boolean SUPPORTS_GFF = false;
	
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
		
		private Map<String,String> extractGffProperties(final String str) {
			final Map<String,String> properties = new HashMap<>();
			int i=0;
			final StringBuilder key = new StringBuilder();
			final StringBuilder value = new StringBuilder();
			while(i<str.length()) {
				key.setLength(0);
				value.setLength(0);
				while(i<str.length() && str.charAt(i)!='=') {
					key.append(str.charAt(i));
					i++;
					}
				if(key.length()==0) throw new IllegalArgumentException("Empty key in "+str);
				if(i>=str.length()) throw new IllegalArgumentException("No value for "+key+" in "+str);
				if(str.charAt(i)!='=')  throw new IllegalArgumentException("Expected '=' after key="+key+" in "+str);
				i++;//skip '='
				while(i<str.length() && str.charAt(i)!=';') {
					value.append(str.charAt(i));
					i++;
					}
				final String k = key.toString();
				if(properties.containsKey(k)) throw new IllegalArgumentException("duplicate key="+k+" in "+str);
				properties.put(k, value.toString());
				}
			return properties;
			}
		
		private  void visitGff(final String line) {
			final String tokens[] = CharSplitter.TAB.split(line);
			final String source = tokens[1];
			final String type = tokens[2];
			if(type.equals("chromosome")) return;
			final String contig = contigNameConverter.apply(tokens[0]);
		
			if(StringUtils.isBlank(contig)) return;
			final int start = Integer.parseInt(tokens[3]);
			final int end = Integer.parseInt(tokens[4]);
			if(start<=0) throw new IllegalArgumentException("Bad start in "+line);
			if(end<=0) throw new IllegalArgumentException("Bad start in "+line);

			final Map<String,String> properties = extractGffProperties(tokens[8]);
			
			if(type.equals("gene"))
				{
				String geneid = getRequiredProperty(line,properties, "ID");
				if(geneid.startsWith("gene:")) geneid=geneid.substring(5);
				final GeneImpl g = this.getGene(geneid);
				g.properties.putAll(properties);
				g.contig = contig;
				g.start = start;
				g.end = end;
				g.strand = tokens[6].charAt(0);
				}
			else if(type.equals("transcript") || type.equals("processed_transcript"))
				{
				String geneid = getRequiredProperty(line,properties, "Parent");
				if(geneid.startsWith("gene:")) geneid=geneid.substring(5);

				
				final GeneImpl g = this.getGene(geneid);
				
				String transcriptid = getRequiredProperty(line,properties, "ID");
				if(transcriptid.startsWith("transcript:")) transcriptid=transcriptid.substring(11);

				
				final TranscriptImpl t = this.getTranscript(transcriptid);
				t.properties.putAll(properties);
				t.gene = g;
				t.txStart = start;
				t.txEnd = end;
				t.strand =tokens[6].charAt(0);
				g.transcripts.add(t);
				}
			else if(type.equals("exon"))
				{
				String transcriptid = getRequiredProperty(line,properties, "Parent");
				if(transcriptid.startsWith("transcript:")) transcriptid=transcriptid.substring(11);

				final TranscriptImpl t = this.getTranscript(transcriptid);
				
				List<Coords> coords = this.transcript2exons.get(t.transcript_id);
				if(coords==null) {
					coords = new ArrayList<>();
					this.transcript2exons.put(t.transcript_id,coords);
					}
				final Coords coord = new Coords();
				coord.start = start;
				coord.end = end;
				coords.add(coord);
				}
			throw new UnsupportedOperationException("input looks like gff but gff is currently not supported");
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
	private static String getRequiredProperty(final String line,final Map<String,String> map,String prop)  {
		final String s = map.get(prop);
		if(StringUtils.isBlank(s)) {
			throw new RuntimeIOException("cannot find property \""+prop+"\" in "+line);
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
				filter(S->!StringUtils.isBlank(S)).
				forEach(L->{
					if(L.startsWith("#")) {
						if(this.format.equals(InputFormat.undefined)) {
							if(L.startsWith("##gff-version") && SUPPORTS_GFF) {
								this.format = InputFormat.gff;
								}
							}
						return;
					}

					
					if(this.format.equals(InputFormat.undefined)) {
						final String tokens[] = CharSplitter.TAB.split(L);
						if(tokens.length>6 && (tokens[6].equals(".") || tokens[6].equals("+") || tokens[6].equals("-"))) {
							if(SUPPORTS_GFF && tokens.length>8 && Pattern.compile("^[A-Za-z_][^ \t\"]*=").matcher(tokens[9]).find()) {
								this.format = InputFormat.gff;
								}
							else
								{
								this.format = InputFormat.gtf;
								}
							}
						else
							{
							this.format = InputFormat.knowngene;
							}
						}
					switch(this.format) {
						case gtf: state.visitGtf(L);break;
						case gff: state.visitGff(L);break;
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
	

	
	}
