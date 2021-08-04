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

import java.io.BufferedInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.readers.LineIterator;

public class GFF3Reader extends AbstractGxxReader {
	private static final Logger LOG = Logger.build(GFF3Reader.class).make();

	private final Path path;
	
	private static class PosIter extends AbstractIterator<Integer> {
		final List<Gff3Feature> features;
		final char strand;
		int feat_index;
		int pos;
		PosIter(List<Gff3Feature> features,char strand) {
			this.features = features;
			this.strand = strand;
			if(strand=='+') {
				feat_index=0;
				pos = features.get(0).getStart() - 1;
				}
			else
				{
				feat_index=features.size()-1;
				pos = features.get(feat_index).getEnd() + 1;
				}
			}
		@Override
		protected Integer advance() {
			int ret = pos;
			if(ret<0) return null;
			if(strand=='+') {
				pos++;
				if(pos>features.get(feat_index).getEnd()) {
					feat_index++;
					if(feat_index>=this.features.size()) {
						pos=-1;
						return null;
						}
					pos = features.get(feat_index).getStart();
					}
				}
			else {
				pos--;
				if(pos<features.get(feat_index).getStart()) {
					feat_index--;
					if(feat_index<0) {
						pos=-1;
						return null;
						}
					pos = features.get(feat_index).getEnd();
					}
				}
			return ret;
			}
		}
	
	
	public GFF3Reader(final Path path) {
		this.path = path;
		}
	
	private String getAttribute(final Gff3Feature feat,String key)  {
		final List<String> L = feat.getAttribute(key);
		if(L.isEmpty()) return null;
		if(L.size()==1) return L.get(0);
		LOG.warn("multiple key "+key+" for "+feat.getContig()+":"+feat.getStart()+"-"+feat.getEnd()+" "+String.join(" ", L));
		return null;
 	}
	private String getRequiredAttribute(final Gff3Feature feat,String key)  {
		final String val = getAttribute(feat, key);
		if(StringUtils.isBlank(val)) throw new RuntimeIOException(
				"undefined "+key+" for "+feat.getAttributes());
		return val;
 		}
	
	private void parseFeat(Gff3Feature feat,final List<Gene> genes) {
		if(feat==null) return;
		if(feat.getType().equals("gene")) {
			genes.add(parseGene(feat));
			}
		else
			{
			for(final Gff3Feature c:feat.getChildren()) {
				parseFeat(c,genes);
				}
			}
		return;
		}

	private Gene parseGene(Gff3Feature feat) {
		final GeneImpl gene = new GeneImpl();
		gene.gene_id = getRequiredAttribute(feat,"gene_id");
		String geneName = getAttribute(feat, "gene_name");
		if(StringUtils.isBlank(geneName)) geneName = getAttribute(feat, "Name");
		if(StringUtils.isBlank(geneName)) throw new RuntimeIOException("No @Name or @gene_name for "+feat.getAttributes());
		gene.properties.put("gene_name",geneName);
		gene.contig = feat.getContig();
		gene.start = feat.getStart();
		gene.end = feat.getEnd();
		gene.strand = feat.getStrand().encodeAsChar();
		for(final Gff3Feature c:feat.getChildren()) {
			final TranscriptImpl tr = parseTranscript(gene,c);
			if(tr!=null) {
				gene.transcripts.add(tr);
				}
			else
				{
				LOG.warn("Cannot parse "+c+" of gene "+ gene.gene_id);
				}
			}
		return gene;
		}

	private TranscriptImpl parseTranscript(final GeneImpl gene,Gff3Feature feat) {
		
		if(feat!=null && feat.getAttributes()!=null && 
				feat.getAttributes().entrySet().stream().filter(KV->KV.getValue()!=null).flatMap(KV->KV.getValue().stream()).anyMatch(S->S!=null && S.equals("ENSG00000235954"))) {
				LOG.info("OFUND "+feat.getType()+" "+feat.getStart()+"-"+feat.getEnd()+" "+feat.getAttributes());
				}

		
		
		final TranscriptImpl transcript = new TranscriptImpl();
		for(final String k:feat.getAttributes().keySet()) {
			List<String> values = feat.getAttribute(k);
			if(values.size()>0) transcript.properties.put(k, values.get(0));
		}
		transcript.transcript_id = getRequiredAttribute(feat, "transcript_id");
		transcript.gene = gene;
		transcript.strand = feat.getStrand().encodeAsChar();
		
		final List<Gff3Feature> exonFeat = feat.getChildren().stream().
				filter(F-> F.getType().equals("exon")).
				sorted((A,B)->Integer.compare(A.getStart(), B.getStart())).
				collect(Collectors.toList());
		
		if(exonFeat.isEmpty()) {
			LOG.warn("no exon for "+feat);
			return null;
		}
		
		transcript.exonStarts = exonFeat.stream().mapToInt(E->E.getStart()).toArray();
		transcript.exonEnds = exonFeat.stream().mapToInt(E->E.getEnd()).toArray();
		transcript.txStart = feat.getStart();
		transcript.txEnd = feat.getEnd();
		
		
		final List<Gff3Feature> cdsFeat = feat.getChildren().stream().
				filter(F-> F.getType().equals("CDS")).
				sorted((A,B)->Integer.compare(A.getStart(), B.getStart())).
				collect(Collectors.toList());

		
		if( !cdsFeat.isEmpty()) {
			transcript.coding = true;
			transcript.saw_cds_flag=true;
			transcript.codon_start = transcript.new StartCodonImpl();
			transcript.codon_end = transcript.new StopCodonImpl();
			
			if(transcript.isPositiveStrand()) {
				PosIter iter1 = new PosIter(cdsFeat, '+');
				int x=0;
				while(iter1.hasNext() && x<3) {
					transcript.codon_start.pos[x]=iter1.next();
					x++;
					}
				iter1 = new PosIter(cdsFeat, '-');
				x=0;
				while(iter1.hasNext() && x<3) {
					transcript.codon_end.pos[2-x]=iter1.next();
					x++;
					}
				}
			else
				{
				PosIter iter1 = new PosIter(cdsFeat, '-');
				int x=0;
				while(iter1.hasNext() && x<3) {
					transcript.codon_start.pos[2-x]=iter1.next();
					x++;
					}
				iter1 = new PosIter(cdsFeat, '+');
				x = 0;
				while(iter1.hasNext() && x<3) {
					transcript.codon_end.pos[x]=iter1.next();
					x++;
					}
				}
			} 
		else
			{
			transcript.saw_cds_flag=true;
			transcript.coding = false;
			}
		
		return transcript;
		}
	
	
	private List<Gene> load() throws IOException {
		final List<Gene> genes = new ArrayList<>();
		final Gff3Codec codec = new Gff3Codec(DecodeDepth.DEEP);
		try(BufferedInputStream bufferedInputStream = IOUtil.toBufferedStream(IOUtil.openFileForReading(this.path))) {
			final LineIterator lr = codec.makeSourceFromStream(bufferedInputStream);
			codec.readHeader(lr);
			
			while(!codec.isDone(lr)) {
				final Gff3Feature feat = codec.decode(lr);
				
				parseFeat(feat,genes);
				}
			
			codec.close(lr);
			}
		return genes;
		}
	
	@Override
	public void close()  {
		
	}

	public static void main(String[] args) {
		try {
			GFF3Reader r=new GFF3Reader(Paths.get("/home/lindenb/jeter.gff"));
			 List<Gene> genes=r.load();
			r.close();
			/*
			genes.stream().
				flatMap(G->G.getTranscripts().stream()).forEach(T->System.err.println(T));
				*/
			LOG.debug("genes N="+genes.stream().flatMap(G->G.getTranscripts().stream()).count());
		} catch(Throwable err) {
			LOG.error(err);
		}
	}
	
}
