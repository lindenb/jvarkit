package com.github.lindenb.jvarkit.util.bio.structure;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.structure.AbstractGxxReader.TranscriptImpl.ExonImpl;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.readers.LineIterator;

public class GFF3Reader extends AbstractGxxReader {
	private static final Logger LOG = Logger.build(GFF3Reader.class).make();

	private final Path path;
	
	public GFF3Reader(final Path path) {
		this.path = path;
		}
	
	private String getAttribute(final Gff3Feature feat,String key)  {
		final List<String> L = feat.getAttribute(key);
		if(L.size()==1) return L.get(0);
		LOG.warn("multiple key "+key+" for "+feat);
		return null;
 	}
	private String getRequiredAttribute(final Gff3Feature feat,String key)  {
		final String val = getAttribute(feat, key);
		if(StringUtils.isBlank(val)) throw new RuntimeIOException("undefined "+key+" for "+feat);
		return val;
 		}
	
	private void parseFeat(Gff3Feature feat) {
		if(feat==null) return;
		if(feat.getType().equals("gene")) {
			parseGene(feat);
			return;
		}
		for(final Gff3Feature c:feat.getChildren()) {
			parseFeat(c);
		}
	}
	
	private void parseGene(Gff3Feature feat) {
		final GeneImpl gene = new GeneImpl();
		gene.gene_id = getRequiredAttribute(feat,"gene_id");
		gene.properties.put("gene_name",getRequiredAttribute(feat,"gene_name"));
		gene.contig = feat.getContig();
		gene.start = feat.getStart();
		gene.end = feat.getEnd();
		gene.strand = feat.getStrand().encodeAsChar();
		for(final Gff3Feature c:feat.getChildren()) {
			if(feat.getType().equals("transcript")) {
				gene.transcripts.add(parseTranscript(gene,c));
			}
		}
	}
	
	private TranscriptImpl parseTranscript(final GeneImpl gene,Gff3Feature feat) {
		final TranscriptImpl transcript = new TranscriptImpl();
		transcript.gene = gene;
		transcript.strand = feat.getStrand().encodeAsChar();
		
		final List<Gff3Feature> exonFeat = feat.getChildren().stream().
				filter(F-> F.getType().equals("exon")).
				sorted((A,B)->Integer.compare(A.getStart(), B.getStart())).
				collect(Collectors.toList());
		
		transcript.exonStarts = exonFeat.stream().mapToInt(E->E.getStart()).toArray();
		transcript.exonEnds = exonFeat.stream().mapToInt(E->E.getEnd()).toArray();
		
		return transcript;
		}
	
	
	
	private void load() throws IOException {
		final Gff3Codec codec = new Gff3Codec(DecodeDepth.DEEP);
		try(BufferedInputStream bufferedInputStream = IOUtil.toBufferedStream(IOUtil.openFileForReading(this.path))) {
			final LineIterator lr = codec.makeSourceFromStream(bufferedInputStream);
			codec.readHeader(lr);
			while(!codec.isDone(lr)) {
				final Gff3Feature feat = codec.decode(lr);
				parseFeat(feat);
				}
			
			codec.close(lr);
			}
		}
	
	@Override
	public void close()  {
		
	}

}
