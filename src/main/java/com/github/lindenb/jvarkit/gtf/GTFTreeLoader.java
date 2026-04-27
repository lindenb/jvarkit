package com.github.lindenb.jvarkit.gtf;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

public class GTFTreeLoader implements Consumer<GTFLine> {

	public GTFTreeLoader() {
		// TODO Auto-generated constructor stub
	}

	private final Map<String,GTFTreeNode> gene2node= new HashMap<>(50_000);
	private final Map<String,GTFTreeNode> transcript2node= new HashMap<>(50_000);
	
	void finish() {
		transcript2node.values().stream().forEach(T->T.cleanup());
		gene2node.values().stream().forEach(T->T.cleanup());
		}
	
	@Override
	public void accept(final GTFLine t) {
		if(t==null) return;
		if(t.isGene()) {
			final String gene_id = t.getGeneId();
			if(StringUtils.isBlank(gene_id)) return;
			GTFTreeNode prev = gene2node.get(gene_id);
			if(prev!=null) {
				if(prev.line==null) {
					prev.line = t;
					}
				else
					{
					throw new IllegalArgumentException("Duplicate gene "+gene_id+" "+prev.line+" "+t);
					}
				}
			else
				{
				prev = new  GTFTreeNode();
				prev.line = t;
				this.gene2node.put(gene_id, prev);
				}
			}
		else if(t.isTranscript()) {
			final String transcript_id = t.getTranscriptId();
			if(StringUtils.isBlank(transcript_id)) return;
			final String gene_id = t.getGeneId();
			if(StringUtils.isBlank(gene_id)) return;
			GTFTreeNode gene_node=gene2node.get(gene_id);;
			GTFTreeNode transcript_node = transcript2node.get(gene_id);
			if(gene_node==null) {
				gene_node = new GTFTreeNode();//empty
				this.gene2node.put(gene_id, gene_node);
				}
			
			if(transcript_node!=null) {
				if(transcript_node.line==null) {
					transcript_node.line = t;
					}
				else
					{
					throw new IllegalArgumentException("Duplicate transcript "+transcript_node+" "+transcript_node.line+" "+t);
					}
				}
			else
				{
				transcript_node = new  GTFTreeNode();
				transcript_node.line = t;
				transcript_node.parent = gene_node;
				this.transcript2node.put(transcript_id, transcript_node);
				gene_node.addChild(transcript_node);
				}		
			}
		else if(!StringUtils.isBlank(t.getTranscriptId())) {
			final String transcript_id = t.getTranscriptId();
			final String gene_id = t.getGeneId();
			if(StringUtils.isBlank(gene_id)) return;
			GTFTreeNode gene_node=gene2node.get(gene_id);;
			GTFTreeNode transcript_node = transcript2node.get(gene_id);
			if(gene_node==null) {
				gene_node = new GTFTreeNode();//empty
				this.gene2node.put(gene_id, gene_node);
				}
			
			if(transcript_node==null) {
				transcript_node = new GTFTreeNode();//empty
				this.transcript2node.put(transcript_id, transcript_node);
				gene_node.addChild(transcript_node);
				}
			GTFTreeNode feat  = new GTFTreeNode();
			feat.line = t;
			transcript_node.addChild(feat);
			}
		}

	public static List<GTFTreeNode> slurpGenes(final Iterator<GTFLine> iter) throws IOException {
		final GTFTreeLoader loader=new GTFTreeLoader();
		while(iter.hasNext()) {
			loader.accept(iter.next());
			}
		loader.finish();
		return loader.gene2node.values().stream().filter(G->G.line!=null).collect(Collectors.toList());
		}
	public static List<GTFTreeNode> slurpGenes(final Path path) throws IOException {
		try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
			return slurpGenes(br);
			}
		}
	public static List<GTFTreeNode> slurpGenes(final BufferedReader br) throws IOException {
		return slurpGenes(new GTFCodec(),br);
		}

	public static List<GTFTreeNode> slurpGenes(final GTFCodec codec,final BufferedReader br) throws IOException {
		final GTFTreeLoader loader=new GTFTreeLoader();
		final GTFCodec codec = new GTFCodec();
		String line;
		while((line=br.readLine())!=null) {
			if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
			GTFLine rec = codec.decode(line);
			if(rec==null) continue;
			loader.accept(rec);
			}
		return loader.gene2node.values().stream().filter(G->G.line!=null).collect(Collectors.toList());
		}
}
