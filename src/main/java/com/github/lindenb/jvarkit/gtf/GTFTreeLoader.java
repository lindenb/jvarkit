/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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

/** LOAD a GTFTreeLoader tree */
public class GTFTreeLoader implements Consumer<GTFLine> {

	public GTFTreeLoader() {
		// TODO Auto-generated constructor stub
		}

	private final Map<String,GTFTreeNode> gene2node= new HashMap<>(50_000);
	private final Map<String,GTFTreeNode> transcript2node= new HashMap<>(50_000);
	
	private void clear() {
		gene2node.clear();
		transcript2node.clear();
		}
	
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
			GTFTreeNode prev = this.gene2node.get(gene_id);
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

	private List<GTFTreeNode> _getGenes() {
		this.finish();
		final List<GTFTreeNode> L= this.gene2node.values()
				.stream().filter(G->G.line!=null)
				.collect(Collectors.toList());
		this.clear();
		return L;
		}
	
	public  List<GTFTreeNode> slurpGenes(final Iterator<GTFLine> iter) throws IOException {
		while(iter.hasNext()) {
			this.accept(iter.next());
			}
		return _getGenes();
		}
	public List<GTFTreeNode> slurpGenes(final Path path) throws IOException {
		try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
			return slurpGenes(br);
			}
		}
	public List<GTFTreeNode> slurpGenes(final BufferedReader br) throws IOException {
		return slurpGenes(new GTFCodec(),br);
		}

	public  List<GTFTreeNode> slurpGenes(final GTFCodec codec,final BufferedReader br) throws IOException {
		String line;
		while((line=br.readLine())!=null) {
			if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
			final GTFLine rec = codec.decode(line);
			if(rec==null) continue;
			this.accept(rec);
			}
		return _getGenes();
		}
	
	
}
