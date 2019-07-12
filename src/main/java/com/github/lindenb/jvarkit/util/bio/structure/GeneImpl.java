package com.github.lindenb.jvarkit.util.bio.structure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

class GeneImpl implements Gene {

	private String contig;
	private final List<Transcript> transcripts = new ArrayList<>();
	private final Map<String,String> properties = new HashMap<>();
	
	class TranscriptImpl
		implements Transcript
		{
		private int txStart1;
		private int txEnd1;
		private final List<Exon> exons = new ArrayList<>();
		private final Map<String,String> properties = new HashMap<>();


		@Override
		public Gene getGene() {
			return GeneImpl.this;
			}
		@Override
		public String getContig() {
			return getGene().getContig();
			}
		@Override
		public int getStart() {
			return txStart1;
			}
		@Override
		public int getEnd() {
			return txEnd1;
			}
		

		}
	
	
	
	@Override
	public String getContig() {
		return contig;
		}
	
	public List<Transcript> getTranscripts() {
		return transcripts;
		}
	@Override
	public Map<String, String> getProperties() {
		return this.properties;
		}
	}
