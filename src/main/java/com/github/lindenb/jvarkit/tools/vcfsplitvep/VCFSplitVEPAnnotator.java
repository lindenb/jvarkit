/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfsplitvep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCFSplitVEPAnnotator implements VariantAnnotator {
	private static final Logger LOG = Logger.of(VCFSplitVEPAnnotator.class);
	private String prefix = "VEP_";
	private final List<Tag> tags = new ArrayList<>();
	private VepPredictionParser parser=null;
	private enum Aggregate {none,min,max,uniq,random,first};
	
	private static class Tag {
		String tag;
		VCFInfoHeaderLine info;
		Function<Stream<String>,Object> converter= S->S.collect(Collectors.toList());
		}
	
	public VCFSplitVEPAnnotator() {
		}

	public void setTags(final String tags) {
		final Random random = new Random(System.currentTimeMillis());
		for(String tagstr: CharSplitter.COMMA.split(tags)) {
			if(StringUtils.isBlank(tagstr)) continue;
			
			
			final Tag t = new Tag();
			final String[] tokens = CharSplitter.COLON.split(tagstr);
			t.tag = tokens[0];
			String suffix="";
			VCFHeaderLineType type = VCFHeaderLineType.String;
			VCFHeaderLineCount count = VCFHeaderLineCount.UNBOUNDED;
			
			if(tokens.length>1) {
				final Function<String,Object> converter;
				type = VCFHeaderLineType.valueOf(tokens[1]);
				switch(type) {
					case Character: converter = S->{if(S.length()!=1) throw new IllegalArgumentException("Cannot cast "+S+" to character");return S;};break;
					case Flag: converter = S->(isEmpty(S)?false:true);break;
					case Integer: converter= S->Integer.parseInt(S); break;
					case Float: converter= S->Double.parseDouble(S); break;
					case String:  converter= S->S; break;
					default: throw new IllegalArgumentException();
					}
				
				if(tokens.length>2) {
					if(tokens.length>3) {
						throw new IllegalArgumentException("bad syntax "+tagstr);
						}
					final Aggregate aggregate = Aggregate.valueOf(tokens[2]);
					// set suffix
					switch(aggregate) {
						case first:  case max : case min : case uniq:  case random: suffix="_"+aggregate.name();break;
						case none: break;
						default: throw new IllegalArgumentException();
						}
					// set VCFHeaderLineCount
					switch(aggregate) {
						case first:  case max : case min : case random:count=VCFHeaderLineCount.INTEGER;break;
						case uniq: break;
						case none: break;
						default: throw new IllegalArgumentException();
						}
					// aggregate
					switch(aggregate) {
						case first:  t.converter = STREAM->STREAM.filter(S->!isEmpty(S)).map(converter).findFirst().orElse(null); break;
						case max : t.converter = STREAM->STREAM.filter(S->!isEmpty(S)).map(converter).sorted().reduce((A,B)->B).orElse(null); break;
						case random : t.converter = STREAM->{
								final List<Object> L=STREAM.filter(S->!isEmpty(S)).map(converter).collect(Collectors.toCollection(ArrayList::new)); 
								Collections.shuffle(L, random);
								return L.isEmpty()?null:L.get(0);
								};
							break;
							
						case min : t.converter = STREAM->STREAM.filter(S->!isEmpty(S)).map(converter).sorted().findFirst().orElse(null); break;
						case uniq: t.converter = STREAM->STREAM.filter(S->!isEmpty(S)).map(converter).collect(Collectors.toSet()).stream().collect(Collectors.toList());break;
						case none: t.converter = STREAM->STREAM.collect(Collectors.toList());break;
						default: throw new IllegalArgumentException();
						}
					}
				}
			
			final String description = "Extracted from INFO/CSQ (vep)";
			if(VCFHeaderLineCount.INTEGER.equals(count)) {
				t.info = new VCFInfoHeaderLine(this.prefix+t.tag+suffix,1, type,description);
				}
			else
				{
				t.info = new VCFInfoHeaderLine(this.prefix+t.tag+suffix,count, type,description);
				}
			
			if(this.tags.stream().anyMatch(T->T.info.getID().equals(t.info.getID()))) {
				throw new IllegalArgumentException("duplicate tag "+t.info.getID());
				}
			this.tags.add(t);
			}
		}
	
	public void setPrefix(final String prefix) {
		this.prefix = prefix;
		}
	
	
	private boolean isEmpty(final String s) {
		return StringUtils.isBlank(s)|| s.equals(".");
	}
	
	@Override
	public void fillHeader(VCFHeader header) {
		this.parser =new VepPredictionParserFactory(header).get();
		
		for(Tag t: this.tags) {
			if(!this.parser.getCategories().contains(t.tag)) {
				LOG.warn("no "+t.tag+" in "+String.join("|",this.parser.getCategories()));
				continue;
				}
			}
				
		for(Tag t: this.tags) {
			if(header.getInfoHeaderLine(t.info.getID())!=null) {
				throw new JvarkitException.DuplicateVcfHeaderInfo(header, t.info.getID());
				}
			header.addMetaDataLine(t.info);
			}
		}
	
	
	
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
		if(!ctx.hasAttribute(this.parser.getTag())) return Collections.singletonList(ctx);
		final List<VepPrediction> predictions = this.parser.getPredictions(ctx);
		final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
		for(Tag t:this.tags) {
			vcb.rmAttribute(t.info.getID());//reset if exists
			final Object o = t.converter.apply(predictions.stream().map(PRED->PRED.get(t.tag)));
			if(o==null) continue;
			if(t.info.getType().equals(VCFHeaderLineType.Flag)) {
				if(o.equals(Boolean.FALSE) || o.toString().equals("false")) continue;
				}
			vcb.attribute(t.info.getID(), o);
			}
		return Collections.singletonList(vcb.make());
		}
	
	@Override
	public void close() {
		
		}
}
