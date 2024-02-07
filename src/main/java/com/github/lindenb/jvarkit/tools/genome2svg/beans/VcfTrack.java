package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

public class VcfTrack extends Track {
	private BiPredicate<VCFHeader,VariantContext> filter = (H,V)->true;
	private  BiFunction<VCFHeader, VariantContext, String> variantToString = (H,V)->V.getContig()+":"+V.getStart()+" "+V.getType();
	private  BiFunction<VCFHeader, VariantContext, String> variantToStyle =  (H,V)->"fill:red;stroke:blue;";
	private  BiFunction<VCFHeader, VariantContext, String> variantToURL =  (H,V)->"";
	public VcfTrack() {
		setFeatureHeight(10);
		}
	
	
	@Override
	public void paint(SVGContext ctx) {
		
		final Element g = ctx.element("g");
		ctx.tracksNode.appendChild(g);
		insertTitle(g,ctx);
		final double featureHeight = getFeatureHeight();
		
		try(VCFReader reader = VCFReaderFactory.makeDefault().open(getPath().asPath(), true)) {
			final VCFHeader header= reader.getHeader();
			try(CloseableIterator<VariantContext> iter = reader.query(ctx.loc)) {
				while(iter.hasNext()) {
					final VariantContext vc = iter.next();
					if(!this.filter.test(header, vc)) continue;
					final Element rect = ctx.rect(vc,ctx.y,featureHeight);
					
					final String url = variantToURL.apply(header, vc);
					g.appendChild(ctx.anchor(rect, url));
						
					String title = this.variantToString.apply(header,vc);
					if(!StringUtils.isBlank(title)) rect.appendChild(ctx.title(title));
					
					String style = variantToStyle.apply(header, vc);
					if(!StringUtils.isBlank(style)) rect.setAttribute("style", style);

					}
				}
			ctx.y += featureHeight+1;
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	}
