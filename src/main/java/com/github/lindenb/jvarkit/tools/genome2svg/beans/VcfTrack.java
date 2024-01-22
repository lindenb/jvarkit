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
	private String path = null;
	private BiPredicate<VCFHeader,VariantContext> filter = (H,V)->true;
	private  BiFunction<VCFHeader, VariantContext, String> variantToString = (H,V)->V.getContig()+":"+V.getStart()+" "+V.getType();
	private  BiFunction<VCFHeader, VariantContext, String> variantToStyle =  (H,V)->"fill:red;stroke:blue;";
	private  BiFunction<VCFHeader, VariantContext, String> variantToURL =  (H,V)->"";
	public void setVcf(String path) {
		this.path = path;
		}
	public String getVcf() {
		return path;
		}
	
	@Override
	public void paint(SVGContext ctx) {
		
		final Element g = ctx.element("g");
		ctx.tracksNode.appendChild(g);
		insertTitle(g,ctx);
		double featureHeight = 10;
		
		try(VCFReader reader = VCFReaderFactory.makeDefault().open(getVcf(), true)) {
			final VCFHeader header= reader.getHeader();
			try(CloseableIterator<VariantContext> iter = reader.query(ctx.loc)) {
				while(iter.hasNext()) {
					final VariantContext vc = iter.next();
					if(!this.filter.test(header, vc)) continue;
					Element rect = ctx.element("rect");
					
					String url = variantToURL.apply(header, vc);
					if(StringUtils.isBlank(url)) {
						g.appendChild(rect);
						}
					else
						{
						Element a = ctx.element("a");
						a.setAttribute("href", url);
						g.appendChild(a);
						a.appendChild(rect);
						}
						
					
					rect.setAttribute("x",String.valueOf(ctx.pos2pixel(ctx.trimpos(vc.getStart()))));
					rect.setAttribute("y",String.valueOf(ctx.y));
					rect.setAttribute("width",String.valueOf(
							ctx.pos2pixel(ctx.trimpos(vc.getEnd()+1))-ctx.pos2pixel(ctx.trimpos(vc.getStart()))));
					rect.setAttribute("height",String.valueOf(featureHeight));
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
