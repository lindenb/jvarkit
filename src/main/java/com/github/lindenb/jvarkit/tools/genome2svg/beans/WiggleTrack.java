package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.wig.BigWigReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

public class WiggleTrack extends Track {
	private List<String> wiggleFiles = new ArrayList<>();
	
	protected List<String> getWigPaths() {
		return this.wiggleFiles;
		}
	
	@Override
	public void paint(SVGContext ctx) {
		Float minV=null;
		Float maxV=null;
		final List<String> wigglePaths = getWigPaths();
		/* first, get min/max for each wiggle file */
		for(String path: wigglePaths) {
			try(BigWigReader bbr = new BigWigReader(path)) {
				try(CloseableIterator<BigWigReader.WigItem> iter=bbr.query(ctx.loc)) {
					while(iter.hasNext()) {
						final BigWigReader.WigItem item=iter.next();
						if(!item.overlaps(ctx.loc)) continue;//paranoid
						final float v=item.getValue();
						if(minV==null) {
							minV = v;
							maxV = v;
							}
						else
							{
							minV = Math.min(minV, v);
							maxV = Math.max(maxV, v);
							}
						}
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		final Element g0 = ctx.element("g");
		/* plot each file */
		for(String path: wigglePaths) {
			try(BigWigReader bbr = new BigWigReader(path)) {
				final Element g = ctx.element("g");
				g0.appendChild(g);
				int prev_pos=-1;
				try(CloseableIterator<BigWigReader.WigItem> iter=bbr.query(ctx.loc)) {
					while(iter.hasNext()) {
						final BigWigReader.WigItem item=iter.next();
						if(!item.overlaps(ctx.loc)) continue;//paranoid
						final float v=item.getValue();
						
						}
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
	
		}
		
	

	
	
}
