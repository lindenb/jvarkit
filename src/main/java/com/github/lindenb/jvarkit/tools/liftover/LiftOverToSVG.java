/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.liftover;

import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## Screenshot

https://twitter.com/yokofakun/status/525014376387215360

![svg ](https://pbs.twimg.com/media/B0k5ytFCIAAiBPJ.png)

## Examples

### Examples 1

hg16 to hg38 http://cardioserve.nantes.inserm.fr/~lindenb/liftover2svg/hg16ToHg38.svg

BIG FILE : hg38 to panTro4 http://cardioserve.nantes.inserm.fr/~lindenb/liftover2svg/hg38ToPanTro4.svg



### Example 2

create a **SVG** file for **hg16** to **hg38**.

```makefile
.PHONY:all 
all: test.svg


test.svg: hg16ToHg17.over.chain hg17ToHg18.over.chain hg18ToHg19.over.chain hg19ToHg38.over.chain dist/liftover2svg.jar
	java -jar dist/liftover2svg.jar -b hg16 -b hg17 -b hg18 -b hg19 -b hg38 \
			 hg16ToHg17.over.chain hg17ToHg18.over.chain hg18ToHg19.over.chain hg19ToHg38.over.chain > $@

hg16ToHg17.over.chain : 
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg16/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@

hg17ToHg18.over.chain : 
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg17/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@

hg18ToHg19.over.chain : 
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg18/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@
	
hg38ToPanTro4.over.chain :
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@

hg19ToHg38.over.chain hg19ToPanTro4.over.chain :
	curl -s "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/liftOver/$@.gz" | gunzip -c |\
	grep chain > $@


dist/liftover2svg.jar: ./src/main/java/com/github/lindenb/jvarkit/tools/liftover/LiftOverToSVG.java
	$(MAKE)  liftover2svg

```

END_DOC
 *
 */
@Program(name="liftover2svg",
	description="Convert LiftOver chain files to animated SVG",
	keywords={"svg","liftover","ucsc","xml"}	
	)
public class LiftOverToSVG extends Launcher
	{
	private final Logger LOG=Logger.build(LiftOverToSVG.class).make();

	private class ChromName
		{
		String name;
		ChromName(String name)
			{
			this.name=name;
			}
		@Override
		public int hashCode()
			{
			return this.name.hashCode();
			}
		@Override
		public boolean equals(Object obj)
			{
			if(obj==null) return false;
			if(obj==this) return true;
			return ChromName.class.cast(obj).name.equals(name);
			}
		
		List<Rectangle2D.Double> getBounds()
			{
			List<Rectangle2D.Double> rects=new ArrayList<Rectangle2D.Double>(LiftOverToSVG.this.builds.size());
			for(Build build: LiftOverToSVG.this.builds)
				{
				Chromosome c=build.name2chrom.get(this.name);
				rects.add(c==null?null:c.getBounds());
				}
			return rects;
			}
		
		public int getMaxLength()
			{
			int L=0;
			for(Build build: LiftOverToSVG.this.builds)
				{
				Chromosome c=build.name2chrom.get(this.name);
				if(c==null) continue;
				L=Math.max(L, c.length);
				}
			return L;
			}
		
		@Override
		public String toString()
			{
			return name;
			}
		}	
	
	private class Chromosome implements Comparable<Chromosome>
		{
		@SuppressWarnings("unused")
		Build build;
		String name;
		int length;
		int index=-1;
		
		Rectangle2D.Double getBounds()
			{
			return new Rectangle2D.Double(
					0,
					this.index*LiftOverToSVG.this.featureHeight,
					baseToPos(this.length),
					LiftOverToSVG.this.chromosomeHeigh
					);
			}
		
		@Override
		public int compareTo(Chromosome o)
			{
			int i= o.length-this.length;
			if(i!=0) return i;
			return this.name.compareTo(o.name);
			}
		}
	
	private class Build
		{
		String name=null;
		int index=-1;
		Map<String,Chromosome> name2chrom=new HashMap<String,Chromosome>();
		List<Chromosome> ordered=null;
		
		@SuppressWarnings("unused")
		float getTimeBegin()
			{
			return this.index* this.getTimeDuration();
			}
		
		float getTimeDuration()
			{
			return LiftOverToSVG.this.secondsPerStep;
			}
		
		Chromosome getChromosomeByName(String name,int length)
			{
			Chromosome c=this.name2chrom.get(name);
			if(c==null)
				{
				c=new Chromosome();
				ChromName cN=new ChromName(name);
				LiftOverToSVG.this.chromNames.add(cN);
				c.name=name;
				c.length=length;
				c.build=this;
				this.name2chrom.put(name, c);
				
				LiftOverToSVG.this.maxCountChrom=Math.max(this.name2chrom.size(), LiftOverToSVG.this.maxCountChrom);
				LiftOverToSVG.this.maxLength=Math.max(length, LiftOverToSVG.this.maxLength);
				}
			return c;
			}
		public void calcChromosomes()
			{
			this.ordered=new ArrayList<Chromosome>(this.name2chrom.values());
			Collections.sort(this.ordered);
			for(int i=0;i< this.ordered.size();++i) this.ordered.get(i).index=i;
			}
		
		public String getName()
			{
			return this.name==null?"Build "+this.index:this.name;
			}
		
		@Override
		public String toString() {
			return getName();
			}
		}
	
	private class Segment
		{
		Chromosome chrom;
		int chromStart;
		int chromEnd;
		char strand;
		
		Rectangle2D.Double getBounds()
			{
			double x = baseToPos(chromStart);
			double width= baseToPos(chromEnd)-x;
			return new Rectangle2D.Double(
					x,
					chrom.index*LiftOverToSVG.this.featureHeight+5,
					width,
					LiftOverToSVG.this.chromosomeHeigh-10
					);
			}
		@Override
		public String toString()
			{
			return chrom.name+":"+chromStart+"-"+chromEnd;
			}
		}
	
	private class Chain
		{
		long score;
		Segment start;
		Segment end;
		}
	
	private Set<ChromName> chromNames=new HashSet<ChromName>();
	private List<Build> builds=new ArrayList<Build>();
	private List<List<Chain>> chains=new ArrayList<List<Chain>>();
	/** max chromosome length */
	private int maxLength=0;
	/** max number of chromosome per build */
	private int maxCountChrom=0;
	/** svg width */
	@Parameter(names="-w",description="width")
	private int width=1000;
	/** distance between two chroms */
	private int featureHeight=30;
	/** chromosome height */
	private int chromosomeHeigh=25;
	/** duration */
	@Parameter(names="-t",description=" seconds per step")
	private float secondsPerStep=20f;
	
	@Parameter(names="-b",description="add this build name")
	private List<String> buildNames=new ArrayList<>();

	
	private LiftOverToSVG()
		{
		
		}
	
	
	
    
    private void readChain(String uri) throws IOException
    	{
    	LOG.info("Reading chain file "+uri);
    	BufferedReader in=IOUtils.openURIForBufferedReading(uri);
    	readChain(in);
    	in.close();
    	}
    private void readChain(BufferedReader in) throws IOException
    	{
    	Build build1;
    	if(this.builds.isEmpty())
    		{
    		build1=new Build();
    		build1.index=0;
    		this.builds.add(build1);
    		}
    	else
    		{
    		build1=this.builds.get(this.builds.size()-1);
    		}
    	Build build2=new Build();
    	build2.index=this.builds.size();
    	this.builds.add(build2);
    	List<Chain> chain=new ArrayList<Chain>();
    	this.chains.add(chain);
    	
    	final Pattern tab=Pattern.compile("[ \t]+");
    	String line;
    	while((line=in.readLine())!=null)
    		{
    		if(line.isEmpty() || !line.startsWith("chain")) continue;
    		String tokens[]=tab.split(line);
    		if(!tokens[0].equals("chain")) continue;
    		
    		
    		long score=Long.parseLong(tokens[1]);
    		String chromName = tokens[2];
    		Chromosome chr1 = build1.getChromosomeByName(chromName,Integer.parseInt(tokens[3]));
    		Segment seg1=new Segment();
    		seg1.chrom=chr1;
    		seg1.strand= tokens[4].charAt(0);
    		seg1.chromStart=Integer.parseInt(tokens[5]);
    		seg1.chromEnd=Integer.parseInt(tokens[6]);
    		
    		chromName = tokens[7];
    		Chromosome chr2 = build2.getChromosomeByName(chromName,Integer.parseInt(tokens[8]));
    		
    		Segment seg2=new Segment();
    		seg2.chrom=chr2;
    		seg2.strand= tokens[9].charAt(0);
    		seg2.chromStart=Integer.parseInt(tokens[10]);
    		seg2.chromEnd=Integer.parseInt(tokens[11]);
    		
    		
    		
    		
    		Chain ch=new Chain();
    		ch.score= score;
    		ch.start=seg1;
    		ch.end=seg2;
    		
    		chain.add(ch);
    		}
    	LOG.info("chain.size: "+chain.size());
    	}
    private double baseToPos(int pos)
    	{
    	return (pos/(double)this.maxLength)*this.width;
    	}
    
    private void write(XMLStreamWriter w) throws XMLStreamException
    	{
    	final String repeatCount="1";
    	final String blues[]={"#EEEEEE","#DDDDDD"};
    	w.writeStartDocument("UTF-8", "1.0");
    	w.writeStartElement("svg");
    	w.writeAttribute("style","fill:white;stroke:black;");
    	w.writeAttribute("width",String.valueOf(this.width+100));
    	double svgHeight=100+this.featureHeight*this.maxCountChrom;
    	w.writeAttribute("height",String.valueOf(svgHeight));
    	w.writeDefaultNamespace(SVG.NS);
    	
    	
		w.writeStartElement(SVG.NS,"title");
		w.writeCharacters("liftover2svg");
		w.writeEndElement();
		
		w.writeStartElement(SVG.NS,"description");
		w.writeCharacters("Cmd:"+getProgramCommandLine()+"\n");
		w.writeCharacters("Version:"+getVersion()+"\n");
		w.writeEndElement();

    	
		w.writeStartElement("style");
		w.writeCharacters(
				"rect.chain {fill:#99AAAA;stroke:blue;opacity:0;stroke-width:0.5;}\n"+
				"text.chromName {stroke:black;}\n"+
				"text.buildname {fill: #999999; stroke:#000000;font-size: 48px;opacity:0;text-anchor:end;}\n"
				);
		w.writeEndElement();
				
    	
    	w.writeStartElement("g");
		w.writeAttribute("transform","translate(50,50)");
    	
		//build-names
		w.writeStartElement("g");
		w.writeComment("BEGIN build names");
		for(int step=0;step+1<builds.size();++step)
			{
			final float opacityTimeFraction=0.1f;
			final float invisibleDuration2 = (this.secondsPerStep*opacityTimeFraction)/2f;
			final float visibleDuration = this.secondsPerStep - invisibleDuration2*2f;
			final float stepTimeStart= step * this.secondsPerStep;
			final float stepVisibleStart= stepTimeStart + invisibleDuration2 ;
			
			Build b= this.builds.get(step);
			w.writeStartElement("text");
			w.writeAttribute("x", String.valueOf(this.width));
			w.writeAttribute("y", String.valueOf(svgHeight/4.0));
			w.writeAttribute("class","buildname");
			
			
			w.writeCharacters(b.getName()+" > " + this.builds.get(step+1).getName() );
			
			//show title
			w.writeEmptyElement("animate");
			w.writeAttribute("attributeType","CSS");
			w.writeAttribute("attributeName","opacity");
			w.writeAttribute("begin",String.valueOf(stepVisibleStart));
			w.writeAttribute("dur",String.valueOf(invisibleDuration2));
			w.writeAttribute("from","0");
			w.writeAttribute("to","1");
			w.writeAttribute("repeatCount",repeatCount);
			w.writeAttribute("fill","freeze");
				
			//hide title
			w.writeEmptyElement("animate");
			w.writeAttribute("attributeType","CSS");
			w.writeAttribute("attributeName","opacity");
			w.writeAttribute("begin",String.valueOf(stepVisibleStart+visibleDuration-invisibleDuration2));
			w.writeAttribute("dur",String.valueOf(invisibleDuration2));
			w.writeAttribute("from","1");
			w.writeAttribute("to","0");
			w.writeAttribute("repeatCount",repeatCount);
			w.writeAttribute("fill","freeze");
			
			
			w.writeEndElement();//text
			}
		w.writeComment("END build names");
		w.writeEndElement();
		//end buidl name
		
		
		w.writeStartElement("g");
    	for(ChromName c: this.chromNames)
			{
			List<Rectangle2D.Double> bounds=c.getBounds();
			int step=0;
			Rectangle2D.Double first_rect=null;
			while((first_rect=bounds.get(step))==null) step++;
			
			w.writeStartElement("g");
			w.writeAttribute("title",c.name);
			
			for(int shape=0;shape<2;++shape)
				{
				w.writeStartElement(shape==0?"rect":"text");
	    		w.writeAttribute("x",String.valueOf(shape==0?first_rect.x:10+baseToPos(c.getMaxLength())));
	    		
	    		w.writeAttribute("y", String.valueOf(first_rect.y+(shape==0?0:this.chromosomeHeigh/2f)));
	    		w.writeAttribute("opacity",(step==0?"1":"0"));

	    		if(shape==0)
	    			{
		    		w.writeAttribute("width",String.valueOf(first_rect.width));
		    		w.writeAttribute("height",String.valueOf(first_rect.height));
	    			}
	    		else
	    			{
		    		w.writeAttribute("class","chromName");
		    		w.writeCharacters(c.name);
	    			}
				
				for(step=0;step+1 <bounds.size();++step)
					{
					Rectangle2D.Double rect= bounds.get(step);
					Rectangle2D.Double next= bounds.get(step+1);
					//opacity
					w.writeEmptyElement("animate");
					w.writeAttribute("attributeType","CSS");
					w.writeAttribute("attributeName","opacity");
					w.writeAttribute("begin",String.valueOf(step*secondsPerStep));
					w.writeAttribute("dur",String.valueOf(secondsPerStep));
					w.writeAttribute("from",(rect==null?"0":"1"));
					w.writeAttribute("to",(next==null?"0":"1"));
					w.writeAttribute("repeatCount",repeatCount);
					w.writeAttribute("fill","freeze");

					
					//fill
					if(shape==0)
						{
						w.writeEmptyElement("animate");
						w.writeAttribute("attributeType","CSS");
						w.writeAttribute("attributeName","fill");
						w.writeAttribute("begin",String.valueOf(step*secondsPerStep));
						w.writeAttribute("dur",String.valueOf(secondsPerStep));
						w.writeAttribute("from",(step%2==0?blues[0]:blues[1]));
						w.writeAttribute("to",(step%2==0?blues[1]:blues[0]));
						w.writeAttribute("repeatCount",repeatCount);
						w.writeAttribute("fill","freeze");
						}
					
					if(rect!=null && next!=null)
						{
						w.writeEmptyElement("animate");
						w.writeAttribute("attributeType","XML");
						w.writeAttribute("attributeName","y");
						w.writeAttribute("begin",String.valueOf(step*secondsPerStep));
						w.writeAttribute("dur",String.valueOf(secondsPerStep));
						w.writeAttribute("from",String.valueOf(rect.y+(shape==0?0:this.chromosomeHeigh/2)));
						w.writeAttribute("to",String.valueOf(next.y+(shape==0?0:this.chromosomeHeigh/2)));
						w.writeAttribute("repeatCount",repeatCount);
						w.writeAttribute("fill","freeze");
	
						if(shape==0)
							{
							w.writeEmptyElement("animate");
							w.writeAttribute("attributeType","XML");
							w.writeAttribute("attributeName","width");
							w.writeAttribute("begin",String.valueOf(step*secondsPerStep));
							w.writeAttribute("dur",String.valueOf(secondsPerStep));
							w.writeAttribute("from",String.valueOf(rect.width));
							w.writeAttribute("to",String.valueOf(next.width));
							w.writeAttribute("repeatCount",repeatCount);
							w.writeAttribute("fill","freeze");
							}
						}
					
					}
				w.writeEndElement();//rect
				}
			w.writeEndElement();//g			
    		}
    	w.writeEndElement();//g
    	
		w.writeStartElement("g");
		for(int step=0;step < this.chains.size();++step)
			{
			final float opacityTimeFraction=0.2f;
			final float invisibleDuration2 = (this.secondsPerStep*opacityTimeFraction)/2f;
			final float visibleDuration = this.secondsPerStep - invisibleDuration2*2f;
			final float stepTimeStart= step * this.secondsPerStep;
			final float stepVisibleStart= stepTimeStart + invisibleDuration2 ;
			final String chainOpacity="0.7";
			
			w.writeComment("BEGIN CHAIN["+step+"]");

			List<Chain> chainList = this.chains.get(step);
			//print longest first (higher score?)
			Collections.sort(chainList,new Comparator<Chain>()
				{
				@Override
				public int compare(Chain o1, Chain o2)
					{
					if (o1.score>o2.score) return -1;
					if (o1.score<o2.score) return 1;
					return 0;
					}
				});
			/* print all chains */
			for(Chain chain:chainList)
				{
				/* starting rectangle */
				Rectangle2D.Double rect=chain.start.getBounds();
				/* end rectangle */
				Rectangle2D.Double next=chain.end.getBounds();
				
				/* write initial rectangle */
				w.writeStartElement("rect");
				w.writeAttribute("x",String.valueOf(rect.x));
				w.writeAttribute("y",String.valueOf(rect.y));
				w.writeAttribute("width",String.valueOf(rect.width));
				w.writeAttribute("height",String.valueOf(rect.height));
				w.writeAttribute("class","chain");
				w.writeAttribute("title","score "+chain.score+" "+chain.start+" to  "+chain.end);

				/* at begin make it visible */
				w.writeEmptyElement("animate");
				w.writeAttribute("attributeType","CSS");
				w.writeAttribute("attributeName","opacity");
				w.writeAttribute("begin",String.valueOf(stepTimeStart));
				w.writeAttribute("dur",String.valueOf(invisibleDuration2));
				w.writeAttribute("from","0");
				w.writeAttribute("to",chainOpacity);
				w.writeAttribute("repeatCount",repeatCount);
				w.writeAttribute("fill","freeze");
				
				
				/* move rectangle.x  to destination.x */
				w.writeEmptyElement("animate");
				w.writeAttribute("attributeType","XML");
				w.writeAttribute("attributeName","x");
				w.writeAttribute("begin",String.valueOf(stepVisibleStart));
				w.writeAttribute("dur",String.valueOf(visibleDuration));
				w.writeAttribute("from",String.valueOf(rect.x));
				w.writeAttribute("to",String.valueOf(next.x));
				w.writeAttribute("repeatCount",repeatCount);
				w.writeAttribute("fill","freeze");
				
				/* move rectangle.y  to destination.y */
				w.writeEmptyElement("animate");
				w.writeAttribute("attributeType","XML");
				w.writeAttribute("attributeName","y");
				w.writeAttribute("begin",String.valueOf(stepVisibleStart));
				w.writeAttribute("dur",String.valueOf(visibleDuration));
				w.writeAttribute("from",String.valueOf(rect.y));
				w.writeAttribute("to",String.valueOf(next.y));
				w.writeAttribute("repeatCount",repeatCount);
				w.writeAttribute("fill","freeze");
				
				/* change rectangle.width  to destination.width */
				w.writeEmptyElement("animate");
				w.writeAttribute("attributeType","XML");
				w.writeAttribute("attributeName","width");
				w.writeAttribute("begin",String.valueOf(stepVisibleStart));
				w.writeAttribute("dur",String.valueOf(visibleDuration));
				w.writeAttribute("from",String.valueOf(rect.width));
				w.writeAttribute("to",String.valueOf(next.width));
				w.writeAttribute("repeatCount",repeatCount);
				w.writeAttribute("fill","freeze");

				/* if !=strand rotate around center.x,center.y */
				if(chain.start.strand!=chain.end.strand)
					{
					w.writeEmptyElement("animateTransform");
					w.writeAttribute("attributeType","XML");
					w.writeAttribute("attributeName","transform");
					w.writeAttribute("type","rotate");
					w.writeAttribute("begin",String.valueOf(stepVisibleStart));
					w.writeAttribute("dur",String.valueOf(visibleDuration));
					w.writeAttribute("repeatCount",repeatCount);

					w.writeAttribute("from","0 "+rect.getCenterX()+" "+rect.getCenterY());
					w.writeAttribute("to","180 "+next.getCenterX()+" "+next.getCenterY());
					w.writeAttribute("fill","freeze");
					}
				
				/* at begin make it invisible  */
				w.writeEmptyElement("animate");
				w.writeAttribute("attributeType","CSS");
				w.writeAttribute("attributeName","opacity");
				w.writeAttribute("begin",String.valueOf(stepTimeStart+this.secondsPerStep-invisibleDuration2));
				w.writeAttribute("dur",String.valueOf(invisibleDuration2));
				w.writeAttribute("from",chainOpacity);
				w.writeAttribute("to","0");
				w.writeAttribute("repeatCount",repeatCount);
				w.writeAttribute("fill","freeze");


				
				w.writeEndElement();//rect
				}
			w.writeComment("END CHAIN["+step+"]");
			}
		w.writeEndElement();//g

    	
    	
    	w.writeEndElement();//svg
    	w.writeEndDocument();
    	w.flush();
    	w.close();
    	}
    
	@Override
	public int doWork(List<String> args) {
		try {
			if(args.isEmpty())
				{
				readChain(new BufferedReader(new InputStreamReader(stdin())));
				}
			else
				{
				for(String fname:args)
					{
					readChain(fname);
					}
				}
		
			
			for(Build b:this.builds)b.calcChromosomes();
			
			for(int i=0;i< buildNames.size() && i< this.builds.size() ;++i)
				{
				this.builds.get(i).name=buildNames.get(i);
				}
			
			LOG.info("Max count chroms: "+this.maxCountChrom);
			LOG.info("Max length chroms: "+this.maxLength);
			LOG.info("Transitions: "+this.chains.size());
			
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLStreamWriter w=xof.createXMLStreamWriter(System.out, "UTF-8");
			write(w);
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new LiftOverToSVG().instanceMainWithExit(args);
		}

	}
