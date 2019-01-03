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

/*
 
  ant jeter && cat jeter2.bed | java -jar dist/jeter.jar  -u src/main/resources/img/World_V2.0.svg -o ~/jeter.jpg -R /commun/data/pubdb/broadinstitute.org/bundle/1.5/b37/human_g1k_v37.fasta && eog ~/jeter.jpg
  
   gunzip -c ~/jeter.gz |  gawk -f countries.awk | sed 's/^chr//' | grep -v '_' > jeter2.bed
  
  
 */

package com.github.lindenb.jvarkit.tools.circular;


import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.Source;
import javax.xml.transform.stream.StreamSource;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.circular.CircularContext;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**
BEGIN_DOC

## Example

```
head map.bed

1	141535	143008	taiwan
1	564463	570304	china
1	564463	570304	china
1	564463	564813	canada
1	564510	569779	estonia
1	564510	569779	estonia
1	564510	569170	germany
1	564633	569139	italy
1	564660	565200	iran
1	564733	569130	brazil
1	564776	569431	mexico
```


```bash
$  cat map.bed |\
     java -jar dist/worldmapgenome.jar \
     -u World_V2.0.svg \
     -w 2000 -o ~/ouput.jpg \
     -R human_g1k_v37.fasta

```
![worldmapgenome](https://pbs.twimg.com/media/BfGE0X4CMAAfRAR.jpg)

## See also

* https://twitter.com/yokofakun/status/428269474869817344
* https://twitter.com/GenomeBrowser/status/426398651103997953



END_DOC

 */


@Program(name="worldmapgenome",
	description="Genome/Map of the world. Input is a BED file: chrom/start/end/country.",
	keywords="gis"
	)
public class WorldMapGenome extends Launcher
	{
	private static final Logger LOG = Logger.build(WorldMapGenome.class).make();


	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	File faidx=null;
	@Parameter(names="-w",description=" (int: size) square size")
	private int squareSize=1000;
	@Parameter(names="-T",description=" just list countries and exit.",help=true)
	private  boolean listCountries=false;
	@Parameter(names="-o",description=" (file.jpg) ouput file. Required",required=true)
	private File fileout=null;
		
	private CircularContext context=null;
	/* ouput size */
	private Dimension viewRect=new Dimension(1000, 1000);
    private Map<String,Shape> country2shape=new HashMap<String,Shape>();
    /* SVG map location */
    @Parameter(names="-u",description=" (uri) URL world SVG map.")
    private String svgMapUri="http://upload.wikimedia.org/wikipedia/commons/7/76/World_V2.0.svg";
    /* number of time we saw each country */
    private Counter<String> seen=new Counter<String>();
	
    private WorldMapGenome()
		{
		}
  
		
	private void loadWorld() throws IOException,XMLStreamException
		{
		Source source;
		LOG.warning("openingg "+svgMapUri);
		if(IOUtils.isRemoteURI(svgMapUri))
			{
			source=new StreamSource(svgMapUri);
			}
		else
			{
			source=new StreamSource(new File(svgMapUri));
			}
			
		
		XMLInputFactory xif=XMLInputFactory.newFactory();
		XMLEventReader xef=xif.createXMLEventReader(source);
		while(xef.hasNext())
			{
			XMLEvent evt=xef.nextEvent();
			if(!evt.isStartElement())
				{
				continue;
				}
			StartElement E=evt.asStartElement();
			String localName=E.getName().getLocalPart();
			if(!localName.equals("path")) continue;
			Attribute att=E.getAttributeByName(new QName("id"));
			if(att==null) continue;
			String country=att.getValue().toLowerCase().replaceAll("[ ]+", "");
			att=E.getAttributeByName(new QName("d"));
			if(att==null) continue;
			GeneralPath path=null;
			char op='\0';
			Scanner scanner=new Scanner(att.getValue().replaceAll("[ \t\n\r,]+", " "));
			path=new GeneralPath();

			while(scanner.hasNext())
				{
				if(op=='\0')
					{
					op=scanner.next().charAt(0);
					}
				switch(op)
					{
					case 'M':
							
							path.moveTo(
								scanner.nextDouble(),
								scanner.nextDouble()
								);
							break;
					case 'C':
						path.curveTo(
								scanner.nextDouble(),
								scanner.nextDouble(),
								scanner.nextDouble(),
								scanner.nextDouble(),
								scanner.nextDouble(),
								scanner.nextDouble()
								);
						break;
					case 'Z':
						{
						path.closePath();
						
						break;
						}
					default:throw new IOException("bad operator "+op);
					}
				if(scanner.hasNext("[MCZ]"))
					{
					op=scanner.next().charAt(0);
					}
				}
			scanner.close();
			this.country2shape.put(country, scaleWorld(path));
			}
		xef.close();
		}
	
	private Shape scaleWorld(Shape shape)
		{	
		Dimension svgDim=new Dimension(800,400);

		double diag=Math.sqrt(Math.pow(svgDim.width/2, 2)+Math.pow(svgDim.height/2, 2));
		double ratio=getRadiusInt()/diag;
		
		//move to center of figure
		AffineTransform tr=AffineTransform.getTranslateInstance(-svgDim.width/2, -svgDim.height/2);
		shape=tr.createTransformedShape(shape);
		tr=AffineTransform.getScaleInstance(ratio,ratio);
		shape=tr.createTransformedShape(shape);
		tr=AffineTransform.getTranslateInstance(viewRect.width/2,viewRect.height/2);
		shape=tr.createTransformedShape(shape);
		
		
		
		return  shape;
		}

	@Override
	public int doWork(List<String> args) {
		try
			{
			viewRect.width=squareSize;
			viewRect.height=squareSize;
			loadWorld();
			if(listCountries)
				{
				for(String country:this.country2shape.keySet())
					{
					System.out.println(country);
					}
				return 0;
				}
			if(faidx==null)
				{
				LOG.error("undefined reference file");
				return -1;
				}
			if(fileout==null)
				{
				LOG.error("undefined output file");
				return -1;
				}
			
			LOG.info("loading "+faidx);
			this.context=new CircularContext(SAMSequenceDictionaryExtractor.extractDictionary(faidx));
			this.context.setCenter(viewRect.width/2,viewRect.height/2);
			BufferedImage offscreen=new  BufferedImage(
					viewRect.width,
					viewRect.height, 
					BufferedImage.TYPE_INT_ARGB);
			Graphics2D g=(Graphics2D)offscreen.getGraphics();
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				scan(g,stdin());
				}
			else
				{
				for(String filename :args)
					{
					LOG.info("Reading from "+filename);
					InputStream in=IOUtils.openURIForReading(filename);
					scan(g,in);
					CloserUtil.close(in);
					}
				}
			g.dispose();
			
			
			BufferedImage background=new  BufferedImage(
					viewRect.width,
					viewRect.height, 
					BufferedImage.TYPE_INT_RGB);
			g=(Graphics2D)background.getGraphics();
			//g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			//paint background
			paintGenome(g);
			paintWorld(g);
			g.drawImage(offscreen, 0, 0, null);
			g.dispose();
			ImageIO.write(background, "jpg", fileout);
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	
	private int getRadiusInt()
		{
		return getRadiusExt()-20;
		}
	private int getRadiusExt()
		{
		return Math.min(viewRect.height, viewRect.width)/2-20;
		}
	private void paintWorld(Graphics2D g)
		{
		for(String country:this.country2shape.keySet())
			{
			Shape shape=country2shape.get(country);
			float color=0;
			if(seen.count(country)!=0)
				{
				color=1f-(float)(Math.log(seen.count(country))/Math.log((double)seen.getTotal()));
				}
			g.setColor(seen.count(country)==0?Color.WHITE:new Color(color,0f,0f));
			g.fill(shape);
			}
		for(String country:this.country2shape.keySet())
			{
			Shape shape=country2shape.get(country);
			g.setColor(Color.BLACK);
			g.draw(shape);
				
			}

		}
	private void paintGenome(Graphics2D g)
		{
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, viewRect.width, viewRect.height);
		g.setColor(new Color(0.9f,0.9f,1f));
		g.fill(new Ellipse2D.Double(context.getCenterX()-getRadiusInt(),context.getCenterY()-getRadiusInt(),getRadiusInt()*2, getRadiusInt()*2));

		for(int i=0;i< this.context.getDictionary().getSequences().size();++i)
			{
			SAMSequenceRecord rec=this.context.getDictionary().getSequence(i);
			
			Arc2D outer=context.getArc(rec, 0, rec.getSequenceLength(), 
					getRadiusExt(),
					Arc2D.PIE);
			if(outer.getAngleExtent()==0) continue;
			Area area=new Area(outer);
			Ellipse2D.Double ed=new Ellipse2D.Double(
					context.getCenter().getX()-getRadiusInt(),
					context.getCenter().getY()-getRadiusInt(),
					getRadiusInt()*2,
					getRadiusInt()*2
					);
			area.subtract(new Area(ed));
			
						
			
			g.setColor(i%2==0?Color.LIGHT_GRAY:Color.WHITE);
			g.fill(area);
			g.setColor(Color.BLACK);
			g.draw(area);
			
			if((rec.getSequenceLength()/(double)this.context.getDictionary().getReferenceLength())<0.01) continue;
			String title=rec.getSequenceName();
			double midangle=context.convertPositionToRadian(
					rec,
					rec.getSequenceLength()/2
					);
		
			AffineTransform old=g.getTransform();
			AffineTransform tr=new AffineTransform(old);
			
			
			g.translate(context.getCenterX(),context.getCenterY());
			g.rotate(midangle);
			g.translate(getRadiusExt(),0);
			
			g.drawString(
					title,
					0,
					0
					);
			g.setTransform(tr);
			g.setTransform(old);	
			}
		}
	
	private void scan(Graphics2D g,InputStream input) throws IOException
		{
		Set<String> unknownC=new HashSet<String>();
		Pattern tab=Pattern.compile("[\t]");
		LineIterator in=new LineIteratorImpl(new SynchronousLineReader(input));
		while(in.hasNext())
			{	
			String line=in.next();
			String tokens[]=tab.split(line,5);
			if(tokens.length<4)
				{
				LOG.warning("Ignoring "+line);
				continue;
				}
			SAMSequenceRecord rec=this.context.getDictionary().getSequence(tokens[0]);
			if(rec==null)
				{
				LOG.warning("unknown chromosome "+tokens[0]);
				continue;
				}
			String country=tokens[3].toLowerCase().replaceAll("[ ]", "");
			Shape shape=this.country2shape.get(country);
			if(shape==null)
				{
				if(!unknownC.contains(country))
					{
					unknownC.add(country);
					LOG.warning("unknown country "+country);
					}
				continue;
				}
			seen.incr(country);
			int midpos=(Integer.parseInt(tokens[1])+Integer.parseInt(tokens[2]))/2;
			//country center
			Point2D.Double pt1 =new Point2D.Double(shape.getBounds2D().getCenterX(),shape.getBounds2D().getCenterY());
			//circle point
			Point2D pt3= this.context.convertPositionToPoint(tokens[0],midpos,getRadiusInt());
			double angle= this.context.convertPositionToRadian(rec, midpos);
			double angle2=angle-=Math.PI/10.0;
			
			double distance13= context.getCenter().distance(new Point2D.Double(
					(pt1.getX()+pt3.getX())/2.0,
					(pt1.getY()+pt3.getY())/2.0
					));
			//mid point
			Point2D pt2 =new Point2D.Double(
					context.getCenterX()+distance13*Math.cos(angle2),
					context.getCenterX()+distance13*Math.sin(angle2)
					);
			
			Composite old=g.getComposite();
			Stroke olds=g.getStroke();
			g.setStroke(new BasicStroke(0.8f));
			g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.02f));
			g.setColor(Color.DARK_GRAY);
			GeneralPath p=new GeneralPath();
			p.moveTo(pt1.getX(), pt1.getY());
			p.quadTo(pt2.getX(), pt2.getY(),pt3.getX(), pt3.getY());
			p.closePath();
			g.draw(p);
			g.setComposite(old);
			g.setStroke(olds);
			}
		CloserUtil.close(in);
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new WorldMapGenome().instanceMainWithExit(args);
		}

	}
