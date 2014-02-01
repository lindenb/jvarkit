package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;


import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;

public abstract class AbstractGeneScan extends AbstractCommandLineProgram
	{
	protected Dimension screenSize=new Dimension(1000,300);
    //@Option(shortName="DC",doc="Do not display unseen chromosome",optional=false)	
    //@Option(shortName="DC",doc="Do not display unseen chromosome",optional=false)	
	protected boolean DICARD_UNSEEN_CHROM=false;
    //@Option(shortName="DB",doc="ignore chromosomes left/right side if there's no data there.",optional=false)	
	protected boolean DISCARD_BOUNDS=false;
    private Hershey hershey=new Hershey();
    protected Insets insets=new Insets(20,150, 30, 30);
    protected double max_value=-Double.MAX_VALUE;
    protected double min_value=Double.MAX_VALUE;

    
	protected List<Sample> samples=new ArrayList<Sample>();
	
	protected class Sample
		{
		int sample_id=-1;
		String name="";
		double min_value=Double.MAX_VALUE;
		double max_value=Double.MIN_VALUE;
		double y;
		double weight;
		}

    
    
	protected class ChromInfo 
		{
		int tid=-1;
		String sequenceName=null;
		int sequenceLength=0;
		double x;
		double width;
		int max_pos=Integer.MIN_VALUE;
		int min_pos=Integer.MAX_VALUE;
		
		public double converPosToPixel(int pos)
			{
			if(DISCARD_BOUNDS) pos-=this.min_pos;
			return this.x+(pos)/((double)this.getSequenceLength())*this.width;
			}
		
		
		
		
		private Rectangle2D.Double getFrame()
			{
			return new Rectangle2D.Double(
					x,insets.top,
					width,
					screenSize.getHeight()-(insets.top+insets.bottom)
					);
			}
		
		void draw(Graphics2D g)
			{
			Rectangle2D.Double r=getFrame();
			g.setColor(Color.GRAY);
			g.draw(r);
			
			double optW=Math.min(this.getSequenceName().length()*10,this.width);
			Rectangle2D.Double titleRec=new Rectangle2D.Double(
					(x+width/2)-optW/2.0,
					insets.top-22,
					optW,
					20);
			
			hershey.paint(g, this.getSequenceName(), titleRec);
			}
		
		String getSequenceName()
			{
			return this.sequenceName;
			}
		
		int getSequenceLength()
			{
			if(DISCARD_BOUNDS)
				{
				return max_pos-min_pos;
				}
			else
				{
				return this.getSequenceLength();
				}
			}
		}
	
	
	
	protected AbstractGeneScan()
		{
		
		}
	protected abstract List<ChromInfo> getChromInfos();

	protected List<Sample> getSamples()
		{
		return this.samples;
		}
	
	protected abstract void drawPoints(Graphics2D g);
	
	
	private Rectangle2D.Double getFrame(ChromInfo ci,Sample sample)
		{
		return new Rectangle2D.Double(
				ci.x,
				sample.y,
				ci.width,
				sample.weight//screenSize.getHeight()-(insets.top+insets.bottom)
				);
		}

	
	protected BufferedImage makeImage()
		{
		BufferedImage img=new BufferedImage(
				this.screenSize.width,
				this.screenSize.height,
				BufferedImage.TYPE_INT_RGB
				);
		
		Graphics2D g=img.createGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0,this.screenSize.width, this.screenSize.height);
		
		
			for(int i=0;i< getChromInfos().size();++i)
				{
				g.setColor(i%2==0?Color.WHITE:Color.LIGHT_GRAY);
				g.fill(getChromInfos().get(i).getFrame());
				}
			
		
		//draw axis y
		g.setColor(Color.BLACK);
		for(int i=0;i<= 10;++i)
			{
			double v=min_value+ ((max_value-min_value)/10.0)*i;
			double y=(this.screenSize.height-insets.bottom)-
					((v-min_value)/(max_value-min_value))*(this.screenSize.height-(insets.bottom+insets.top))
					;
			g.setColor(Color.LIGHT_GRAY);
			g.draw(new Line2D.Double(insets.left-5, y, this.screenSize.width-insets.right, y));
			
			Rectangle2D.Double r=new Rectangle2D.Double(0,y-12,insets.left,24);
			g.setColor(Color.BLACK);
			this.hershey.paint(g, String.valueOf(v), r);
			}
		float OPACITY=0.1f;
		Composite oldComposite=g.getComposite();
		g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)OPACITY));
		drawPoints(g);
		g.setComposite(oldComposite);
		
		for(ChromInfo ci:getChromInfos())
			{
			ci.draw(g);
			}
		
		
		g.dispose();	
		return img;
		}
	}