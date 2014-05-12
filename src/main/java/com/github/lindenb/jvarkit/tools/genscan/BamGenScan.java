package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMFileReader.ValidationStringency;
import htsjdk.samtools.util.SequenceUtil;

/**
 * GenScan
 *
 */
public class BamGenScan extends AbstractGeneScan
	{
	private List<Input> inputs=new ArrayList<>();
	
	private class Input
		{
		Sample sample;
		String filenames;
		}
	
	
	protected void drawPoints(Graphics2D g)
		{
		for(int input_id=0;input_id<this.inputs.size();++input_id)
			{		
			Input input=this.inputs.get(input_id);
			List<SAMFileReader> readers=new ArrayList<SAMFileReader>();
			for(String filename:input.filenames.split("[\\:]+"))
				{
				if(filename.isEmpty()) continue;
				info("Reading header for "+filename);
				SAMFileReader sfr=new SAMFileReader(new File(filename));
				sfr.setValidationStringency(ValidationStringency.SILENT);
				readers.add(sfr);
				
				}
			
			for(ChromInfo ci:this.chromInfos)
				{
				if(!ci.visible) continue;
				info("["+input_id+"]alloc "+ci.dictSequenceLength+" for "+ci.getSequenceName());
				int count[]=new int[ci.dictSequenceLength];
				Arrays.fill(count, 0);
				
				
				Shape clip=g.getClip();
				g.setClip(new Rectangle2D.Double(
						ci.x,input.sample.y,
						ci.width,input.sample.height
						));
				
				for(SAMFileReader sfr:readers) 
					{
					SAMRecordIterator sli=sfr.query(ci.getSequenceName(), 0, 0,true);
					
					//sli.setEmitUncoveredLoci(true);
					while(sli.hasNext())
						{
						SAMRecord rec=sli.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(rec.getDuplicateReadFlag()) continue;
						
						for(int pos1=rec.getAlignmentStart();pos1<= rec.getAlignmentEnd();++pos1)
							{
							int pos0=pos1-1;
							if(pos0<0 || pos0>=count.length) continue;
							count[pos0]++;
							}
						}
					sli.close();
					}
				int win_size=(int)Math.ceil(ci.dictSequenceLength.doubleValue()/ci.width);
				if(win_size<0) win_size=1;
				Point2D.Double prev_points[]=new Point2D.Double[4];
				Point2D.Double curr_points[]=new Point2D.Double[4];
				for(int n=0;n+win_size <= count.length;n+=win_size)
					{
					MinMaxDouble minmax_val=new MinMaxDouble();
					double count_2=0;
					double total_2=0;
					double array_median[]=new double[win_size];
					for(int i=0;i< win_size;++i)
						{
						minmax_val.visit(count[i+n]);
						count_2 += count[i+n];
						total_2 += 1;
						array_median[i]=count[i+n];
						}
					if(total_2==0 || !minmax_val.isValid()) continue;
					Arrays.sort(array_median);
					
					for(int k=0;k< curr_points.length;++k)
						{
						curr_points[k]=new Point2D.Double();
						curr_points[k].x= ci.x+(((double)n/ci.minmaxBase.getMax())*ci.width);
						}
					
					curr_points[0].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction((count_2)/total_2);
					curr_points[1].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction(minmax_val.getMin());
					curr_points[2].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction(minmax_val.getMax());
					curr_points[3].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction(array_median[array_median.length/2]);
					
					
					for(int k=0;k< curr_points.length;++k)
						{
						if(prev_points[k]!=null)
							{
							switch(k)
								{
								case 0: g.setColor(Color.GREEN);break;
								case 1: g.setColor(Color.YELLOW);break;
								case 2: g.setColor(Color.RED);break;
								case 3: g.setColor(Color.BLUE);break;
								}
							g.draw(new Line2D.Double(prev_points[k],curr_points[k]));
							}
						prev_points[k]=curr_points[k];
						}
					}	
				g.setClip(clip);
				} //end for ci
			for(SAMFileReader sfr:readers) sfr.close();
			}//end for input
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BamGenScan";
		}

	
	@Override
	protected List<ChromInfo> getChromInfos() {
		return this.chromInfos;
		}
	
	@Override
	public String getProgramName()
		{
		return "BamGenScan";
		}
	
	@Override
	public String getProgramDescription() {
		return "Paint a Genome Scan picture from a BAM.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (file.jpg) picture filename out. if undefined, show a Window");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);

		File filout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
			{
				case 'o': filout=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		try
			{
			
			if(opt.getOptInd()==args.length)
				{
				error("illegal number of arguments");
				return -1;
				}
			SAMSequenceDictionary dict=null;
			for(int i=opt.getOptInd();i<args.length;++i)
				{
				Input input=new Input();
				input.filenames=args[i];
				input.sample=new Sample();
				input.sample.name=input.filenames;
				input.sample.sample_id=super.samples.size();
				input.sample.minmax.setMin(0);
				input.sample.minmax.setMax(50);
				super.samples.add(input.sample);
				this.inputs.add(input);
				
				for(String filename:input.filenames.split("[\\:]+"))
					{
					if(filename.isEmpty()) continue;
					info("["+inputs.size()+"]Reading header for "+filename);
					SAMFileReader sfr=new SAMFileReader(new File(filename));
					sfr.setValidationStringency(ValidationStringency.SILENT);
					SAMFileHeader h=sfr.getFileHeader();
					sfr.close();
					SAMSequenceDictionary d=h.getSequenceDictionary();
					if(dict==null)
						{
						dict=d;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(d, dict))
						{
						error("Sequence dictionaries are not the same.");
						return -1;
						}
					}
				
				}
			info("number of input:"+inputs.size());
			if(dict==null )
				{
				error("No dictionary found");
				return -1;
				}
		
			for(SAMSequenceRecord rec: dict.getSequences())
				{
				ChromInfo ci=new ChromInfo();
				ci.dictSequenceLength=rec.getSequenceLength();
				ci.sequenceName=rec.getSequenceName();
				ci.minmaxBase.setMin(0);
				ci.minmaxBase.setMax(ci.dictSequenceLength.doubleValue());
				
				ci.tid=chromInfos.size();
				if(ci.tid!=rec.getSequenceIndex()) throw new IllegalStateException("boum");
				chromInfos.add(ci);
				}
				
			
			
			BufferedImage img=makeImage();

			
			if(filout==null)
				{
				showGui(img);
				}
			else
				{
				ImageIO.write(img, "JPG", filout);
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			}
		
		}
	public static void main(String[] args)
		{
		new BamGenScan().instanceMainWithExit(args);
		}
	
	}
