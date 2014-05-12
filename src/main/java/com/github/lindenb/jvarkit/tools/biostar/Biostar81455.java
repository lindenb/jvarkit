package com.github.lindenb.jvarkit.tools.biostar;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Iterator;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Log;


public class Biostar81455 extends AbstractCommandLineProgram
	{
	
	
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Defining precisely the genomic context based on a position .";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="input process (or stdin): chromosome(tab)pos0.",
    		optional=true)
	public File IN=null;

 
	private Log LOG=Log.getInstance(Biostar81455.class);
	
	
    @Option(shortName="KG",doc="KnownGene data URI/File. should look like http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz . Beware chromosome names are formatted the same as your REFERENCE.",optional=false)
	public String kgUri="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";
	
    private static int distance(int pos0,KnownGene kg)
		{
		if(pos0<kg.getTxStart()) return kg.getTxStart()-pos0;
		if(pos0>kg.getTxEnd())  return kg.getTxEnd()-pos0;
		return 0;
		}
    
    private static int distance(int pos0,KnownGene.Exon exon)
    	{
    	if(pos0<exon.getStart()) return exon.getStart()-pos0;
    	if(pos0>exon.getEnd())  return exon.getEnd()-pos0;
    	return 0;
    	}
    
    
    
	@Override
	protected int doWork()
		{
		IntervalTreeMap<KnownGene> kgMap=new IntervalTreeMap<KnownGene>();
		try
			{
			LOG.info("readubf "+kgUri);
			String line;
			Pattern tab=Pattern.compile("[\t]");
			BufferedReader r=IOUtils.openURIForBufferedReading(this.kgUri);
			while((line=r.readLine())!=null)
				{
				String tokens[]=tab.split(line);
				
				KnownGene kg=new KnownGene();
				kg.setName(tokens[0]);
				kg.setChrom(tokens[1]);
				kg.setStrand(tokens[2].charAt(0));
				kg.setTxStart(Integer.parseInt(tokens[3]));
				kg.setTxEnd(Integer.parseInt(tokens[4]));
				kg.setCdsStart(Integer.parseInt(tokens[5]));
				kg.setCdsEnd(Integer.parseInt(tokens[6]));
				kg.setExonBounds(Integer.parseInt(tokens[7]), tokens[8], tokens[9]);
				if(kg.getExonCount()==0) continue;
				
				kgMap.put(new Interval(kg.getChromosome(), kg.getTxStart()+1, kg.getTxEnd()), kg);
				}
			r.close();
			
			
			if(IN!=null)
				{
				LOG.info("opening stdin "+this.IN);
				r=IOUtils.openFileForBufferedReading(this.IN);
				}
			else
				{
				LOG.info("opening stdin");
				r=new BufferedReader(new InputStreamReader(System.in));
				}
			
			while((line=r.readLine())!=null)
				{
				if(line.startsWith("#"))
					{
					System.out.println(line);
					continue;
					}
				boolean found=false;
				String tokens[]=tab.split(line);
				int pos0=Integer.parseInt(tokens[1]);
			    IntervalTree<KnownGene> kgs=kgMap.debugGetTree(tokens[0]);
			    if(kgs==null)
			    	{
			    	LOG.info("no gene found in chromosome "+tokens[0]+" (check chrom prefix?)");
			    	}
			    else
					{
			    	KnownGene bestGene=null;
					
					for(Iterator<IntervalTree.Node<KnownGene>> iter=kgs.iterator();iter.hasNext();)
						{
						KnownGene kg=iter.next().getValue();
						if(bestGene==null|| Math.abs(distance(pos0,kg))< Math.abs(distance(pos0,bestGene)))
							{
							bestGene=kg;
							}
						}
					if(bestGene!=null)
						{
						//get all transcritpts
						Collection<KnownGene> overlapKg=kgMap.getOverlapping(new Interval(
								bestGene.getChromosome(),
								bestGene.getTxStart()+1,
								bestGene.getTxEnd()
								));
						for(KnownGene kg:overlapKg)
							{
							KnownGene.Exon bestExon=null;
							for(KnownGene.Exon exon:kg.getExons())
								{
								if(bestExon==null || Math.abs(distance(pos0, exon))< Math.abs(distance(pos0,bestExon)))
									{
									bestExon=exon;
									}
								}
							if(bestExon!=null)
								{
								System.out.println(
										line+"\t"+
										kg.getName()+"\t"+
										kg.getTxStart()+"\t"+kg.getTxEnd()+"\t"+kg.getStrand()+"\t"+
										bestExon.getName()+"\t"+bestExon.getStart()+"\t"+bestExon.getEnd()+"\t"+
										distance(pos0,bestExon)
										);
								found=true;
								}
							}
						}
					
					}
				if(!found)
					{
					System.out.println(line+"\tNULL");
					}
				}
			r.close();
			LOG.info("done");
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar81455().instanceMain(args);
		}
	}
