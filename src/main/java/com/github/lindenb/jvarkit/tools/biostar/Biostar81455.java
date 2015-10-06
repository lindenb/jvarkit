/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.biostar;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;


public class Biostar81455 extends AbstractBiostar81455
	{
	
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
    
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar78400.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar81455.AbstractBiostar81455Command
		{    

		@Override
			public Collection<Throwable> call(String inputFile) throws Exception {
			IntervalTreeMap<KnownGene> kgMap=new IntervalTreeMap<KnownGene>();
			BufferedReader r=null;
			PrintStream out=null;
			try
				{
				LOG.info("readubf "+kgUri);
				String line;
				Pattern tab=Pattern.compile("[\t]");
				r=IOUtils.openURIForBufferedReading(this.kgUri);
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
				
				
				if(inputFile!=null)
					{
					LOG.info("opening stdin "+inputFile);
					r=IOUtils.openURIForBufferedReading(inputFile);
					}
				else
					{
					LOG.info("opening stdin");
					r=new BufferedReader(new InputStreamReader(stdin()));
					}
				
				if(getOutputFile()==null)
					{
					out=stdout();
					}
				else
					{
					out= new PrintStream(getOutputFile());
					}
				
				while((line=r.readLine())!=null)
					{
					if(line.startsWith("#"))
						{
						out.println(line);
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
									out.println(
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
						out.println(line+"\tNULL");
						}
					}
				r.close();
				out.flush();
				LOG.info("done");
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(r);
				CloserUtil.close(out);
				}
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar81455().instanceMain(args);
		}
	}
