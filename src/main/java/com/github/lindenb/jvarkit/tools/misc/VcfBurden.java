/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;



public class VcfBurden extends AbstractKnimeApplication
	{
	private Map<String,Boolean> gene2seen=null;
	private String dirName="burden";
	public VcfBurden()
		{
		}
	
	
	private static class GeneTranscript
		{
		String geneName;
		String transcriptName;
		GeneTranscript(String geneName,String transcriptName)
			{
			this.geneName = geneName;
			this.transcriptName = transcriptName;
			}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((geneName == null) ? 0 : geneName.hashCode());
			result = prime * result
					+ ((transcriptName == null) ? 0 : transcriptName.hashCode());
			return result;
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			GeneTranscript other = (GeneTranscript) obj;
			if (geneName == null) {
				if (other.geneName != null)
					return false;
			} else if (!geneName.equals(other.geneName))
				return false;
			if (transcriptName == null) {
				if (other.transcriptName != null)
					return false;
			} else if (!transcriptName.equals(other.transcriptName))
				return false;
			return true;
			}
		
		}

	@Override
	public String getProgramDescription() {
		return "Solena: vcf to (chrom/pos/ref/alt/individus(G:0/1/2/-9)";
		}
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"SolenaVcfToRaw";
		}
	
	@Override
	public int initializeKnime() {
		return super.initializeKnime();
		}

	@Override
	public void disposeKnime() {
		super.disposeKnime();
		}
	
	
	@Override
	public int executeKnime(List<String> args)
		{
		ZipOutputStream zout=null;
		FileOutputStream fout=null;
		VcfIterator in=null;
		try
			{
			if(args.isEmpty())
				{
				info("reading from stdin.");
				in = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				String filename=args.get(0);
				info("reading from "+filename);
				in = VCFUtils.createVcfIterator(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			if(getOutputFile()==null)
				{
				error("undefined output");
				return -1;
				}
			else if(!getOutputFile().getName().endsWith(".zip"))
				{
				error("output "+getOutputFile()+" should end with .zip");
				return -1;
				}
			else
				{
				fout = new FileOutputStream(getOutputFile());
				zout = new ZipOutputStream(fout);
				}
			List<String> samples= in.getHeader().getSampleNamesInOrder();
			VCFHeader header=in.getHeader();
			String prev_chrom = null;
			VepPredictionParser vepPredParser=new VepPredictionParser(header);
			Map<GeneTranscript,List<VariantContext>> gene2variants=new HashMap<>();
			SequenceOntologyTree soTree= SequenceOntologyTree.getInstance();
			Set<SequenceOntologyTree.Term> acn=new HashSet<>();
			for(String acns:new String[]{
					"SO:0001589", "SO:0001587", "SO:0001582", "SO:0001583",
					"SO:0001575", "SO:0001578", "SO:0001574", "SO:0001889",
					"SO:0001821", "SO:0001822", "SO:0001893"
					})
				{
				acn.addAll(soTree.getTermByAcn(acns).getAllDescendants());
				}
			
		
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader());
			for(;;)
				{
				VariantContext ctx1=null;
				if(in.hasNext())
					{
					ctx1= progress.watch(in.next());
					if(ctx1.getAlternateAlleles().size()!=1)
						{
						//info("count(ALT)!=1 in "+ctx1.getChr()+":"+ctx1.getStart());
						continue;
						}
					}
				
				if(ctx1==null || !ctx1.getContig().equals(prev_chrom))
					{
					info("DUMP to zip n="+gene2variants.size());
					for(GeneTranscript gene_transcript:gene2variants.keySet() )
						{
						ZipEntry ze = new ZipEntry(this.dirName+"/"+gene_transcript.geneName+"_"+gene_transcript.transcriptName+".txt");
						zout.putNextEntry(ze);
						PrintWriter pw = new PrintWriter(zout);
						pw.print("CHROM\tPOS\tREF\tALT");
						for(String sample:samples)
							{
							pw.print("\t");
							pw.print(sample);
							}
						pw.println();
						for(VariantContext ctx:gene2variants.get(gene_transcript))
							{
							pw.print(ctx.getContig());
							pw.print("\t");
							pw.print(ctx.getStart());
							pw.print("\t");
							pw.print(ctx.getReference().getDisplayString());
							pw.print("\t");
							pw.print(ctx.getAlternateAlleles().get(0).getDisplayString());
							for(String sample:samples)
								{
								Genotype g=ctx.getGenotype(sample);
								pw.print("\t");
								if(g.isHomRef())
									{
									pw.print("0");
									}
								else if(g.isHomVar())
									{
									pw.print("2");
									}
								else if(g.isHet())
									{
									pw.print("1");
									}
								else
									{
									pw.print("-9");
									}
								}
							pw.println();
							}
						pw.flush();
						zout.closeEntry();
						}
					
					if(ctx1==null) break;
					gene2variants.clear();
					prev_chrom = ctx1.getContig();
					}
				Set<GeneTranscript> seen_names=new HashSet<>();
				for(VepPredictionParser.VepPrediction pred: vepPredParser.getPredictions(ctx1))
					{
					String geneName= pred.getSymbol();
					if(geneName==null || geneName.trim().isEmpty()) continue;
					
					
					if(this.gene2seen!=null)
						{
						if(!this.gene2seen.containsKey(geneName)) continue;
						
						}
					
					
					String transcriptName = pred.getFeature();
					if(transcriptName==null || transcriptName.trim().isEmpty())
						{
						info("No transcript in "+ctx1);
						continue;
						}
					
					GeneTranscript geneTranscript = new GeneTranscript(geneName, transcriptName);
					
					if(seen_names.contains(geneTranscript)) continue;
					boolean ok=false;
					for(SequenceOntologyTree.Term so:pred.getSOTerms())
						{
						if(acn.contains(so))
							{
							ok=true;
							}
						}
					if(!ok) continue;
					
					List<VariantContext> L = gene2variants.get(geneTranscript);
					if(L==null)
						{
						L=new ArrayList<>();
						gene2variants.put(geneTranscript,L);
						}
					L.add(ctx1);
					seen_names.add(geneTranscript);
					if(this.gene2seen!=null)
						{
						this.gene2seen.put(geneTranscript.geneName, Boolean.TRUE);
						}
					}
				}
			
			if(this.gene2seen!=null)
				{
				for(String gene:this.gene2seen.keySet())
					{
					if(this.gene2seen.get(gene).equals(Boolean.TRUE)) continue;
					warning("Gene not found : "+gene);
					ZipEntry ze = new ZipEntry(this.dirName+"/"+gene+"_000000000000000.txt");
					zout.putNextEntry(ze);
					PrintWriter pw = new PrintWriter(zout);
					pw.print("CHROM\tPOS\tREF\tALT");
					for(String sample:samples)
						{
						pw.print("\t");
						pw.print(sample);
						}
					pw.println();
					pw.flush();
					zout.closeEntry();
					}
				}
			
			progress.finish();
			
			zout.finish();
			fout.flush();
			zout.flush();
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(zout);
			CloserUtil.close(fout);
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (file) output file (default stdout)");
		out.println(" -g (file) optional list of gene names (restrict genes, print genes without data)");
		out.println(" -d (dir) base zip dir default:"+this.dirName);
		super.printOptions(out);
		}
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:d:g:"))!=-1)
			{
			switch(c)
				{
				case 'o': setOutputFile(opt.getOptArg()); break;
				case 'd': this.dirName=opt.getOptArg();break;
				case 'g': 
					{
						BufferedReader in=null;
					try {
						if(this.gene2seen==null) this.gene2seen=new HashMap<>();
						in = IOUtils.openURIForBufferedReading(opt.getOptArg());
						String line;
						while((line=in.readLine())!=null)
							{
							line=line.trim();
							if(line.isEmpty()||line.startsWith("#")) continue;
							this.gene2seen.put(line, Boolean.FALSE);
							}
						}
					catch (Exception e) {
						error(e);
						return -1;
						}
					finally
						{
						CloserUtil.close(in);
						}
					break;
					}
					
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new VcfBurden().instanceMainWithExit(args);
		}
	}
