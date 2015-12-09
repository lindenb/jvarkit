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
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
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
	private boolean highdamage=false;
	private Map<String,Boolean> _gene2seen=null;
	private String dirName="burden";
	private boolean printSOTerms=false;//mail matilde 11/10/2015 11:38 AM
	public VcfBurden()
		{
		}
	
	private static class VariantAndCsq
		{
		final VariantContext variant;
		final Set<SequenceOntologyTree.Term> terms;
		VariantAndCsq(VariantContext variant,Set<SequenceOntologyTree.Term> terms)
			{
			this.variant = variant;
			this.terms = new HashSet<>(terms);
			}
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
	
	
	private void dump(
			ZipOutputStream zout,
			final String contig,
			String filename,
			List<String> samples,
			List<VariantAndCsq> variants
			) throws IOException
		{
		info(this.dirName+"/"+filename+".txt");
		ZipEntry ze = new ZipEntry(this.dirName+"/"+contig+"_"+filename+".txt");
		zout.putNextEntry(ze);
		PrintWriter pw = new PrintWriter(zout);
		pw.print("CHROM\tPOS\tREF\tALT");
		if(printSOTerms)
			{
			pw.print("\t");
			pw.print("SO");//TODO
			}
		for(String sample:samples)
			{
			pw.print("\t");
			pw.print(sample);
			}
		pw.println();
		for(VariantAndCsq vac:variants)
			{
			pw.print(vac.variant.getContig());
			pw.print("\t");
			pw.print(vac.variant.getStart());
			pw.print("\t");
			pw.print(vac.variant.getReference().getDisplayString());
			pw.print("\t");
			pw.print(vac.variant.getAlternateAlleles().get(0).getDisplayString());
			if(printSOTerms)
				{
				pw.print("\t");
				for(SequenceOntologyTree.Term t:vac.terms)
					{
					pw.print(t.getAcn());
					pw.print(";");
					}
				}
			for(String sample:samples)
				{
				Genotype g=vac.variant.getGenotype(sample);
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
			Map<GeneTranscript,List<VariantAndCsq>> gene2variants=new HashMap<>();
			final SequenceOntologyTree soTree= SequenceOntologyTree.getInstance();
			final Set<SequenceOntologyTree.Term> acn=new HashSet<>();
			/* mail solena  *SO en remplacement des SO actuels (VEP HIGH + MODERATE) - pas la peine de faire retourner les analyses mais servira pour les futures analyses burden */
			String acn_list[]=new String[]{
					"SO:0001893",  "SO:0001574",  "SO:0001575", 
					"SO:0001587",  "SO:0001589",  "SO:0001578", 
					"SO:0002012",  "SO:0001889",  "SO:0001821", 
					"SO:0001822",  "SO:0001583",  "SO:0001818"
					
					/*
					"SO:0001589", "SO:0001587", "SO:0001582", "SO:0001583",
					"SO:0001575", "SO:0001578", "SO:0001574", "SO:0001889",
					"SO:0001821", "SO:0001822", "SO:0001893"*/
					};
			
			if(this.highdamage)
				{
				acn_list=new String[]{
						"SO:0001893",  "SO:0001574",  "SO:0001575",
						"SO:0001587",  "SO:0001589",  "SO:0001578", 
						"SO:0002012",  "SO:0001889"
						};
				}
			
			for(final String acns:acn_list)
				{
				final SequenceOntologyTree.Term  tacn = soTree.getTermByAcn(acns);
				if(tacn==null)
					{
					in.close();
					throw new NullPointerException("tacn == null pour "+acns);
					}
				acn.addAll(tacn.getAllDescendants());
				}
			
		
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader());
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
					Set<String> geneNames= new HashSet<>();
					for(GeneTranscript gene_transcript:gene2variants.keySet() )
						{
						geneNames.add(gene_transcript.geneName);
						dump(zout,
								prev_chrom,
								gene_transcript.geneName+"_"+gene_transcript.transcriptName,
								samples,
								gene2variants.get(gene_transcript)
								);
						}
					
					for(String geneName : geneNames)
						{
						Comparator<VariantAndCsq> cmp = new Comparator<VariantAndCsq>()
									{
									@Override
									public int compare(
											final VariantAndCsq v1,
											final VariantAndCsq v2) {
										final VariantContext o1 = v1.variant;
										final VariantContext o2 = v2.variant;
										int i = o1.getContig().compareTo(o2.getContig());
										if(i!=0) return i;
										i = o1.getStart() - o2.getStart();
										if(i!=0) return i;
										i = o1.getReference().compareTo(o2.getReference());
										if(i!=0) return i;
										i = o1.getAlternateAllele(0).compareTo(o2.getAlternateAllele(0));
										if(i!=0) return i;
										return 0;
										}
									};
						SortedSet<VariantAndCsq> lis_nm = new TreeSet<>(cmp);
						SortedSet<VariantAndCsq> lis_all = new TreeSet<>(cmp);
						SortedSet<VariantAndCsq> lis_refseq = new TreeSet<>(cmp);
						SortedSet<VariantAndCsq> lis_enst = new TreeSet<>(cmp);
						
						for(GeneTranscript gene_transcript:gene2variants.keySet() )
							{
							if(!geneName.equals(gene_transcript.geneName)) continue;
							lis_all.addAll(gene2variants.get(gene_transcript));
							if(gene_transcript.transcriptName.startsWith("NM_"))
								{
								lis_nm.addAll(gene2variants.get(gene_transcript));
								}
							if(! gene_transcript.transcriptName.startsWith("ENST"))
								{
								lis_refseq.addAll(gene2variants.get(gene_transcript));
								}
							if( gene_transcript.transcriptName.startsWith("ENST"))
								{
								lis_enst.addAll(gene2variants.get(gene_transcript));
								}
							}
						dump(zout,
								prev_chrom,
								geneName+"_ALL_TRANSCRIPTS",
								samples,
								new ArrayList<VariantAndCsq>(lis_all)
								);
						dump(zout,
								prev_chrom,
								geneName+"_ALL_NM",
								samples,
								new ArrayList<VariantAndCsq>(lis_nm)
								);
						dump(zout,
								prev_chrom,
								geneName+"_ALL_REFSEQ",
								samples,
								new ArrayList<VariantAndCsq>(lis_refseq)
								);
						dump(zout,
								prev_chrom,
								geneName+"_ALL_ENST",
								samples,
								new ArrayList<VariantAndCsq>(lis_enst)
								);
						}
					
					if(ctx1==null) break;
					gene2variants.clear();
					System.gc();
					prev_chrom = ctx1.getContig();
					}
				Set<GeneTranscript> seen_names=new HashSet<>();
				for(VepPredictionParser.VepPrediction pred: vepPredParser.getPredictions(ctx1))
					{
					
					String geneName= pred.getSymbol();
					if(geneName==null || geneName.trim().isEmpty()) continue;
					
					
					if(this._gene2seen!=null)
						{
						if(!this._gene2seen.containsKey(geneName)) continue;
						
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
					if(!printSOTerms && !ok) continue;
					
					List<VariantAndCsq> L = gene2variants.get(geneTranscript);
					if(L==null)
						{
						L=new ArrayList<>();
						gene2variants.put(geneTranscript,L);
						}
					L.add(new VariantAndCsq(ctx1,pred.getSOTerms()));
					seen_names.add(geneTranscript);
					if(this._gene2seen!=null)
						{
						this._gene2seen.put(geneTranscript.geneName, Boolean.TRUE);
						}
					}
				}
			
			if(this._gene2seen!=null)
				{
				final List<VariantAndCsq> emptylist = Collections.emptyList();
				for(String gene:this._gene2seen.keySet())
					{
					if(this._gene2seen.get(gene).equals(Boolean.TRUE)) continue;
					warning("Gene not found : "+gene);
					dump(zout,
							"UNDEFINED",
							gene+"_000000000000000.txt",
							samples,
							emptylist
							);
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
		out.println(" -H only high damage");
		out.println(" -f print ALL consequence+ colum, SO terms (mail matilde 11/10/2015 11:38 AM)");
		super.printOptions(out);
		}
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:d:g:Hf"))!=-1)
			{
			switch(c)
				{
				case 'f': this. printSOTerms = true; break;
				case 'H':  this.highdamage=true;break;
				case 'o': setOutputFile(opt.getOptArg()); break;
				case 'd': this.dirName=opt.getOptArg();break;
				case 'g': 
					{
						BufferedReader in=null;
					try {
						if(this._gene2seen==null) this._gene2seen=new HashMap<>();
						in = IOUtils.openURIForBufferedReading(opt.getOptArg());
						String line;
						while((line=in.readLine())!=null)
							{
							line=line.trim();
							if(line.isEmpty()||line.startsWith("#")) continue;
							this._gene2seen.put(line, Boolean.FALSE);
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
