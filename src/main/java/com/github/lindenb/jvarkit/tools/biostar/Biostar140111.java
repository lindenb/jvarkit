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


History:
* 2015 creation 

*/
package com.github.lindenb.jvarkit.tools.biostar;

import gov.nih.nlm.ncbi.dbsnp.gt.SnpInfo;
import gov.nih.nlm.ncbi.dbsnp.gt.SnpInfo.SnpLoc;
import gov.nih.nlm.ncbi.dbsnp.gt.SnpInfo.SsInfo;
import gov.nih.nlm.ncbi.dbsnp.gt.SnpInfo.SsInfo.ByPop;
import gov.nih.nlm.ncbi.dbsnp.gt.SnpInfo.SsInfo.ByPop.GTypeByInd;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

@Program(name="biostar140111",
	description="How to obtain human genotype data from dpSNP ftp?",
	biostars=140111
	)
public class Biostar140111 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar130456.class).make();
	
	@SuppressWarnings("unused")
	private static final gov.nih.nlm.ncbi.dbsnp.gt.ObjectFactory _fool_javac1=null;	
	/** transforms XML/DOM to GBC entry */
	private Unmarshaller unmarshaller;
	


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	/** parses NCBI GT output */
	private void parseGenotypes(XMLEventReader r,OutputStream outStream)  throws Exception
		{
		final Pattern dnaRegex=Pattern.compile("[ATGCatgc]+");
		final QName attInId=new QName("indId");
		Set<String> samples=new TreeSet<>();
		Set<VCFHeaderLine> metaData = new HashSet<>();
		metaData.add(new VCFFormatHeaderLine(
				"GT",
				1,
				VCFHeaderLineType.String,
				"Genotype"));
		metaData.add(new VCFFormatHeaderLine(
				"DP",
				1,
				VCFHeaderLineType.Integer,
				"Depth"));
		
		metaData.add(new VCFInfoHeaderLine(
				"SnpInfoObserved",
				1,
				VCFHeaderLineType.String,
				"SnpInfo : Oberved"));

		
		VCFHeader header=null;
		VariantContextWriter w = null;
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
			if(evt.isStartElement()  )
				{
				StartElement startE=evt.asStartElement();
				String localName= startE.getName().getLocalPart();
				if(localName.equals("SnpInfo"))
					{
					if(header==null)
						{
						header=new VCFHeader(metaData, samples);
						w = VCFUtils.createVariantContextWriterToOutputStream(outStream);
						w.writeHeader(header);
						}
					SnpInfo snpInfo = this.unmarshaller.unmarshal(r, SnpInfo.class).getValue();
					
					if(snpInfo.getSnpLoc().isEmpty())
						{
						LOG.warning("no snploc for rs"+snpInfo.getRsId());
						}
					
						
					for(SnpLoc snpLoc:snpInfo.getSnpLoc())
						{
						int chromStart;
						try {
							chromStart = Integer.parseInt(snpLoc.getStart());
						} catch (Exception e) {
							LOG.warning("bad start in rs"+snpInfo.getRsId()+" "+snpLoc.getStart());
							continue;
							}
						chromStart++;
						
						
						
						String contigAllele=snpLoc.getContigAllele();
						if(contigAllele==null || !dnaRegex.matcher(contigAllele).matches())
							{
							LOG.warning("bad contigAllele in rs"+snpInfo.getRsId()+" "+contigAllele);
							continue;
							}
						if(!"fwd".equals(snpLoc.getRsOrientToChrom()))
							{
							contigAllele = AcidNucleics.reverseComplement(contigAllele);
							}
						Allele ref= Allele.create(contigAllele, true);
						Map<String,Genotype> sample2genotype=new HashMap<>();
						for(SsInfo ssinfo:snpInfo.getSsInfo())
							{
							boolean revcomp = !"fwd".equals(snpLoc.getRsOrientToChrom());
							if(!"fwd".equals(ssinfo.getSsOrientToRs())) revcomp=!revcomp;
							for(ByPop byPop:ssinfo.getByPop())
								{
								for(GTypeByInd gt:byPop.getGTypeByInd())
									{
									String sample=String.valueOf(gt.getIndId());
									if(!samples.contains(sample))
										{
										LOG.warning("Undefined sample:"+sample);
										continue;
										}
									boolean ok=true;
									String tokens[]= gt.getGtype().split("[/]");
									if(tokens.length==1)
										{
										tokens=new String[]{tokens[0],tokens[0]};
										}
									else if(tokens.length!=2)
										{
										LOG.warning("Bad genotypes in sample:"+sample+" "+gt.getGtype());
										continue;
										}
									List<Allele> sampleAlleles = new ArrayList<>(2);
									for(int i=0;i< tokens.length;++i)
										{
										if(revcomp) tokens[i]=AcidNucleics.reverseComplement(tokens[i]);
										if( !dnaRegex.matcher(tokens[i]).matches())
											{
											ok=false;
											break;
											}
										sampleAlleles.add(tokens[i].equalsIgnoreCase(contigAllele)?
												ref:
												Allele.create(tokens[i],false)
												);
										}
									if(!ok) continue;
									
									
									GenotypeBuilder gb=new GenotypeBuilder(sample,sampleAlleles);
									sample2genotype.put(sample, gb.make());
									}
								}
							}
						Set<Allele> alleles=new HashSet<>();
						alleles.add(ref);

						for(String sample:samples)
							{
							if(!sample2genotype.containsKey(sample))
								{
								sample2genotype.put(sample, GenotypeBuilder.createMissing(sample, 2));
								}
							else
								{
								alleles.addAll(sample2genotype.get(sample).getAlleles());
								}
							}

						
						VariantContextBuilder vcb=new VariantContextBuilder("dbsnp",
								snpLoc.getChrom(),
								chromStart,
								chromStart+ref.getBaseString().length()-1,
								alleles
								);
						if(snpInfo.getObserved()!=null)
							{
							vcb.attribute("SnpInfoObserved", VCFUtils.escapeInfoField(snpInfo.getObserved()));
							}

						vcb.genotypes(sample2genotype.values());
						vcb.id("rs"+snpInfo.getRsId());
						w.add(vcb.make());
						}
						
					}
				else if(localName.equals("Individual"))
					{
					
					if(header!=null) throw new XMLStreamException(
							"Error got "+localName+" after genotypes",evt.getLocation());
					Attribute  att = startE.getAttributeByName(attInId);
					if(att==null ) throw new XMLStreamException(
							"Cannot get "+attInId,evt.getLocation());
					samples.add(att.getValue());
					
					r.next();//consumme
					}
				else
					{
					r.next();
					}
				}
			else
				{
				r.next();//consumme
				}
			}
		if(w==null) throw new IOException("No Genotype was found");
		w.close();
		}
	
	@Override
	public int doWork(final List<String> args) {
		InputStream inputStream=null;
		XMLEventReader r=null;
		OutputStream pw=null;
		try
			{
			
			
			//create a Unmarshaller for genbank
			JAXBContext jc = JAXBContext.newInstance(
					"gov.nih.nlm.ncbi.dbsnp.gt");
			this.unmarshaller=jc.createUnmarshaller();
	
			XMLInputFactory xif=XMLInputFactory.newFactory();
			xif.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String publicID,
						String systemID, String baseURI, String namespace)
						throws XMLStreamException {
							return new ByteArrayInputStream(new byte[0]);
						}
				});
			
			if(args.isEmpty())
				{
				r= xif.createXMLEventReader(stdin(), "UTF-8");
				}
			
			else if(args.size()==1)
				{
				inputStream  = IOUtils.openURIForReading(args.get(0));
				r= xif.createXMLEventReader(inputStream);
				}
			else
				{
				LOG.error("Illegal number of arguments.");
				return -1;
				}
			
			pw = super.openFileOrStdoutAsStream(this.outputFile);
			
			this.parseGenotypes(r,pw);
			pw.flush();
			return 0;
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(this.outputFile!=null) CloserUtil.close(pw);
			CloserUtil.close(inputStream);
			CloserUtil.close(r);
			}
		}
	


	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar140111().instanceMainWithExit(args);
		}

	}
