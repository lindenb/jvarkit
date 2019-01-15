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
package com.github.lindenb.jvarkit.tools.evs2bed;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.XMLEvent;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import edu.washington.gs.evs.SnpData;

@Program(name="evs2vcf",description= "Download data from EVS http://evs.gs.washington.edu/EVS as a BED chrom/start/end/XML For later use, see VCFTabixml.")
public class EvsToVcf extends Launcher
	{
	private static final Logger LOG=Logger.build(EvsToVcf.class).make();

	private EvsToVcf()
		{
		
		}
	
	private void _fillDict(SAMSequenceDictionary dict,String name,int sequenceLength)
		{
		dict.addSequence(new SAMSequenceRecord(name, sequenceLength));
		}
	
	private boolean isEmpty(String s)
		{
		return s==null  || s.trim().isEmpty();
		}
	private Float parseDouble(String s)
		{
		if(isEmpty(s)) return null;
		try {
			Float d=new Float(s);
			return d;
		} catch (Exception e) {
			return null;
			}
		}
	@Override
	public int doWork(List<String> args) {
		
		VariantContextWriter out=null;
		try
			{
			if(!args.isEmpty())
				{
				LOG.error("Illegal number of arguments");
				return -1;
				}
			JAXBContext jc = JAXBContext.newInstance(SnpData.class);
	        Unmarshaller unmarshaller = jc.createUnmarshaller();		

			
            out=VCFUtils.createVariantContextWriterToStdout();
            SAMSequenceDictionary dict=new SAMSequenceDictionary();
            _fillDict(dict,"1",249250621);
            _fillDict(dict,"2",243199373);
            _fillDict(dict,"3",198022430);
            _fillDict(dict,"4",191154276);
            _fillDict(dict,"5",180915260);
            _fillDict(dict,"6",171115067);
            _fillDict(dict,"7",159138663);
            _fillDict(dict,"8",146364022);
            _fillDict(dict,"9",141213431);
            _fillDict(dict,"10",135534747);
            _fillDict(dict,"11",135006516);
            _fillDict(dict,"12",133851895);
            _fillDict(dict,"13",115169878);
            _fillDict(dict,"14",107349540);
            _fillDict(dict,"15",102531392);
            _fillDict(dict,"16",90354753);
            _fillDict(dict,"17",81195210);
            _fillDict(dict,"18",78077248);
            _fillDict(dict,"19",59128983);
            _fillDict(dict,"20",63025520);
            _fillDict(dict,"21",48129895);
            _fillDict(dict,"22",51304566);
            _fillDict(dict,"X",155270560);
            _fillDict(dict,"Y",59373566);
            _fillDict(dict,"MT",16569);

            VCFHeader header=new VCFHeader();
            header.setSequenceDictionary(dict);
            header.addMetaDataLine(new VCFInfoHeaderLine("CONS",VCFHeaderLineCount.INTEGER,VCFHeaderLineType.Float,"conservationScore"));
            header.addMetaDataLine(new VCFInfoHeaderLine("GERP",VCFHeaderLineCount.INTEGER,VCFHeaderLineType.Float,"conservationScoreGERP"));
            header.addMetaDataLine(new VCFInfoHeaderLine("uaMAF",VCFHeaderLineCount.INTEGER,VCFHeaderLineType.Float,"conservationScoreGERP"));
            header.addMetaDataLine(new VCFInfoHeaderLine("aaMAF",VCFHeaderLineCount.INTEGER,VCFHeaderLineType.Float,"conservationScoreGERP"));
            header.addMetaDataLine(new VCFInfoHeaderLine("totalMAF",VCFHeaderLineCount.INTEGER,VCFHeaderLineType.Float,"conservationScoreGERP"));
            header.addMetaDataLine(new VCFInfoHeaderLine("DP",VCFHeaderLineCount.INTEGER,VCFHeaderLineType.Integer,"conservationScoreGERP"));

            header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
            header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
            JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header);
            
            out.writeHeader(header);
            Pattern comma= Pattern.compile("[,]");
            XMLInputFactory xif=XMLInputFactory.newFactory();
            xif.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, false);
            XMLEventReader xmlr=xif.createXMLEventReader(System.in);
            while(xmlr.hasNext() && !System.out.checkError())
	        	{
		    	XMLEvent evt=xmlr.peek();
		    	if(!evt.isStartElement() ||
		    		!evt.asStartElement().getName().getLocalPart().equals("snpList"))
		    		{
		    		xmlr.nextEvent();
		    		continue;
		    		}
		    	SnpData snpData= unmarshaller.unmarshal(xmlr,SnpData.class).getValue();
		    	VariantContextBuilder vcb=new VariantContextBuilder();
		    	Set<Allele> alleles=new HashSet<Allele>();
		    	alleles.add(Allele.create(snpData.getRefAllele(), true));
		    	for(String s: comma.split(snpData.getAltAlleles()))
		    		{
		    		if(isEmpty(s)) continue;
		    		alleles.add(Allele.create(s,false));
		    		}
		    	
		    	vcb.chr(snpData.getChromosome());
		    	vcb.start(snpData.getChrPosition());
		    	vcb.stop(snpData.getChrPosition()+snpData.getRefAllele().length()-1);
		    	if(!isEmpty(snpData.getRsIds()) && !snpData.getRsIds().equals("none"))
			    	{
				    vcb.id(snpData.getRsIds());
			    	}
		    	vcb.alleles(alleles);
		    	Float d=parseDouble(snpData.getConservationScore());
		    	if(d!=null)
			    	{
			    	vcb.attribute("CONS", d);
			    	}
		    	d=parseDouble(snpData.getConservationScoreGERP());
		    	if(d!=null)
			    	{
			    	vcb.attribute("GERP", d);
			    	}
		    	
			    	
		    	vcb.attribute("uaMAF",(float)snpData.getUaMAF());
		    	vcb.attribute("aaMAF",(float)snpData.getAaMAF());
		    	vcb.attribute("totalMAF",(float)snpData.getTotalMAF());
		    	vcb.attribute("DP",snpData.getAvgSampleReadDepth());
			    	
		    	out.add(vcb.make());
		    	
	        	}
            xmlr.close();
            
            out.close();

			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}

	public static void main(String[] args) {
		new EvsToVcf().instanceMainWithExit(args);
		}
	}
