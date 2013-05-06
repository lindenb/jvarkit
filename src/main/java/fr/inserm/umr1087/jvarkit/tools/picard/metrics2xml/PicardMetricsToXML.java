package fr.inserm.umr1087.jvarkit.tools.picard.metrics2xml;

import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import net.sf.picard.metrics.Header;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Histogram;

public class PicardMetricsToXML
	{
	private static final String NS="http://picard.sourceforge.net/";
	private static class MinMax
		{
		Double min=null;
		Double max=null;
		double sum=0.0;
		int count=0;
		}
	private boolean print_sums=false;
	private MetricsFile<MetricBase, Comparable<?>> metricsFile;
	private XMLStreamWriter out;
	
	private MetricsFile<MetricBase, Comparable<?>> getFile() {
		return metricsFile;
		}
	
	private void printHeaders() throws XMLStreamException
		{
		out.writeStartElement("headers");
		for (final Header h : getFile().getHeaders())
			{
			out.writeStartElement("header");
            out.writeAttribute("class", h.getClass().getName());
            out.writeCharacters(h.toString());
            out.writeEndElement();
			}
		out.writeEndElement();
		}
	
	private void printBeanMetrics() throws XMLStreamException
		{
		out.writeStartElement("metrics");
		if(getFile().getMetrics()!=null && !getFile().getMetrics().isEmpty())
			{
			Class<?> beanType=getFile().getMetrics().get(0).getClass();
			Field fields[]=beanType.getFields();
			List<MinMax> minmaxList=new ArrayList<MinMax>(fields.length);
			for(int i=0;i< fields.length;++i) minmaxList.add(new MinMax());
			out.writeStartElement("thead");
			out.writeAttribute("class",beanType.getName());
			for(Field field:fields )
				{
				out.writeStartElement("th");
				out.writeAttribute("class",field.getType().getName());
				out.writeCharacters(field.getName());
				out.writeEndElement();//th
				}
			out.writeEndElement();//thead
			out.writeStartElement("tbody");
			for(MetricBase bean:getFile().getMetrics())
				{
				out.writeStartElement("tr");
				for(int i=0; i<fields.length; ++i)
					{
					try{
						final Object value = fields[i].get(bean);
						
						if(value==null)
							{
							out.writeEmptyElement("td");
							out.writeAttribute("xsi",XMLConstants.W3C_XML_SCHEMA_INSTANCE_NS_URI,"nil","true");
							}
						else
							{
							if(print_sums && value instanceof Number)
								{
								double numeric=Number.class.cast(value).doubleValue();
								MinMax minmax=minmaxList.get(i);
								if(!(Double.isNaN(numeric) || Double.isInfinite(numeric)))
									{
									if(minmax.count==0)
										{
										minmax.max=numeric;
										minmax.min=numeric;										}
									else
										{
										minmax.max=Math.max(numeric,minmax.max);
										minmax.min=Math.min(numeric,minmax.min);
										}
									minmax.sum+=numeric;
									minmax.count++;
									}
								}
							out.writeStartElement("td");
							out.writeCharacters(String.valueOf(value));
							out.writeEndElement();//td
							}
						}
					catch (IllegalAccessException e)
						{
						out.writeEmptyElement("td");
						out.writeAttribute("error",e.getClass().getSimpleName());
						}
					
					}
				out.writeEndElement();//tr
				}
			out.writeEndElement();//tbody
			
			if(print_sums)
				{
				out.writeStartElement("tfoot");
				for(int i=0; i<fields.length; ++i)
					{
					MinMax minmax=minmaxList.get(i);
					out.writeEmptyElement("measurement");
					out.writeAttribute("for", fields[i].getName());
					if(minmax.count>0)
						{
						out.writeAttribute("count", String.valueOf(minmax.count));
						out.writeAttribute("min", String.valueOf(minmax.min));
						out.writeAttribute("max", String.valueOf(minmax.max));
						out.writeAttribute("mean", String.valueOf(minmax.sum/minmax.count));
						}		
					}
				out.writeEndElement();//tfoot
				}
			}	
		out.writeEndElement();
		}
	@SuppressWarnings({"unchecked","rawtypes"})
	private void printHistogram()throws XMLStreamException
		{
		
		final List<Histogram<?>> nonEmptyHistograms = new ArrayList<Histogram<?>>();
		for(Histogram<?> histo: getFile().getAllHistograms())
			{
			if(histo!=null && !histo.isEmpty()) nonEmptyHistograms.add(histo);
			}
		if(nonEmptyHistograms.isEmpty()) return;
		SortedSet keys = new TreeSet(nonEmptyHistograms.get(0).comparator());
		for(Histogram<?> histo:nonEmptyHistograms)
			{
			keys.addAll(histo.keySet());
			}
		out.writeStartElement("histogram");
		out.writeAttribute("class",  nonEmptyHistograms.get(0).keySet().iterator().next().getClass().getName());
		out.writeStartElement("thead");
		out.writeStartElement("th");
		out.writeCharacters(String.valueOf(nonEmptyHistograms.get(0).getBinLabel()));
		out.writeEndElement();//th
		for (final Histogram<?> histo : nonEmptyHistograms)
		 	{
			out.writeStartElement("th");
			out.writeCharacters(String.valueOf(histo.getValueLabel()));
			out.writeEndElement();//th
		 	}
		out.writeEndElement();//thead
		out.writeStartElement("tbody");
		for(Object key:keys)
			{
			out.writeStartElement("tr");
			out.writeStartElement("td");
			out.writeCharacters(String.valueOf(key));
			out.writeEndElement();//td
			for (final Histogram histo : nonEmptyHistograms)
				{
				Histogram.Bin bin = (Histogram.Bin)histo.get(key);
				if(bin==null)
					{
					out.writeEmptyElement("td");
					out.writeAttribute("xsi",XMLConstants.W3C_XML_SCHEMA_INSTANCE_NS_URI,"nil","true");
					}
				else
					{
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(bin.getValue()));
					out.writeEndElement();//td
					}
				}
			out.writeEndElement();//tr
			}
		out.writeEndElement();
		out.writeEndElement();
		}
	
	private void parse(String filename,Reader in) throws IOException,XMLStreamException
		{
		out.writeStartElement("metrics-file");
		out.writeAttribute("file", filename);
		this.metricsFile = new MetricsFile<MetricBase, Comparable<?>>();
		this.metricsFile.read(in);
		printHeaders();
		printBeanMetrics();
		printHistogram();
		out.writeEndElement();
		}
	
	private void parse(File file) throws IOException,XMLStreamException
		{
		parse(file.toString(),new FileReader(file));
		}
	public void run(String[] args) throws IOException,XMLStreamException
		{
		
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				return;
				}
			else if(args[optind].equals("-s"))
				{
				print_sums=true;
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return;
				}
			else
				{
				break;
				}
			++optind;
			}
		XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		this.out= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
		this.out.setDefaultNamespace(NS);
		out.writeStartDocument("UTF-8","1.0");
		out.writeStartElement("picard-metrics");
		out.writeDefaultNamespace(NS);
		out.writeAttribute(XMLConstants.XMLNS_ATTRIBUTE, XMLConstants.XML_NS_URI,"xsi", XMLConstants.W3C_XML_SCHEMA_INSTANCE_NS_URI);
		if(optind==args.length)
			{
			parse("stdin",new InputStreamReader(System.in));
			}
		else
			{
			while(optind< args.length)
				{
				String inputName=args[optind++];
				parse(new File(inputName));
				}
			}
		out.writeEndElement();
		out.writeEndDocument();
		out.flush();
		out.close();
		}
	public static void main(String[] args)
		throws IOException,XMLStreamException
		{
		new PicardMetricsToXML().run(args);
		}
	}
