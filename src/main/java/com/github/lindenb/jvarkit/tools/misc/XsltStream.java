package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.xml.XMLConstants;
import javax.xml.namespace.QName;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.Namespace;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

/**
BEGIN_DOC

## Example:

Dumping the Orcid from pubmed:

```
 java  -jar dist/pubmeddump.jar 'orcid[AUID]' |\
 	java -jar dist/xsltstream.jar -t pubmed2orcid.xsl -n "PubmedArticle" 
```

The XSLT stylesheet:

```xml
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
    xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
    version='1.0'
	>

<xsl:output method="text" />


<xsl:template match="/">
<xsl:apply-templates select="PubmedArticle"/>
</xsl:template>

<xsl:template match="PubmedArticle">
<xsl:apply-templates select="MedlineCitation/Article/AuthorList/Author[Identifier/@Source='ORCID']"/>
</xsl:template>

<xsl:template match="Author">
<xsl:value-of select="LastName"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="ForeName"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="Initials"/>
<xsl:text>	</xsl:text>
<xsl:call-template name="orcid"><xsl:with-param name="s" select="Identifier[@Source='ORCID']"/></xsl:call-template>
<xsl:text>	</xsl:text>
<xsl:for-each select="Affiliation"><xsl:text> </xsl:text></xsl:for-each>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../../PMID[1]"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../../DateCreated/Year"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../Journal/ISOAbbreviation"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../ArticleTitle"/>
<xsl:text>
</xsl:text>
</xsl:template>

<xsl:template name="orcid">
<xsl:param name="s"/>
<xsl:choose>
	<xsl:when test="starts-with($s,'http://orcid.org/')">
		<xsl:call-template name="orcid">
			<xsl:with-param name="s" select="substring($s,18)"/>
		</xsl:call-template>
	</xsl:when>
	<xsl:when test="starts-with($s,'https://orcid.org/')">
		<xsl:call-template name="orcid">
			<xsl:with-param name="s" select="substring($s,19)"/>
		</xsl:call-template>
	</xsl:when>
	<xsl:when test="starts-with($s,'https://')">
		<xsl:call-template name="orcid">
			<xsl:with-param name="s" select="substring($s,9)"/>
		</xsl:call-template>
	</xsl:when>
	<xsl:otherwise>
		<xsl:value-of select="translate($s,'-','')"/>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

</xsl:stylesheet>
```

output:

```
Kerkis	I	I	0000000344337580		28618452	2017	Cell Prolif.	Murine melanoma cells incomplete reprogramming using non-viral vector.
Zhang	Shuijun	S	0000000205993289		28618450	2017	Cell Prolif.	SAV1 represses the development of human colorectal cancer by regulating the Akt-mTOR pathway in a YAP-dependent manner.
Nguyen	Ha Trong	HT	0000000222408942		28618448	2017	Health Econ	Out of sight but not out of mind: Home countries' macroeconomic volatilities and immigrants' mental health.
Lee	Jeongmi	J	0000000299487554		28618213	2017	J Sep Sci	Solid-phase-extraction-assisted dispersive liquid-liquid microextraction based on solidification of floating organic droplet to determine sildenafil and its analogues in dietary supplements.
Kwon	Sung Won	SW	0000000171614737		28618213	2017	J Sep Sci	Solid-phase-extraction-assisted dispersive liquid-liquid microextraction based on solidification of floating organic droplet to determine sildenafil and its analogues in dietary supplements.
Villaverde	Juan J	JJ	000000025911792X		28618212	2017	Pest Manag. Sci.	Quantum chemistry in environmental pesticide risk assessment.
Pollard	Thomas D	TD	0000000217852969		28618211	2017	Cytoskeleton (Hoboken)	Tribute to Fumio Oosawa the pioneer in actin biophysics.
Xiao	Bingxiu	B	0000000285929251		28618205	2017	J. Clin. Lab. Anal.	Reduced expression of circRNA hsa_circ_0003159 in gastric cancer and its clinical significance.
Heal	M Elisabeth	ME	0000000150571141		28618202	2017	Congenit Heart Dis	Effects of persistent Fontan fenestration patency on cardiopulmonary exercise testing variables.
Hrubec	Terry C	TC	0000000239619201		28618200	2017	Birth Defects Res	Ambient and dosed exposure to quaternary ammonium disinfectants causes neural tube defects in rodents.
Somri	Mostafa	M	0000000238141402		28618198	2017	Int J Paediatr Dent	Effect of intravenous paracetamol as pre-emptive compared to preventive analgesia in a pediatric dental setting: a prospective randomized study.
Haggblom	Max M	MM	0000000163077863		28618195	2017	Environ Microbiol Rep	Novel Reductive Dehalogenases from the Marine Sponge Associated Bacterium Desulfoluna spongiiphila.
Ehl	Stefan	S	0000000162861234		28618194	2017	Insect Sci.	Sexual dimorphism in the alpine butterflies Boloria pales and Boloria napaea: Differences in movement and foraging behaviour (Lepidoptera: Nymphalidae).
Gautam	Nischal K	NK	0000000224916705		28618193	2017	Paediatr Anaesth	Introduction of color-flow injection test to confirm intravascular location of peripherally placed intravenous catheters.
Kaymaz	Dicle	D	0000000179512065		28618190	2017	Clin Respir J	RELATION BETWEEN UPPER-LIMB MUSCLE STRENGTH WITH EXERCISE CAPACITY, QUALITY OF LIFE, AND DYSPNEA IN PATIENTS WITH SEVERE CHRONIC OBSTRUCTIVE PULMONARY DISEASE.
Erill	Ivan	I	0000000272807191		28618189	2017	Environ. Microbiol.	Comparative genomics of the DNA-damage inducible network in the Patescibacteria.
Ii	Satoshi	S	0000000254285385		28618187	2017	Int J Numer Method Biomed Eng	Physically consistent data assimilation method based on feedback control for patient-specific blood flow analysis.
Arzi	Boaz	B	0000000272898994		28618186	2017	Stem Cells Transl Med	Therapeutic Efficacy of Fresh, Allogeneic Mesenchymal Stem Cells for Severe Refractory Feline Chronic Gingivostomatitis.
Clark	Kaitlin C	KC	0000000260959382		28618186	2017	Stem Cells Transl Med	Therapeutic Efficacy of Fresh, Allogeneic Mesenchymal Stem Cells for Severe Refractory Feline Chronic Gingivostomatitis.
Friedrich	Anja	A	0000000297356286		28618185	2017	J Sleep Res	Let's talk about sleep: a systematic review of psychological interventions to improve sleep in college students.
Lee	Yun Hee	YH	0000000150273988		28618180	2017	Clin Respir J	Neutrophil-lymphocyte ratio and a dosimetric factor for predicting symptomatic radiation pneumonitis in non-small-cell lung cancer patients treated with concurrent chemoradiotherapy.
Dhatariya	Ketan	K	0000000336199579		28618177	2017	Int. J. Clin. Pract.	Assessing the quality of primary care referrals to surgery of patients with diabetes in the East of England: A multi-centre cross-sectional cohort study.
Tougeron	Kevin	K	0000000348973787		28618174	2017	Insect Sci.	Intraspecific maternal competition induces summer diapause in insect parasitoids.
De Paepe	Kim	K	0000000279486765		28618173	2017	Environ. Microbiol.	Inter-individual differences determine the outcome of wheat bran colonization by the human gut microbiome.
Daria	Dzema	D	0000000194181022		28618171	2017	J Sep Sci	Highly fluorinated polymers with sulfonate, sulfamide and N,N-diethylamino groups for the capillary electromigration separation of proteins and steroid hormones.
Tinoco	Adelita	A	0000000221905824		28618169	2017	Ann Noninvasive Electrocardiol	ECG-derived Cheyne-Stokes respiration and periodic breathing in healthy and hospitalized populations.
Doyle	Zelda	Z	0000000186481383		28618161	2017	Aust J Rural Health	Prevention of osteoporotic refractures in regional Australia.
Locker	Jacomine Krijnse	JK	0000000186582977		28618160	2017	Cell. Microbiol.	VACCINIA VIRUS A11 IS REQUIRED FOR MEMBRANE RUPTURE AND VIRAL MEMBRANE ASSEMBLY.
Davis	Adam S	AS	0000000271961197		28618159	2017	Pest Manag. Sci.	Are herbicides a once in a century method of weed control?
Nelson	C E	CE	0000000325253496		28618153	2017	Environ. Microbiol.	Cascading influence of inorganic nitrogen sources on DOM production, composition, lability and microbial community structure in the open ocean.
(...)
```

END_DOC
*/
@Program(name="xsltstream",
	description="XSLT transformation for large XML files. xslt is only applied on a given subset of nodes.",
	keywords={"xml","xslt","xsl","stylesheet"}
	)
public class XsltStream extends Launcher {
	private static final Logger LOG = Logger.build(XsltStream.class).make();

	private static class QNameConverter implements IStringConverter<QName>
		{
		@Override
		public QName convert( String s) {
			String ns=null;
			String prefix =null;
			String localName=null;
			if(s.startsWith("{"))
				{
				int par=s.indexOf('}');
				if(par==-1) throw new IllegalArgumentException("closing '}' missing in '"+s+"'.");
				ns=s.substring(1, par);
				s=s.substring(par+1);
				}
			int colon = s.lastIndexOf(':');
			if(colon==-1)
				{
				localName=s;
				if(ns==null || ns.trim().isEmpty()) return new QName(localName);
				return new QName(ns, localName);
				}
			else
				{
				prefix=s.substring(0,colon);
				localName=s.substring(colon+1);
				if(ns==null || ns.trim().isEmpty()) return new QName(XMLConstants.NULL_NS_URI,localName,prefix);
				return new QName(ns, localName,prefix);
				}
			}
		}
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-t","-template"},description="XSLT template file.",required=true)
	private File templateFile = null;
	@Parameter(names={"-n","--tag","--name","-tag","-name"},description="XML node name. name has syntax '{ns}prefix:localName' or 'prefix:localName' or 'localName' or '{ns}localName' ",converter=QNameConverter.class)
	private Set<QName> targetQnames = new HashSet<>();
	@Parameter(names={"-skip","--skip"},description="Ignore those names",converter=QNameConverter.class)
	private Set<QName> skipQNames = new HashSet<>();

	
	private static String toQNAME(final QName q)
		{
		return
			(q.getPrefix()==null|| q.getPrefix().isEmpty()?"":q.getPrefix()+":")+q.getLocalPart();
		}
	
	private void skip(final StartElement startE,final XMLEventReader r)throws XMLStreamException
		{
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement())
				{
				skip(evt.asStartElement(),r);
				}
			else if(evt.isEndElement())
				{
				return;
				}
			}
			
		}
	
	private Node parseDocument(final Document owner,final StartElement startE,final XMLEventReader r) throws XMLStreamException
		{
		final QName qname=startE.getName();
		final Element root;
		
		if(qname.getNamespaceURI()==null || qname.getNamespaceURI().isEmpty())
			{
			root = owner.createElement(toQNAME(qname));
			}
		else
			{
			root = owner.createElementNS(
					qname.getNamespaceURI(),
					toQNAME(qname)
					);
			}
		Iterator<?> it= startE.getNamespaces();
		while(it.hasNext())
			{
			final Namespace att = Namespace.class.cast(it.next());
			root.setAttribute("xmlns:"+att.getPrefix(),att.getNamespaceURI());
			}
		
		it= startE.getAttributes();
		while(it.hasNext())
			{
			final Attribute att = Attribute.class.cast(it.next());
			final QName attName=att.getName();
			if(attName.getNamespaceURI()==null || attName.getNamespaceURI().isEmpty())
				{
				root.setAttribute(toQNAME(attName),att.getValue());
				}
			else
				{	
				root.setAttributeNS(
						attName.getNamespaceURI(), 
						toQNAME(attName),
						att.getValue());
				}
			}
		
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			
			if(evt.isCharacters())
				{
				root.appendChild(owner.createTextNode(evt.asCharacters().getData()));
				}
			else if(evt.isProcessingInstruction())
				{
				}
			else if(evt.isEndElement())
				{
				break;
				}
			else if(evt.isStartElement())
				{
				if(this.skipQNames.contains(evt.asStartElement().getName()))
					{
					skip(evt.asStartElement(),r);
					}
				else
					{
					root.appendChild(parseDocument(owner,evt.asStartElement(),r));
					}
				}
			else 
				{
				LOG.warn("Cannot handle "+evt);
				}
			}
		return root;
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.targetQnames.isEmpty())
			{
			LOG.error("No target name defined");
			return -1;
			}
		XMLEventReader xmlReader=null;
		OutputStream outputStream = null;
		try {
			final String inputSource = oneFileOrNull(args);
			final TransformerFactory transformerFactory = TransformerFactory.newInstance();
			final Templates template = transformerFactory.newTemplates(new StreamSource(this.templateFile));
			final Transformer transformer = template.newTransformer();
			
			outputStream = openFileOrStdoutAsStream(this.outputFile);
			
			final XMLInputFactory xif = XMLInputFactory.newFactory();
			xmlReader = (inputSource==null?
					xif.createXMLEventReader(stdin()):
					xif.createXMLEventReader(new StreamSource(new File(inputSource)))	
					);
			
			final Document dom = DocumentBuilderFactory.newInstance().
					newDocumentBuilder().
					newDocument();
			
			while(xmlReader.hasNext())
				{
				final XMLEvent evt = xmlReader.nextEvent();
				if(!evt.isStartElement())
					{
					continue;
					}
				final StartElement SE = evt.asStartElement();
				if(!this.targetQnames.contains(SE.getName()))
					{
					continue;
					}
				
				dom.appendChild(parseDocument(dom,SE,xmlReader));
				transformer.transform(
						new DOMSource(dom),
						new StreamResult(outputStream)
						);
				
				while(dom.hasChildNodes()) dom.removeChild(dom.getFirstChild());
				}
			outputStream.flush();
			outputStream.close();
			outputStream=null;
			xmlReader.close();
			xmlReader=null;
			return 0;
			}
		catch (final Exception e) {
			LOG.error(e);
			return -1;
			}
		finally {
			CloserUtil.close(xmlReader);
			CloserUtil.close(outputStream);
			}
		}
	
	public static void main(final String[] args) {
		new XsltStream().instanceMainWithExit(args);
	}
}
