# XsltStream

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

XSLT transformation for large XML files. xslt is only applied on a given subset of nodes.


## Usage

```
Usage: xsltstream [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -n, --tag, --name, -tag, -name
      XML node name. name has syntax '{ns}prefix:localName' or 
      'prefix:localName' or 'localName' or '{ns}localName'
      Default: []
    -o, --output
      Output file. Optional . Default: stdout
    -skip, --skip
      Ignore those names
      Default: []
    --version
      print version and exit
  * -t, -template
      XSLT template file.

```


## Keywords

 * xml
 * xslt
 * xsl
 * stylesheet



## See also in Biostars

 * [https://www.biostars.org/p/270498](https://www.biostars.org/p/270498)
 * [https://www.biostars.org/p/280581](https://www.biostars.org/p/280581)
 * [https://www.biostars.org/p/282545](https://www.biostars.org/p/282545)
 * [https://www.biostars.org/p/282602](https://www.biostars.org/p/282602)
 * [https://www.biostars.org/p/335867](https://www.biostars.org/p/335867)
 * [https://www.biostars.org/p/343432](https://www.biostars.org/p/343432)
 * [https://www.biostars.org/p/365479](https://www.biostars.org/p/365479)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew xsltstream
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/XsltStream.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/XsltStream.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **xsltstream** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


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

## Example:

Sample identifiers from NCBI biosamples


the xslt:

```xslt
<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0'>
<xsl:output method="text"  encoding="UTF-8"/>
<xsl:template match="BioSample">
<xsl:copy>
<xsl:apply-templates select="Ids"/>
</xsl:copy>
</xsl:template>

<xsl:template match="Ids">
<xsl:value-of select="Id[@db='BioSample']/text()"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="Id[@db='SRA']/text()"/>
<xsl:text>
</xsl:text>
</xsl:template>

</xsl:stylesheet>

```

execute:

```
curl -s "ftp://ftp.ncbi.nlm.nih.gov/biosample/biosample_set.xml.gz" |\
  gunzip -c |\
  java -jar dist/xsltstream.jar -n BioSample -t transform.xsl |\
  nl


(....)
7211789	SAMN07945461	SRS2643051
7211790	SAMN07945462	SRS2643052
7211791	SAMN07945463	SRS2643049
7211792	SAMN07945464	SRS2643050
7211793	SAMN07945465	
7211794	SAMN07945466	
7211795	SAMN07945467	
7211796	SAMN07945468	
```


## Example:

rs / ss list from dbsnp

```xsl
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0' xmlns:r="https://www.ncbi.nlm.nih.gov/SNP/docsum">
<xsl:output method="text"/>
<xsl:template match="/">
<xsl:apply-templates select="//r:Rs/r:Ss"/>
</xsl:template>

<xsl:template match="r:Ss">rs<xsl:value-of select="../@rsId"/> ss<xsl:value-of select="@ssId"/><xsl:text>
</xsl:text>
</xsl:template>
</xsl:stylesheet>
```

execute:

```
 wget -O - -q  "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/XML/ds_ch1.xml.gz" |\
  gunzip -c |\
  java -jar dist/xsltstream.jar -t transform.xsl --tag '{https://www.ncbi.nlm.nih.gov/SNP/docsum}Rs'
```

output:

```
rs171 ss41715810
rs171 ss43026199
rs171 ss96405203
rs242 ss242
rs242 ss287669350
rs242 ss326012704
rs242 ss326012781
rs242 ss498801024
rs242 ss550913725
rs242 ss552749651
(...)
```

## Example:

rough exploration of non-coding variants with pathogenic consequences in clinvar:

```xsl
<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml">
<xsl:output method="text"/>
<xsl:template match="/"><xsl:apply-templates/></xsl:template>
<xsl:template match="ReleaseSet"><xsl:apply-templates/></xsl:template>
<xsl:template match="ClinVarSet">
<xsl:if test="count(.//Attribute[@Type='MolecularConsequence'])=1">
<xsl:variable name="csq" select=".//Attribute[@Type='MolecularConsequence']/text()"/>
<xsl:if test="($csq = 'non-coding transcript variant' or $csq  = 'intergenic variant'  or $csq  = '2kb upstream variant'  or $csq  = '5 prime utr variant' ) and .//Description[contains(text(),'athogenic')]">
<xsl:value-of select="ReferenceClinVarAssertion/ClinVarAccession/@Acc"/>
<xsl:text>	</xsl:text>
<xsl:for-each select=".//Citation/ID[@Source = 'Pubmed']">pmid<xsl-value-of select="text()"/>;</xsl:for-each>
<xsl:text>
</xsl:text>
</xsl:if>
</xsl:if>
</xsl:template>
</xsl:stylesheet>
```

invoke:

```
$ curl -s  "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz" | gunzip -c |\
   java -jar dist/xsltstream.jar  -t transform.xsl -n ClinVarSet
   
RCV000203227	
RCV000256194	
RCV000256207	
RCV000000913	
RCV000000914	
RCV000006518
```

## Example

convert drugbank xml to TSV

```xslt
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet xmlns:d="http://www.drugbank.ca" xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0'>
<xsl:output method="text"/>

<xsl:template match="d:drugbank">
<xsl:apply-templates select="d:drug"/>
</xsl:template>

<xsl:template match="d:drug">
<xsl:value-of select="d:name/text()"/>
<xsl:text>	</xsl:text>
<xsl:for-each select="d:groups/d:group">
	<xsl:if test='position()>1'>-&gt;</xsl:if>
	<xsl:value-of select="./text()"/>
</xsl:for-each>
<xsl:text>	</xsl:text>
<xsl:for-each select="d:calculated-properties/d:property[d:kind/text()='InChIKey']/d:value">
	<xsl:if test='position()>1'> </xsl:if>
	<xsl:value-of select="./text()"/>
</xsl:for-each>
<xsl:text>	</xsl:text>
<xsl:for-each select="d:external-identifiers/d:external-identifier[d:resource/text()='ChEMBL']/d:identifier">
	<xsl:if test='position()>1'> </xsl:if>
	<xsl:value-of select="./text()"/>
</xsl:for-each>
<xsl:text>	</xsl:text>
<xsl:for-each select="d:external-identifiers/d:external-identifier[d:resource/text()='PubChem Compound']/d:identifier">
	<xsl:if test='position()>1'> </xsl:if>
	<xsl:value-of select="./text()"/>
</xsl:for-each>
<xsl:text>	</xsl:text>
<xsl:for-each select="d:external-identifiers/d:external-identifier[d:resource/text()='PubChem Substance']/d:identifier">
	<xsl:if test='position()>1'> </xsl:if>
	<xsl:value-of select="./text()"/>
</xsl:for-each>
<xsl:text>
</xsl:text>
</xsl:template>

</xsl:stylesheet>
```

and

```
$ java -jar dist/xsltstream.jar \
	-n '{http://www.drugbank.ca}drug' \
	-t drugbank2tsv.xsl \
	full_database.xml
```

## Example

https://www.biostars.org/p/335867/#335885

How to download database of Human protein sequences with sub cellular locations?

see https://gist.github.com/lindenb/b3c726adecde90e37acd92bc940dfdd5

## Example

https://www.biostars.org/p/365479/ "Bioinformatics word cloud to use in classes bioinformatics"

```xslt
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet  xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0' >
<xsl:output method="text" encoding="UTF-8"/>

<xsl:template match="/">
<xsl:apply-templates select="*"/>
</xsl:template>

<xsl:template match="*">
<xsl:apply-templates select="PubmedArticle"/>
</xsl:template>

<xsl:template match="PubmedArticle">
<xsl:variable name="year" select="MedlineCitation/Article/Journal/JournalIssue/PubDate/Year/text()"/>
<xsl:for-each select="MedlineCitation/MeshHeadingList/MeshHeading/DescriptorName">
<xsl:value-of select="$year"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="./text()"/>
<xsl:text>
</xsl:text>
</xsl:for-each>
</xsl:template>

</xsl:stylesheet>
```


see https://gist.github.com/lindenb/5d7773a93d8c2b0edbd4c01bf8834919


