# KnimeVariantHelper

A java library to be used in the java nodes of [http://knime.org](http://knime.org). This library allow to use the htsjdk library into knime.


## Requirements

* Tested with knime 3.3.2
* Check you're running with java 1.8 ( http://www.oracle.com/technetwork/java/javase/downloads/index.html ) .

```
Knime -> Help -> About Knime Analytics Platform -> Installation details -> Configuration : 

java.version=1.8.*

```


## Download or Compile:

A version of the library might be available at: [https://github.com/lindenb/jvarkit/releases][https://github.com/lindenb/jvarkit/releases).

Compiling:

```
$ make knimehelper
```

will generate a jar file in `dist/knimehelper.jar`

## In Knime

* create a new Node `java Filter`
* in the tab 'additional libraries', add 'knimehelper.jar'.

## See also:

* https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/knime/KnimeVariantHelper.java
* https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java


## Examples

# configuration Knime.ini

```
-Dhttp.proxyHost=cache.ha...
-Dhttps.proxyHost=cache.ha...
-Dhttp.proxyPort=3128
-Dhttps.proxyPort=3128
-Dhttp.nonProxyHosts=IP1,IP2,IP3
-Duser.language=en
-Duser.country=US
-Duser.variant=EN
```

dans File->Prefs->General->network


(Manual)
Http/Https: cache.ha...
Fill: ProxybyPass



## Example

Julien B. 2017-05-31 (16 samples, two families  ). **Java Snippet Row Filter**: 

In Tab **Additional Libararies** , add `knimehelper.jar`.


In Tab **Java Snippet**  :

Section **Global Variable Déclaration** :


```java
/** KnimeVariantHelper : it's a bridge between the htsjdk java library for HTS data and the tables in knime */
final com.github.lindenb.jvarkit.knime.KnimeVariantHelper helper = new com.github.lindenb.jvarkit.knime.KnimeVariantHelper();

/** an associative map between the family names their sample names */
final Map<String,Set<String>> fam2samples = new TreeMap<>();
```

Section **Method Body** :

```java
try	{
	/** fam2samples is empty: this is the first time we're scanning the table: let's initialize fam2samples by filling the associative map 'family-name' -> sample-set */
	if( fam2samples.isEmpty())
		{
		/* family FAM1 */
		fam2samples.put("FAM1", new java.util.HashSet<>(java.util.Arrays.asList("S1","S2","S3")));
		/* family FAM2 */
		fam2samples.put("FAM2", new java.util.HashSet<>(java.util.Arrays.asList("S4","S5","S6")));
		}

	/** re-build a java object 'VariantContext' ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html ) from the table.  */

	htsjdk.variant.variantcontext.VariantContext variant = helper.build().
		contig($#CHROM$).
		pos($POS$).
		id($ID$).
		ref($REF$).
		alts($ALT$).
		qual($QUAL$).
		filter($FILTER$).
		info($INFO$).
		format($FORMAT$).
		genotype("S1",$S1$).
		genotype("S2",$S2$).
		genotype("S3",$S3$).
		genotype("S4",$S4$).
		genotype("S5",$S5$).
		genotype("S6",$S6$).
		genotype("S7",$S7$).
		genotype("S8",$S8$).
		genotype("S9",$S9$).
		genotype("S10",$S10$).
		genotype("S11",$S11$).
		genotype("S12",$S12$).
		genotype("S13",$S13$).
		genotype("S14",$S14$).
		genotype("S15",$S15$).
		genotype("S16",$S16$).build();
	
	/** reject the variant if it's not annotated with Sequence Ontology term http://www.sequenceontology.org/miso/current_svn/term/SO:0001818 */
	if( !helper.hasSequenceOntologyLabel(variant,"protein_altering_variant") ) return false;


	/** look at a few attribute from gnomad */
	for(final String tag: new String[]{"gnomad_genome_AF_NFE","gnomad_exome_AF_NFE"})
		{
		/* if there is any frequency that is greater than 0.1 , reject the variant */
		if( variant.getAttributeAsStringList(tag,"").stream().
			       filter(S->!(S.isEmpty()) || S.equals(".")).
			       map(S->Double.parseDouble(S)).
			       filter(AF->AF>0.1).findAny().isPresent() ) return false;
		}


	/** should we keep the variant  */
	boolean keep=false;

	/** loop over the  familles */
	for(final String family: fam2samples.keySet())
		{
		/* number of samples IN the current family carrying a variant */
		int count_in = 0; 
		/* number of samples OUT of the current family carrying a variant */
		int count_out = 0;
		/** affected: samples in the current family */
		final Set<String> affected =fam2samples.get(family);

		/** loop over all the genotypes in the variant */
		for(int j=0;j< variant.getNSamples();++j)
			{
			/* get the j-th genotype */
			htsjdk.variant.variantcontext.Genotype genotype = variant.getGenotype(j);
			/* ignore this genotype if it's  '0/0' or './.' */
			if(genotype.isHomRef() || genotype.isNoCall()) continue;
			/* the sample linked to the genotype belongs to the current family */
			if(affected.contains(genotype.getSampleName()))
				{
				count_in++;
				}
			else /* the sample linked to the genotype doesn't belong to the current family */
				{
				count_out++;
				}
			}
		/* keep the variant if count_in == number-of-individuals-in-the-family and if there is not individual affected out of the current_family */
		if(count_in == affected.size() && count_out==0)
			{
			keep = true;
			}
		}

	/** last line ? cleanup things. */
	if( $$ROWINDEX$$ +1 == $$ROWCOUNT$$ )
		{
		helper.dispose();
		}

	return keep;	
	} 
catch(final Throwable err)
	{
	System.err.println("################### ERROR with "+ $#CHROM$ +" "+$POS$);
	err.printStackTrace();
	return false;
	}
```




## Example

```
 try {
	final com.github.lindenb.jvarkit.knime.KnimeVariantHelper helper = new com.github.lindenb.jvarkit.knime.KnimeVariantHelper();
	
	htsjdk.variant.variantcontext.VariantContext ctx = helper.build().
	   contig($#CHROM$).
	   pos($POS$).
	   id($ID$).
	   ref($REF$).
	   alts($ALT$).
	   filter($FILTER$).
	   info($INFO$).
	   format($FORMAT$).
	   genotype("S1",$S1$).
	   genotype("S2",$S2$).
	   build()
	   ;
	
	helper.initSnpEffParser("Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|RefSeq|SIFT|PolyPhen");
	
	
	
	return  ctx.hasID()==false  && ctx.isIndel() &&
	        helper.hasSequenceOntologyLabel(ctx,"protein_altering_variant")  && 
	        ctx.getAlternateAlleles().size()==1 && 
			ctx.getAttributeAsDouble("AF",1.0) > 0.0001  &&
	        !ctx.getGenotype("S1").sameGenotype( ctx.getGenotype("S2") )
	    ;
	}
catch(Throwable err)
	{
	System.err.println("################### ERROR with "+ $#CHROM$ +" "+$POS$);
	err.printStackTrace();
	return false;
	}
	 
```

## generating the script from the VCF header

![http://i.imgur.com/3SO6Mic.png](http://i.imgur.com/3SO6Mic.png)


* Node Read File :VCF
* Extract Column Header;: prefix column
* Unpivoting: Value Columns: (all), Retained columns (nothing/empty)
* JavasSnippet: append String column 'js'

```java
String s="";

switch($$ROWINDEX$$)
	{
	case 0: s="htsjdk.variant.variantcontext.VariantContext ctx = helper.build().contig($"+ $ColumnValues$ +"$)."; break;
	case 1: s="pos($"+ $ColumnValues$ +"$)."; break;
	case 2: s="id($"+ $ColumnValues$ +"$)."; break;
	case 3: s="ref($"+ $ColumnValues$ +"$)."; break;
	case 4: s="alts($"+ $ColumnValues$ +"$)."; break;
	case 5: s="qual($"+ $ColumnValues$ +"$)."; break;
	case 6: s="filter($"+ $ColumnValues$ +"$)."; break;
	case 7: s="info($"+ $ColumnValues$ +"$)."; break;
        case 8: s="format($"+ $ColumnValues$ +"$)."; break;
	default: s= "genotype(\""+$ColumnValues$ +"\",$" +$ColumnValues$ +"$).";if( $$ROWINDEX$$ +1 == $$ROWCOUNT$$) s+="build();"; break;
	}

return s;
```

* ColumnFilter: keep 'js' column
* CSV Writer: overwrite, quote Mode: never, no header, no row-id


