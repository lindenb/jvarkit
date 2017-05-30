How to use the KnimeVariantHelper.

## Compile:

```
$ make knimehelper
```

## In Knime

* create a new Node `java Filter`
* in the tab 'additional libraries', add 'knimehelper.jar'.

## See also:

* https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/knime/KnimeVariantHelper.java
* https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java


### Example


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
	   genotype("CD12100",$CD12100$).
	   genotype("CD12218",$CD12218$).
	   build()
	   ;
	
	helper.initSnpEffParser("Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|RefSeq|SIFT|PolyPhen");
	
	
	
	return  ctx.hasID()==false  && ctx.isIndel() &&
	        helper.hasSequenceOntologyLabel(ctx,"protein_altering_variant")  && 
	        ctx.getAlternateAlleles().size()==1 && 
			ctx.getAttributeAsDouble("AF",1.0) > 0.0001  &&
	        !ctx.getGenotype("CD12100").sameGenotype( ctx.getGenotype("CD12218") )
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


