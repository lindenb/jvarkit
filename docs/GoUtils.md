# GoUtils

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Gene Ontology Utils. Retrieves terms from Gene Ontology


## Usage

```
Usage: goutils [options] Files
  Options:
    -A, --accession
      User Go Terms accession numbers or name
      Default: []
    -af, --accession-file
      File containing accession numbers. One per line. After the first white 
      space one can define optional attributes for 
      gexf:`color=<COLOR>;size=<SIZE> 
    -action, --action
      What shoud I do ? default is dump as table. 'goa' only keeps GOA 
      elements in GOA input in GAF format (e.g 
      http://geneontology.org/gene-associations/goa_human.gaf.gz). 
      Default: dump_table
      Possible Values: [dump_table, dump_gexf, goa]
    -go, --go, --gene-ontology
      Gene ontology URI. Formatted as OWL+XML. Can be gzipped.
      Default: http://purl.obolibrary.org/obo/go.owl
    -go-divisions, --go-divisions
      limit the gene ontology tree to those divisions. empty: all possible 
      divisions. 
      Default: [biological_process, cellular_component, molecular_function]
    -go-relations, --go-relations
      limit the gene ontology tree to those relationships. empty: all possible 
      relationships. 
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --inverse
      inverse the result
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * geneontology
 * go
 * gexf



## See also in Biostars

 * [https://www.biostars.org/p/488538](https://www.biostars.org/p/488538)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew goutils
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20180130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/go/GoUtils.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/go/GoUtils.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **goutils** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

children of  GO:0005216 'ion channel activity' 

```
$ java -jar dist/goutils.jar -R is_a  -A 'GO:0005216' 

#ACN	NAME	DEFINITION
GO:1905030	voltage-gated ion channel activity involved in regulation of postsynaptic membrane potential	Any voltage-gated ion channel activity that is involved in regulation of postsynaptic membrane potential.
GO:1905057	voltage-gated calcium channel activity involved in regulation of postsynaptic cytosolic calcium levels	Any voltage-gated calcium channel activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:1905054	calcium-induced calcium release activity involved in regulation of presynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of presynaptic cytosolic calcium ion concentration.
GO:1905058	calcium-induced calcium release activity involved in regulation of postsynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:0016286	small conductance calcium-activated potassium channel activity	Enables the transmembrane transfer of potassium by a channel with a unit conductance of 2 to 20 picoSiemens that opens in response to stimulus by internal calcium ions. Small conductance calcium-activated potassium channels are more sensitive to calcium than are large conductance calcium-activated potassium channels. Transport by a channel involves catalysis of facilitated diffusion of a solute (by an energy-independent process) involving passage through a transmembrane aqueous pore or channel, without evidence for a carrier-mediated mechanism.
GO:0043855	cyclic nucleotide-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0043854	cyclic nucleotide-gated mechanosensitive ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens in response to a mechanical stress and when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0099142	intracellularly ATP-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when ATP has been bound by the channel complex or one of its constituent parts on the intracellular side of the plasma membrane.
GO:0099101	G-protein gated potassium channel activity	A potassium channel activity that is gated by binding of a G-protein beta-gamma dimer.
(...)
```

## Example

Use GO annotation to retrieve genes associated to GO:0005216 'ion channel activity' 

```
join -t $'\t' -1 1 -2 2 \
	<(java -jar dist/goutils.jar -A 'GO:0005216' | cut -f 1 | sort | uniq) \
	<(wget -q -O - "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD" | gunzip -c | grep -v '^!' | cut -f3,5 | uniq | LC_ALL=C sort -t $'\t' -k2,2) |\
sort -t $'\t' -k2,2 |\
grep SCN5A -A 10 -B 10
```

```
(...)
GO:0086006	SCN2B
GO:0005244	SCN3A
GO:0005248	SCN3A
GO:0005248	SCN3A
GO:0005248	SCN3B
GO:0086006	SCN3B
GO:0005248	SCN4A
GO:0005248	SCN4A
GO:0005248	SCN4B
GO:0086006	SCN4B
GO:0005244	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0086006	SCN5A
GO:0086060	SCN5A
GO:0086061	SCN5A
GO:0086062	SCN5A
GO:0086063	SCN5A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN8A
GO:0005248	SCN8A
GO:0005248	SCN9A
GO:0005248	SCN9A
GO:0005248	SCN9A
GO:0005272	SCNN1A
(...)
```

```
$ wget -q -O - "http://geneontology.org/gene-associations/goa_human.gaf.gz" | gunzip -c | java -jar dist/goutils.jar -go go.rdf.gz --action goa -A 'ion transmembrane transport'  | head
UniProtKB	A0A1W2PN81	CHRNA7	enables	GO:0022848	GO_REF:0000002	IEA	InterPro:IPR002394	FNeuronal acetylcholine receptor subunit alpha-7	CHRNA7	protein	taxon:9606	20210612	InterPro
UniProtKB	A0PJK1	SLC5A10	enables	GO:0015370	Reactome:R-HSA-8876283	TAS		F	Sodium/glucose cotransporter 5	SLC5A10|SGLT5	protein	taxon:9606	20200515	Reactome
UniProtKB	A0PJK1	SLC5A10	involved_in	GO:0035725	GO_REF:0000108	IEA	GO:0015370	P	Sodium/glucose cotransporter 5	SLC5A10|SGLT5	protein	taxon:9606	20210613	GOC
UniProtKB	A1A4F0	SLC66A1L	involved_in	GO:1903401	GO_REF:0000108	IEA	GO:0015189	PPutative uncharacterized protein SLC66A1L	SLC66A1L|C3orf55|PQLC2L	protein	taxon:9606	20210613	GOC
UniProtKB	A1A4F0	SLC66A1L	involved_in	GO:1903826	GO_REF:0000108	IEA	GO:0015181	PPutative uncharacterized protein SLC66A1L	SLC66A1L|C3orf55|PQLC2L	protein	taxon:9606	20210613	GOC
UniProtKB	A1A5B4	ANO9	enables	GO:0005229	PMID:22946059	IMP		F	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20140424	UniProt
UniProtKB	A1A5B4	ANO9	enables	GO:0005229	Reactome:R-HSA-2684901	TAS		F	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20200515	Reactome
UniProtKB	A1A5B4	ANO9	involved_in	GO:0034220	Reactome:R-HSA-983712	TAS		P	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20210310	Reactome
UniProtKB	A5X5Y0	HTR3E	enables	GO:0022850	PMID:17392525	IDA		F	5-hydroxytryptamine receptor 3E	HTR3E	protein	taxon:9606	20130210	CACAO
UniProtKB	A6NJY1	SLC9B1P1	enables	GO:0015299	GO_REF:0000002	IEA	InterPro:IPR006153	FPutative SLC9B1-like protein SLC9B1P1	SLC9B1P1	protein	taxon:9606	20210612	InterPro
```

