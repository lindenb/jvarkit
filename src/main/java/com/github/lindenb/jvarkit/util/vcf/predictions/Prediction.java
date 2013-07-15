package com.github.lindenb.jvarkit.util.vcf.predictions;

public interface Prediction {
public String getGeneName();
public String getEnsemblGene();
public String getEnsemblTranscript();
public String getEnsemblProtein();
public String getReferenceAminoAcid();
public Integer getAminoAcidPosition();
public String getAltAminoAcid();
}
