package com.github.lindenb.jvarkit.tools.ws.impl;

import java.util.HashSet;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;

import com.github.lindenb.jvarkit.tools.ws.WSProject;

@XmlRootElement(name=WSProjectImpl.xmlName)
@XmlAccessorType(XmlAccessType.PROPERTY)
public class WSProjectImpl implements WSProject
	{
	private static final long serialVersionUID = 1L;
	public static final String xmlName="project";
	private String id;
	private String label;
	private String description;
	private Set<String> bamIds=null;
	private Set<String> vcfIds=null;
	
	
	public WSProjectImpl()
		{
		this(null,null,null,null,null);
		}
	
	public WSProjectImpl(String id, String label, String description,
			Set<String> bamIds,
			Set<String> vcfIds
			)
		{
		super();
		this.id = id;
		this.label = label;
		this.description = description;
		this.bamIds=bamIds;
		this.vcfIds=vcfIds;
		}

	@Override
	public String getId() {
		return id;
	}

	@Override
	public String getLabel() {
		return label;
	}

	@Override
	public String getDescription() {
		return description;
	}

	public void setId(String id) {
		this.id = id;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public void setDescription(String description) {
		this.description = description;
	}
	@Override
	@XmlElementWrapper(name = "bams")
	@XmlElement(name = "bam-id")
	public Set<String> getBamIds() {
		if(bamIds==null) bamIds=new HashSet<String>();
		return bamIds;
	}

	public void setBamIds(Set<String> bamIds) {
		this.bamIds = bamIds;
	}
	@Override
	@XmlElementWrapper(name = "vcfs")
	@XmlElement(name = "vcf-id")
	public Set<String> getVcfIds() {
		if(vcfIds==null) vcfIds=new HashSet<String>();
		return vcfIds;
	}

	public void setVcfIds(Set<String> vcfIds) {
		this.vcfIds = vcfIds;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		WSProjectImpl other = (WSProjectImpl) obj;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}

	}
