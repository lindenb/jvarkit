package com.github.lindenb.jvarkit.tools.ws.impl;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import com.github.lindenb.jvarkit.tools.ws.WSVcf;

@XmlRootElement(name="vcf")
@XmlAccessorType(XmlAccessType.PROPERTY)
public class WSVcfImpl implements WSVcf
	{
	private static final long serialVersionUID = 1L;
	private String id;
	private String path;
	private String referenceId=null;
	
	
	public WSVcfImpl()
		{
		this(null,null,null);
		}
	
	public WSVcfImpl(String id,String path,String referenceId)
		{
		this.id=id;
		this.path=path;
		this.referenceId=referenceId;
		}
	
	@Override
	@XmlAttribute(name="id",required=true)
	public String getId()
		{
		return id;
		}

	@Override
	public String getReferenceId()
		{
		return referenceId;
		}

	@Override
	public String getPath()
		{
		return path;
		}

	public void setId(String id) {
		this.id = id;
	}

	public void setPath(String path) {
		this.path = path;
	}

	public void setReferenceId(String referenceId) {
		this.referenceId = referenceId;
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
		WSVcfImpl other = (WSVcfImpl) obj;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "WSVcfImpl [id=" + id + ", path=" + path + ", referenceId="
				+ referenceId + "]";
	}
	
}
