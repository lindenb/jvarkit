package com.github.lindenb.jvarkit.tools.ws.server;

import java.util.List;

import javax.jws.WebMethod;
import javax.jws.WebService;

@WebService
public interface NGSService
	{
	@WebMethod(operationName="listBams")
	public List<String> getBams();
	@WebMethod
	public List<String> getVCFs();
	@WebMethod
	public Bam getBam(String filename,String chrom,int start,int end);
	@WebMethod
	public Vcf getVcf(String filename,String chrom,int start,int end);
	}

