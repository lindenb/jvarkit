package com.github.lindenb.jvarkit.tools.springbatch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import org.mortbay.io.RuntimeIOException;
import org.springframework.batch.core.StepExecution;
import org.springframework.batch.core.annotation.AfterStep;
import org.springframework.batch.item.ItemWriter;

import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextSplitWriter implements ItemWriter<List<VariantContext>>
	{
	
	
	private class VcfCache
		{
		final String key;
		private long last_time_insert=System.currentTimeMillis();
		private PrintWriter writer=null;
		private VCFEncoder encoder;
		VcfCache(final String key,final VCFHeader header) throws IOException {
			this.key=key;
			this.encoder=new VCFEncoder(header, false, false);
			this.writer = new PrintWriter(getFile());
			VCFUtils.convertVCFHeaderToList(header).
				stream().
				forEach(L->writer.println(L))
				;
			}
		File getFile() { return new File(this.key+".vcf");}
		boolean isOpen() { return this.writer!=null;}
		VcfCache reOpen() throws IOException {
			if(this.writer==null)
				{
				this.writer=new PrintWriter(new FileWriter(getFile(), true));
				}
			return this;
			}
		VcfCache close() throws IOException  {
				if(this.writer!=null)
					{
					this.writer.flush();
					this.writer.close();
					this.writer=null;
					}
				return this;
				}
		VcfCache add(VariantContext ctx) {
			this.writer.println(encoder.encode(ctx));
			this.last_time_insert = System.currentTimeMillis();;
			return this;}
			}
	private final Map<String,VcfCache> key2vcf = new HashMap<>();
	private Function<VariantContext, Set<String>> variant2keys = V->Collections.singleton(V.getContig());
	private int max_open=50;
	
	
	
	@Override
	public void write(final List<? extends List<VariantContext>> variants)
			throws Exception {
		variants.stream().flatMap(L->L.stream()).forEach(V->add(V));
		}
	private void add(final VariantContext ctx) {
		final Set<String> keys = this.variant2keys.apply(ctx);
		if(keys==null || keys.isEmpty()) return;
		for(final String key:keys) 
			{
			if(StringUtil.isBlank(key)) continue;
			try {
				getVcfCacheForKey(key).add(ctx);
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}
	private VcfCache getVcfCacheForKey(final String key) throws IOException{
		VcfCache vcfCache= this.key2vcf.get(key);
		if(vcfCache==null)
			{
			removeOneInCache();
			final VCFHeader header=null;//TODO
			vcfCache = new VcfCache(key,header);
			this.key2vcf.put(key,vcfCache);
			return vcfCache;
			}
		if(vcfCache.isOpen()) return vcfCache;
		removeOneInCache();
		return vcfCache.reOpen();
		}
	
	private void removeOneInCache()throws IOException
		{
		if(key2vcf.values().stream().filter(V->V.isOpen()).count() > this.max_open) {
			key2vcf.values().stream().filter(V->V.isOpen()).
				sorted((V1,V2)-> Long.compare(V1.last_time_insert,V2.last_time_insert)).
				findFirst().
				get().
				close();
			}
		}
	@AfterStep
	public void saveStepExecution(final StepExecution stepExecution) {
		this.key2vcf.values().stream().forEach(V->{try {V.close();}catch(Exception err){}});;
		stepExecution.getExecutionContext().putString("file", "out");
		}
	}
