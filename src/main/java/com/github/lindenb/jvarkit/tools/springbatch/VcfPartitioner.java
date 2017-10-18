package com.github.lindenb.jvarkit.tools.springbatch;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.core.partition.support.Partitioner;
import org.springframework.batch.item.ExecutionContext;

public class VcfPartitioner implements Partitioner {
private static final Log LOG = LogFactory.getLog(VcfPartitioner.class);    
@Override
public Map<String, ExecutionContext> partition(final int gridSize) {
	LOG.info("creating partition for gridsize="+gridSize);
	final Map<String, ExecutionContext> map = new HashMap<>(gridSize);
	for(int i=0;i<10;i++)
		{
		final ExecutionContext exec= new ExecutionContext();
		exec.put("contig", i);
		map.put(String.valueOf(i), exec);
		
		}
	return map;
	}
}
