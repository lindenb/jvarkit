package com.github.lindenb.jvarkit.tools.springbatch;

import java.util.HashMap;
import java.util.Map;

import org.springframework.batch.core.partition.support.Partitioner;
import org.springframework.batch.item.ExecutionContext;

public class VcfPartitioner implements Partitioner {
@Override
public Map<String, ExecutionContext> partition(final int gridSize) {
	Map<String, ExecutionContext> map = new HashMap<>(gridSize);
	return map;
	}
}
