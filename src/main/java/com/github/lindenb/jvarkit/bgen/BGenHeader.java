package com.github.lindenb.jvarkit.bgen;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.tribble.TribbleException;

public class BGenHeader extends BGenUtils {
	public enum Compression {
		e_NoCompression,
		e_ZlibCompression,
		e_ZstdCompression
		};
		
	public enum Layout {
		e_Layout1,
		e_Layout2
		};
    private List<String> samples;
    private Map<String,Integer> sample2idx;
    private BitSet headerFlags;
    private long num_variants_blocks;
    
    public Compression getCompression() {
    	final int op = getCompressedSNPBlocks();
    	switch(op) {
			case 0: return Compression.e_NoCompression;
			case 1: return Compression.e_ZlibCompression;
			case 2: return Compression.e_ZstdCompression;
			default: throw new IllegalArgumentException("bad compression ("+op+")");
    		}
    	}
    
    private int getCompressedSNPBlocks() {
        return (this.headerFlags.get(1)?2:0)+
               (this.headerFlags.get(0)?1:0);
        }
    
    public long getNVariants() {
    	return num_variants_blocks;
    	}
    
    public boolean hasAnonymousSamples() {
    	return  !this.headerFlags.get(31);
    	}
    
   public Layout getLayout() {
    	final int op= getLayoutOpCode();
    	switch(op) {
			case 1: return Layout.e_Layout1;
			case 2: return Layout.e_Layout2;
			default: throw new IllegalArgumentException("undefined layout ("+op+")");
			}
    	}
    
   private int getLayoutOpCode() {
        return (this.headerFlags.get(5)?8:0)+
                (this.headerFlags.get(4)?4:0)+
                (this.headerFlags.get(3)?2:0)+
                (this.headerFlags.get(2)?1:0);
        }
    
    public int getSampleIndex(final String sn) {
    	return sn==null?-1: sample2idx.getOrDefault(sn, -1);
    	}
    
    public List<String> getSamples() {
    	return this.samples;
    	}
    public int getNSamples() {
    	return this.samples.size();
    	}
    public String toString() {
    	return "BGenHeaderImpl";
    	}
    
   private static List<String> readSampleIdentifierBlock(BinaryCodec binaryCodec,int expect_n_samples)throws IOException {
	   /* sample_block_size */ binaryCodec.readUInt();        
       final int  number_of_samples = longToUnsignedInt(binaryCodec.readUInt());
       
       if(number_of_samples!=expect_n_samples) {
    	throw new IOException("expected "+expect_n_samples+" samples but got "+number_of_samples);   
       	}
       final List<String> samples= new ArrayList<>(number_of_samples);
       for(int i=0;i< number_of_samples;++i) {
    	   final String sn = readStringUInt16(binaryCodec);
           samples.add(sn);
           }      
       return Collections.unmodifiableList(samples);
   		}
    
  public  static BGenHeader of(final BinaryCodec binaryCodec) throws IOException {
	   final BGenHeader header=new BGenHeader();
	   
       final int header_block_size =  longToUnsignedInt(binaryCodec.readUInt());
       
       header.num_variants_blocks =  binaryCodec.readUInt();
      
       
       final int num_samples =  longToUnsignedInt(binaryCodec.readUInt());

      
       
       final byte magic[]=new byte[4];
       binaryCodec.readBytes(magic);
       if(!(Arrays.equals(magic, BGenUtils.BGEN_MAGIC) ||Arrays.equals(magic, new byte[] {0,0,0,0}))) {
           throw new TribbleException("Cannot decode magic "+magic[0]+"/"+magic[1]+"/"+magic[2]+"/"+magic[3]);
           }
      
       binaryCodec.getInputStream().skipNBytes(header_block_size-20);
       //Free data area. This could be used to store, for example, identifying information about the file
       //byte[] free_area = new byte[];
       //binaryCodec.readBytes(free_area);
       //free_area = null;
       
       // A set of flags, with bits numbered as for an unsigned integer. See below for flag definitions.
       byte flags[]=new byte[4];
       binaryCodec.readBytes(flags);
       header.headerFlags = BitSet.valueOf(flags);
     
   
       
       if( !header.hasAnonymousSamples()) {
           /*final int sample_block_size = unsignedIntToInt(binaryCodec.readUInt());
           if(Log.isEnabled(LogLevel.DEBUG)) {
        	   LOG.debug("sample_block_size :"+sample_block_size);
           		}*/
           
           header.samples= readSampleIdentifierBlock(binaryCodec,num_samples);
           /*
           try(BufferOutputStream bos = new BufferOutputStream(sample_block_size)) {
	           bos.copyFully(stream, sample_block_size);
	           }*/
           }
       else {
    	   
       		header.samples = new java.util.AbstractList<String>() {
	       		@Override
	       		public int size() { return num_samples;}
	       		@Override
	       		public String get(int index) {
	       			return String.format("anonymous_sample_%06d",(index+1));
	       			}
	       		};
       	}
       final Map<String,Integer> hash = new HashMap<>(header.samples.size());
       for(String sn:header.samples) {
       	if(hash.containsKey(sn)) {
       		throw new TribbleException("duplicate sample "+sn);
       		}
       	hash.put(sn, hash.size());
        header.sample2idx = Collections.unmodifiableMap(hash);
       	}
       
       return header;
   		}
    
   

}
