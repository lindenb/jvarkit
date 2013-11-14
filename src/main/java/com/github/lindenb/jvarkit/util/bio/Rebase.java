package com.github.lindenb.jvarkit.util.bio;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Rebase
	implements Iterable<Rebase.Enzyme>
	{
	private List<Rebase.Enzyme> enzymes=new ArrayList<Rebase.Enzyme>();
	public interface Enzyme
		{
		public String getName();
		public String getDecl();
		public String getBases();
		public boolean isPalindromic();
		public int size();
		public char at(int index);
		}
	
	public  Rebase.Enzyme getEnzymeByName(String name)
		{
		for(Rebase.Enzyme E:this)
			{
			if(E.getName().equals(name)) return E;
			}
		return null;
		}
	
	public Rebase.Enzyme get(int index)
		{
		return getEnzymes().get(index);
		}
	
	public int size()
		{
		return getEnzymes().size();
		}
	
	public List<Rebase.Enzyme> getEnzymes()
		{
		return enzymes;
		}
	
	@Override
	public Iterator<Enzyme> iterator()
		{
		return getEnzymes().iterator();
		}
	
	public static boolean compatible(
			char plasmid,
			char enzyme
			)
		{
		switch(plasmid)
			{
			case 'A': case 'a' : return "ANDHVMWR".indexOf(enzyme)!=-1;
            case 'T': case 't' : return "TNYBDHKW".indexOf(enzyme)!=-1;
            case 'G': case 'g' : return "GNRBDVKS".indexOf(enzyme)!=-1;
            case 'C': case 'c' : return "CNYBHVMS".indexOf(enzyme)!=-1;
            case 'N': case 'n' : return false;
			default: return false;
			}
		}

	
	
	public class EnzymeImpl implements Enzyme
		{
		private String name;
		private String seq;
		private String decl;
		private boolean palindromic;
		public EnzymeImpl(String name,String decl)
			{
			this.name=name;
			this.decl=decl;
			this.seq=decl.replaceAll("[^A-Za-z]","").toUpperCase();
			this.palindromic=decl.contains("^");
			}
		public String getName() { return this.name;}
		public String getBases() { return this.seq;}
		public String getDecl() { return this.decl;}
		public int size() { return this.seq.length();}
		public boolean isPalindromic() { return this.palindromic;}
		public char at(int index) { return this.seq.charAt(index); }
		@Override
		public boolean equals(Object obj)
			{
			if(this==obj) return true;
			if(obj==null || !(obj instanceof EnzymeImpl)) return false;
			return EnzymeImpl.class.cast(obj).name.equals(this.name);
			}
		@Override
		public int hashCode()
			{
			return name.hashCode();
			}
		@Override
		public String toString()
			{
			return name+"["+decl+"]";
			}
		}
	
	private void add(String name,String decl)
		{
		this.enzymes.add(new EnzymeImpl(name, decl));
		}
	public static Rebase createDefaultRebase()
		{
		Rebase rebase=new Rebase();
		/* curl -s "ftp://ftp.neb.com/pub/rebase/allenz.txt" |grep -E '^<[1257]>' | awk '/^<[^7]>/ {printf("%s\t",substr($0,4)); next;} /^<7>/ {print substr($0,4); next;}' | awk -F '       ' '($2=="" && $3!="?" && $4!="?" && $4!="")'  | cut  -f 1,3 | sort -k1,1 | awk '{printf("rebase.add(\"%s\",\"%s\");\n",$1,$2);}' */
		rebase.add("AarI","CACCTGC(4/8)");
		rebase.add("AatII","GACGT^C");
		rebase.add("AbsI","CC^TCGAGG");
		rebase.add("AccI","GT^MKAC");
		rebase.add("AciI","CCGC(-3/-1)");
		rebase.add("AclI","AA^CGTT");
		rebase.add("AcyI","GR^CGYC");
		rebase.add("AflII","C^TTAAG");
		rebase.add("AflIII","A^CRYGT");
		rebase.add("AgeI","A^CCGGT");
		rebase.add("AgsI","TTS^AA");
		rebase.add("AjuI","(7/12)GAANNNNNNNTTGG(11/6)");
		rebase.add("AlfI","(10/12)GCANNNNNNTGC(12/10)");
		rebase.add("AloI","(7/12)GAACNNNNNNTCC(12/7)");
		rebase.add("AluI","AG^CT");
		rebase.add("AlwNI","CAGNNN^CTG");
		rebase.add("ApaI","GGGCC^C");
		rebase.add("ApaLI","G^TGCAC");
		rebase.add("ApoI","R^AATTY");
		rebase.add("ArsI","(8/13)GACNNNNNNTTYG(11/6)");
		rebase.add("AscI","GG^CGCGCC");
		rebase.add("AsuII","TT^CGAA");
		rebase.add("AvaI","C^YCGRG");
		rebase.add("AvaII","G^GWCC");
		rebase.add("AvrII","C^CTAGG");
		rebase.add("BaeI","(10/15)ACNNNNGTAYC(12/7)");
		rebase.add("BalI","TGG^CCA");
		rebase.add("BamHI","G^GATCC");
		rebase.add("BarI","(7/12)GAAGNNNNNNTAC(12/7)");
		rebase.add("BbvCI","CCTCAGC(-5/-2)");
		rebase.add("BbvI","GCAGC(8/12)");
		rebase.add("BccI","CCATC(4/5)");
		rebase.add("BcgI","(10/12)CGANNNNNNTGC(12/10)");
		rebase.add("BciVI","GTATCC(6/5)");
		rebase.add("BclI","T^GATCA");
		rebase.add("BglI","GCCNNNN^NGGC");
		rebase.add("BglII","A^GATCT");
		rebase.add("BisI","GC^NGC");
		rebase.add("BplI","(8/13)GAGNNNNNCTC(13/8)");
		rebase.add("Bpu10I","CCTNAGC(-5/-2)");
		rebase.add("BsaAI","YAC^GTR");
		rebase.add("BsaBI","GATNN^NNATC");
		rebase.add("BsaXI","(9/12)ACNNNNNCTCC(10/7)");
		rebase.add("BseMII","CTCAG(10/8)");
		rebase.add("BsePI","G^CGCGC");
		rebase.add("BseRI","GAGGAG(10/8)");
		rebase.add("BseSI","GKGCM^C");
		rebase.add("BseYI","CCCAGC(-5/-1)");
		rebase.add("BsgI","GTGCAG(16/14)");
		rebase.add("BsmAI","GTCTC(1/5)");
		rebase.add("BsmI","GAATGC(1/-1)");
		rebase.add("Bsp1407I","T^GTACA");
		rebase.add("BspHI","T^CATGA");
		rebase.add("BspMI","ACCTGC(4/8)");
		rebase.add("BsrBI","CCGCTC(-3/-3)");
		rebase.add("BsrDI","GCAATG(2/0)");
		rebase.add("BsrI","ACTGG(1/-1)");
		rebase.add("BstEII","G^GTNACC");
		rebase.add("BstXI","CCANNNNN^NTGG");
		rebase.add("BtgZI","GCGATG(10/14)");
		rebase.add("BtrI","CACGTC(-3/-3)");
		rebase.add("BtsI","GCAGTG(2/0)");
		rebase.add("BtsIMutI","CAGTG(2/0)");
		rebase.add("Cac8I","GCN^NGC");
		rebase.add("Cfr10I","R^CCGGY");
		rebase.add("ClaI","AT^CGAT");
		rebase.add("CspCI","(11/13)CAANNNNNGTGG(12/10)");
		rebase.add("CviJI","RG^CY");
		rebase.add("DdeI","C^TNAG");
		rebase.add("DpnI","GA^TC");
		rebase.add("DraIII","CACNNN^GTG");
		rebase.add("DrdI","GACNNNN^NNGTC");
		rebase.add("Eam1105I","GACNNN^NNGTC");
		rebase.add("EciI","GGCGGA(11/9)");
		rebase.add("Eco31I","GGTCTC(1/5)");
		rebase.add("Eco47III","AGC^GCT");
		rebase.add("Eco57I","CTGAAG(16/14)");
		rebase.add("EcoNI","CCTNN^NNNAGG");
		rebase.add("EcoP15I","CAGCAG(25/27)");
		rebase.add("EcoRI","G^AATTC");
		rebase.add("EcoRII","^CCWGG");
		rebase.add("EcoRV","GAT^ATC");
		rebase.add("Esp3I","CGTCTC(1/5)");
		rebase.add("FaiI","YA^TR");
		rebase.add("FalI","(8/13)AAGNNNNNCTT(13/8)");
		rebase.add("FauI","CCCGC(4/6)");
		rebase.add("Fnu4HI","GC^NGC");
		rebase.add("FokI","GGATG(9/13)");
		rebase.add("FseI","GGCCGG^CC");
		rebase.add("FspAI","RTGC^GCAY");
		rebase.add("FspEI","CC(12/16)");
		rebase.add("GlaI","GC^GC");
		rebase.add("GsuI","CTGGAG(16/14)");
		rebase.add("HaeII","RGCGC^Y");
		rebase.add("HaeIII","GG^CC");
		rebase.add("HgaI","GACGC(5/10)");
		rebase.add("HhaI","GCG^C");
		rebase.add("HindII","GTY^RAC");
		rebase.add("HindIII","A^AGCTT");
		rebase.add("HinfI","G^ANTC");
		rebase.add("HpaI","GTT^AAC");
		rebase.add("HpaII","C^CGG");
		rebase.add("HphI","GGTGA(8/7)");
		rebase.add("Hpy188I","TCN^GA");
		rebase.add("Hpy99I","CGWCG^");
		rebase.add("KpnI","GGTAC^C");
		rebase.add("KroI","G^CCGGC");
		rebase.add("LpnPI","CCDG(10/14)");
		rebase.add("MaeI","C^TAG");
		rebase.add("MaeII","A^CGT");
		rebase.add("MaeIII","^GTNAC");
		rebase.add("MauBI","CG^CGCGCG");
		rebase.add("MboI","^GATC");
		rebase.add("MboII","GAAGA(8/7)");
		rebase.add("MfeI","C^AATTG");
		rebase.add("MluI","A^CGCGT");
		rebase.add("MmeI","TCCRAC(20/18)");
		rebase.add("MnlI","CCTC(7/6)");
		rebase.add("MseI","T^TAA");
		rebase.add("MslI","CAYNN^NNRTG");
		rebase.add("MspJI","CNNR(9/13)");
		rebase.add("MwoI","GCNNNNN^NNGC");
		rebase.add("NaeI","GCC^GGC");
		rebase.add("NarI","GG^CGCC");
		rebase.add("NcoI","C^CATGG");
		rebase.add("NdeI","CA^TATG");
		rebase.add("NheI","G^CTAGC");
		rebase.add("NlaIII","CATG^");
		rebase.add("NlaIV","GGN^NCC");
		rebase.add("NmeAIII","GCCGAG(21/19)");
		rebase.add("NotI","GC^GGCCGC");
		rebase.add("NruI","TCG^CGA");
		rebase.add("NspI","RCATG^Y");
		rebase.add("OliI","CACNN^NNGTG");
		rebase.add("PacI","TTAAT^TAA");
		rebase.add("PasI","CC^CWGGG");
		rebase.add("PcsI","WCGNNNN^NNNCGW");
		rebase.add("PflMI","CCANNNN^NTGG");
		rebase.add("PfoI","T^CCNGGA");
		rebase.add("PleI","GAGTC(4/5)");
		rebase.add("PmaCI","CAC^GTG");
		rebase.add("PmeI","GTTT^AAAC");
		rebase.add("PpuMI","RG^GWCCY");
		rebase.add("PshAI","GACNN^NNGTC");
		rebase.add("PsiI","TTA^TAA");
		rebase.add("PspXI","VC^TCGAGB");
		rebase.add("PsrI","(7/12)GAACNNNNNNTAC(12/7)");
		rebase.add("PstI","CTGCA^G");
		rebase.add("PvuI","CGAT^CG");
		rebase.add("PvuII","CAG^CTG");
		rebase.add("RsaI","GT^AC");
		rebase.add("RsrII","CG^GWCCG");
		rebase.add("SacI","GAGCT^C");
		rebase.add("SacII","CCGC^GG");
		rebase.add("SalI","G^TCGAC");
		rebase.add("SapI","GCTCTTC(1/4)");
		rebase.add("ScaI","AGT^ACT");
		rebase.add("ScrFI","CC^NGG");
		rebase.add("SduI","GDGCH^C");
		rebase.add("SetI","ASST^");
		rebase.add("SexAI","A^CCWGGT");
		rebase.add("SfaNI","GCATC(5/9)");
		rebase.add("SfiI","GGCCNNNN^NGGCC");
		rebase.add("SgeI","CNNGNNNNNNNNN^");
		rebase.add("SgfI","GCGAT^CGC");
		rebase.add("SgrAI","CR^CCGGYG");
		rebase.add("SgrDI","CG^TCGACG");
		rebase.add("SmaI","CCC^GGG");
		rebase.add("SmlI","C^TYRAG");
		rebase.add("SnaBI","TAC^GTA");
		rebase.add("SpeI","A^CTAGT");
		rebase.add("SphI","GCATG^C");
		rebase.add("Sse8387I","CCTGCA^GG");
		rebase.add("SspI","AAT^ATT");
		rebase.add("StuI","AGG^CCT");
		rebase.add("StyI","C^CWWGG");
		rebase.add("SwaI","ATTT^AAAT");
		rebase.add("TaqI","T^CGA");
		rebase.add("TaqII","GACCGA(11/9)");
		rebase.add("TatI","W^GTACW");
		rebase.add("TauI","GCSG^C");
		rebase.add("TfiI","G^AWTC");
		rebase.add("TseI","G^CWGC");
		rebase.add("Tsp45I","^GTSAC");
		rebase.add("TspDTI","ATGAA(11/9)");
		rebase.add("TspGWI","ACGGA(11/9)");
		rebase.add("TspRI","CASTGNN^");
		rebase.add("Tth111I","GACN^NNGTC");
		rebase.add("VspI","AT^TAAT");
		rebase.add("XbaI","T^CTAGA");
		rebase.add("XcmI","CCANNNNN^NNNNTGG");
		rebase.add("XhoI","C^TCGAG");
		rebase.add("XhoII","R^GATCY");
		rebase.add("XmnI","GAANN^NNTTC");
		return rebase;
		}
	}
