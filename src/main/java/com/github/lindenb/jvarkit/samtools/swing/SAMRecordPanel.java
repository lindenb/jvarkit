package com.github.lindenb.jvarkit.samtools.swing;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.beans.PropertyChangeListener;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

@SuppressWarnings("serial")
public class SAMRecordPanel extends JPanel {
private SAMRecord samRecord = null;
private final ReadTableModel readTableModel;
private final CigarTableModel cigarTableModel;
private final MetadataTableModel metaTableModel;
private final SAMFlagTableModel flagTableModel;
private final SeqTableModel seqTableModel;
private final OtherAlignTableModel otherAlignTableModel;

static JPanel wrap(String value,final TableModel tm) {
	final JTable t = new JTable(tm);
	t.setAutoResizeMode(JTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
	JScrollPane scr = new JScrollPane(t);
	JPanel jp=new JPanel(new BorderLayout(4, 4));
	jp.setBorder(BorderFactory.createTitledBorder(value));
	jp.add(scr,BorderLayout.CENTER);
	return jp;	
	}

public SAMRecordPanel() {
	super(new GridLayout(2,3,5,5));
	this.add(wrap("Read",this.readTableModel= new ReadTableModel()));
	this.add(wrap("Cigar",this.cigarTableModel = new CigarTableModel()));
	this.add(wrap("Meta",this.metaTableModel = new MetadataTableModel()));
	this.add(wrap("Flags",this.flagTableModel = new SAMFlagTableModel()));
	this.add(wrap("Seq",this.seqTableModel = new SeqTableModel()));
	this.add(wrap("XA:",this.otherAlignTableModel = new OtherAlignTableModel()));
	}


public SAMRecordPanel(final SAMRecord rec) {
	this();
	setSamRecord(rec);
	}

public void setSamRecord(final SAMRecord samRecord) {
	SAMRecord old = this.samRecord;
	this.samRecord = samRecord;
	if(old!=this.samRecord) {
		this.readTableModel.fireTableDataChanged();
		this.cigarTableModel.fireTableDataChanged();
		this.metaTableModel.fireTableDataChanged();
		this.flagTableModel.fireTableDataChanged();
		this.seqTableModel.fireTableDataChanged();
		this.otherAlignTableModel.fireTableDataChanged();
		}
	}

public SAMRecord getSamRecord() {
	return samRecord;
	}


private class ReadTableModel  extends AbstractTableModel
	{
	@Override
	public String getColumnName(final int column) {
		switch(column) {
			case 0: return "Key";
			case 1: return "Value";
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return String.class;
		}
	
	private String flags(int f) {
		return SAMFlag.getFlags(f).stream().map(T->T.name()).collect(Collectors.joining(","));
		}
	private String trim(final String s) {
		if(s==null) return "*";
		if(s.length()>500) return s.substring(0, 500)+"...";
		return s;
	}
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final SAMRecord rec = getSamRecord();
		switch(rowIndex) {
			case 0: return columnIndex==0?"QNAME":(rec==null?null:rec.getReadName());
			case 1: return columnIndex==0?"FLAG":(rec==null?null:String.valueOf(rec.getFlags())+" "+flags(rec.getFlags()));
			case 2: return columnIndex==0?"RNAME":(rec==null?null:rec.getReferenceName());
			case 3: return columnIndex==0?"UPOS":(rec==null || rec.getReadUnmappedFlag()?null:StringUtils.niceInt(rec.getUnclippedStart()));
			case 4: return columnIndex==0?"POS":(rec==null || rec.getReadUnmappedFlag()?null:StringUtils.niceInt(rec.getAlignmentStart()));
			case 5: return columnIndex==0?"END":(rec==null || rec.getReadUnmappedFlag()?null:StringUtils.niceInt(rec.getAlignmentEnd()));
			case 6: return columnIndex==0?"UEND":(rec==null || rec.getReadUnmappedFlag()?null:StringUtils.niceInt(rec.getUnclippedEnd()));
			case 7: return columnIndex==0?"RLEN":(rec==null || rec.getReadUnmappedFlag()?null:StringUtils.niceInt(rec.getLengthOnReference()));
			case 8: return columnIndex==0?"MAPQ":(rec==null || rec.getReadUnmappedFlag()?null:String.valueOf(rec.getMappingQuality()));
			case 9: return columnIndex==0?"CIGAR":(rec==null || rec.getReadUnmappedFlag()?null:rec.getCigarString());
			case 10: return columnIndex==0?"ISIZE":(rec==null || rec.getReadUnmappedFlag()?null:StringUtils.niceInt(rec.getInferredInsertSize()));
			case 11: return columnIndex==0?"MNAME":(rec==null || !rec.getReadPairedFlag() || !rec.getMateUnmappedFlag()?null:rec.getMateReferenceName());
			case 12: return columnIndex==0?"MPOS":(rec==null || !rec.getReadPairedFlag() || !rec.getMateUnmappedFlag()?null:StringUtils.niceInt(rec.getMateAlignmentStart()));
			case 13: return columnIndex==0?"SEQ":(rec==null?null:trim(rec.getReadString()));
			case 14: return columnIndex==0?"QUAL":(rec==null?null:trim(rec.getBaseQualityString()));
			default: return null;
			}
		}
	@Override
	public final int getRowCount() {
		return 15;
		}
	@Override
	public final int getColumnCount() {
		return 2;
		}
	}



private class CigarTableModel  extends AbstractTableModel
	{
	@Override
	public String getColumnName(final int column) {
		switch(column) {
			case 0: return "Op";
			case 1: return "Len";
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	@Override
	public Class<?> getColumnClass(int column) {
			switch(column) {
			case 0: return String.class;
			case 1: return Integer.class;
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final SAMRecord rec = getSamRecord();
		if(rec==null) return null;
		final Cigar cig = rec.getCigar();
		if(cig==null) return null;
		final CigarElement op = cig.getCigarElement(rowIndex);
		switch(columnIndex) {
			case 0: return op.getOperator().name();
			case 1: return StringUtils.niceInt(op.getLength());
			default: return null;
			}
		}
	@Override
	public final int getRowCount() {
		final SAMRecord rec = getSamRecord();
		if(rec==null || rec.getReadUnmappedFlag() || rec.getCigar()==null) return 0;
		return rec.getCigar().numCigarElements();
		}
	@Override
	public final int getColumnCount() {
		return 2;
		}
	}

private class MetadataTableModel  extends AbstractTableModel
	{
	@Override
	public String getColumnName(final int column) {
		switch(column) {
			case 0: return "Name";
			case 1: return "Type";
			case 2: return "Value";
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	@Override
	public Class<?> getColumnClass(int column) {
			switch(column) {
			case 0: return String.class;
			case 1: return String.class;
			case 2: return Object.class;
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final SAMRecord rec = getSamRecord();
		if(rec==null) return null;
		final List<SAMRecord.SAMTagAndValue> L = rec.getAttributes();
		if(L.isEmpty()) return null;
		final SAMRecord.SAMTagAndValue kv = L.get(rowIndex);
		switch(columnIndex) {
			case 0: return kv.tag;
			case 1: return kv.value.getClass().getSimpleName();
			case 2: return kv.value;
			default: return null;
			}
		}
	@Override
	public final int getRowCount() {
		final SAMRecord rec = getSamRecord();
		if(rec==null) return 0;
		return rec.getAttributes().size();
		}
	@Override
	public final int getColumnCount() {
		return 3;
		}
	}


private class SAMFlagTableModel  extends AbstractTableModel
	{
	@Override
	public String getColumnName(final int column) {
		switch(column) {
			case 0: return "Flags";
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	@Override
	public Class<?> getColumnClass(int column) {
			switch(column) {
			case 0: return String.class;
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final SAMRecord rec = getSamRecord();
		if(rec==null) return null;
		for(Iterator<SAMFlag> iter= SAMFlag.getFlags(rec.getFlags()).iterator();iter.hasNext();) {
			SAMFlag flg = iter.next();
			if(rowIndex==0) return flg.name();
			rowIndex--;
			}
		return null;
		}
	@Override
	public final int getRowCount() {
		final SAMRecord rec = getSamRecord();
		if(rec==null) return 0;
		return SAMFlag.getFlags(rec.getFlags()).size();
		}
	@Override
	public final int getColumnCount() {
		return 1;
		}
	}




private class SeqTableModel  extends AbstractTableModel
	{
	@Override
	public String getColumnName(final int column) {
		switch(column) {
			case 0: return "Op";
			case 1: return "REF1";
			case 2: return "READ1";
			case 3: return "BASE";
			case 4: return "QUAL";
			default: throw new IndexOutOfBoundsException(column);
			}
		}
	@Override
	public Class<?> getColumnClass(int column) {
		return String.class;
		}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final SAMRecord rec = getSamRecord();
		final Cigar cigar = rec.getCigar();
		if(cigar==null) return 0;
		final byte[] bases = rec.getReadBases();
		byte[] quals = rec.getBaseQualities();
		int readPos0 = 0;
		int refPos1 = rec.getUnclippedStart();
		final Function<Integer,String> getBase  = N->{
			if(bases==null || bases==SAMRecord.NULL_SEQUENCE) return "*";
			if(N<0 || N>= bases.length) return "?";
			return String.valueOf((char)bases[N]);
			};
		final Function<Integer,String> getQual  = N->{
			if(quals==null ||quals==SAMRecord.NULL_QUALS) return "*";
			if(N<0 || N>= quals.length) return "?";
			return String.valueOf(SAMUtils.phredToFastq(quals[N]));
			};
		for(CigarElement ce:cigar) {
			final CigarOperator op =ce.getOperator();
			for(int i=0;i< ce.getLength();i++) {
				int next_read = readPos0;
				int next_ref = refPos1;
				String base = null;
				String qual = null;
				String basePosStr = null;
				String refPosStr = null;
				switch(op) {
					case P: break;
					case H: {
						base = "*";
						qual = "*";
						refPosStr = StringUtils.niceInt(refPos1);
						next_ref  = refPos1 +1 ;
						break;
						}
					case X: case M: case S: case EQ:{
						base = getBase.apply(readPos0);
						if(op==CigarOperator.X || op==CigarOperator.S) base=base.toLowerCase();
						qual = getQual.apply(readPos0);
						basePosStr = StringUtils.niceInt(readPos0+1);
						refPosStr = StringUtils.niceInt(refPos1);
						next_read= readPos0+1;
						next_ref  = refPos1 +1 ;
						break;
						}
					case D: case N: {
						base = "-";
						qual = null;
						refPosStr = StringUtils.niceInt(refPos1);
						next_ref  = refPos1 +1 ;
						break;
						}
					case I:  {
						base = getBase.apply(readPos0);
						qual = getQual.apply(readPos0);
						basePosStr = StringUtils.niceInt(readPos0+1);
						next_read= readPos0+1;
						break;
						}
					default: throw new IllegalArgumentException();
					}
				if(rowIndex==0) {
					switch(columnIndex) {
						case 0: return op.name();
						case 1: return refPosStr;
						case 2: return basePosStr;
						case 3: return base;
						case 4: return qual;
						default: return null;
						}
					}
				rowIndex--;
				readPos0 = next_read;
				refPos1=  next_ref;
				}
			}
		return null;
		}
	@Override
	public final int getRowCount() {
		final SAMRecord rec = getSamRecord();
		if(rec==null || rec.getReadUnmappedFlag()) return 0;
		final Cigar cigar = rec.getCigar();
		if(cigar==null) return 0;
		return cigar.getCigarElements().stream().mapToInt(CE->CE.getLength()).sum();
		}
	@Override
	public final int getColumnCount() {
		return 5;
		}
	}


private class OtherAlignTableModel  extends AbstractTableModel
{
@Override
public String getColumnName(final int column) {
	switch(column) {
		case 0: return "Contig";
		case 1: return "Pos";
		case 2: return "Cigar";
		default: throw new IndexOutOfBoundsException(column);
		}
	}
@Override
public Class<?> getColumnClass(int column) {
		switch(column) {
		case 0: return String.class;
		default: throw new IndexOutOfBoundsException(column);
		}
	}

@Override
public Object getValueAt(int rowIndex, int columnIndex) {
	final SAMRecord rec = getSamRecord();
	if(rec==null) return null;
	final List<SAMRecord> L= SAMUtils.getOtherCanonicalAlignments(rec);
	final SAMRecord rec2  = L.get(rowIndex);
	switch(columnIndex) {
		case 0: return rec2.getContig();
		case 1: return StringUtils.niceInt(rec2.getAlignmentStart());
		case 2: return rec2.getCigarString();
		}
	return null;
	}
@Override
public final int getRowCount() {
	final SAMRecord rec = getSamRecord();
	if(rec==null) return 0;
	return SAMUtils.getOtherCanonicalAlignments(rec).size();
	}
@Override
public final int getColumnCount() {
	return 3;
	}
}


}
