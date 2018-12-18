package com.github.lindenb.jvarkit.table;

public class AbstractTableExporter {
private boolean printEmptyTable = false;
private boolean removeEmptyColumns = true;

public void setPrintEmptyTable(boolean printEmptyTable) {
	this.printEmptyTable = printEmptyTable;
}
public boolean isPrintEmptyTable() {
	return printEmptyTable;
	}
public void setRemoveEmptyColumns(boolean removeEmptyColumns) {
	this.removeEmptyColumns = removeEmptyColumns;
}
public boolean isRemoveEmptyColumns() {
	return removeEmptyColumns;
	}


}
