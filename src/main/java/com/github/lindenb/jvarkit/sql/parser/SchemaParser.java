/* SchemaParser.java */
/* Generated By:JavaCC: Do not edit this line. SchemaParser.java */
package com.github.lindenb.jvarkit.sql.parser;
import htsjdk.samtools.util.Locatable;
import com.github.lindenb.jvarkit.io.IOUtils;
import java.util.*;
import java.io.*;


/**

java -cp ~/packages/javacc/target/javacc.jar javacc -OUTPUT_DIRECTORY=src/main/java/com/github/lindenb/jvarkit/sql/parser -JDK_VERSION=1.8 src/main/resources/javacc/com/github/lindenb/jvarkit/sql/parser/SchemaParser.jj 


**/


public class SchemaParser implements SchemaParserConstants {
        private final List<SchemaParser.Table> tables=new ArrayList<SchemaParser.Table>();

        /** map a line to an element of the table, implements Locatable if UCSC columns for chrom/start/end were found */
        public interface Row extends Map<String,String>, Locatable {
                /** return table */
                public Table getTable();
                /** return i-th token */
                public String get(int idx);
                }

        /** describe a SQL table */
        public interface Table {
                /** return table name */
                public String getName();
                /** return columns for that table */
                public List<Column> getColumns();
                /** get column by it's name */
                public Optional<Column> getColumnByName(final String colName);
                /** convert line to Row */
                public Row toMap(final String[] tokens);
                /** convert line to Row */
                public Row toMap(final List<String> tokens);
                /**  implements Locatable if UCSC columns for chrom/start/end were found */
                public boolean isLocatable();
                }

        /** describe a column of a SQL table */
        public interface Column {
                /** get associated table */
                public Table getTable();
                /** get O-based index of the column */
                public int getIndex();
                /** return column name */
                public String getName();
                }

        private static class RowImpl extends AbstractMap<String,String> implements Row
                {
                /** cached value for Locatable::getStart */
                private int loc_start =-1;
                /** cached value for Locatable::getEnd */
                private int loc_end =-1;
                /** owner table */
                private final TableImpl table;
                /** original line */
                private final List<String> tokens;
                RowImpl( final TableImpl table, final List<String> tokens) {
                        this.table = table;
                        this.tokens = tokens;
                        }
                @Override
                public String getContig() {
                        if(table.contig_col==-1) throw new IllegalStateException("table "+table.getName()+" doesn't support locatable");
                        return get(table.contig_col);
                        }

        @Override
        public int getStart() {
                        if(loc_start!=-1) return loc_start;
                if(table.start_col==-1) throw new IllegalStateException("table "+table.getName()+" doesn't support locatable");
                loc_start = Integer.parseInt(get(table.start_col))+1; // UCSC IS 0-BASED
                return loc_start;
                }
        @Override
        public int getEnd() {
                        if(loc_end!=-1) return loc_end;
                if(table.end_col==-1) throw new IllegalStateException("table "+table.getName()+" doesn't support locatable");
                loc_end =  Integer.parseInt(get(table.end_col));
                return loc_end;
                }

                @Override
                public Table getTable() {
                        return this.table;
                        }

        @Override
        public int size() {
                return table.getColumns().size();
                }

        @Override
        public Set<Entry<String, String>> entrySet() {
                final Set<Entry<String, String>> set = new HashSet<>(size());
                for(int i=0;i< size();i++) {
                        set.add(new AbstractMap.SimpleEntry<String,String>(
                                        this.table.columns.get(i).getName(),get(i)));
                        }
                return set;
                }

        @Override
        public boolean containsKey(Object key) {
                if(!(key instanceof String)) return false;
                return getTable().getColumnByName(String.class.cast(key)).isPresent();
                }

        @Override
        public String getOrDefault(Object key, String defaultValue) {
                if(!(key instanceof String)) return defaultValue;
                final Optional<Column> col= getTable().getColumnByName(String.class.cast(key));
                if(!col.isPresent()) return defaultValue;
                return get(col.get().getIndex());
                }

        @Override
        public String get(Object key) {
                return getOrDefault(key,null);
                }

        @Override
        public String get(int idx) {
                        return idx>=this.tokens.size()?"":this.tokens.get(idx);
                }

                @Override
                public int hashCode() {
                        return tokens.hashCode();
                        }
                @Override
                public String toString() {
                        return String.join("\t", tokens);
                        }
                }

        private class ColumnImpl implements Column
                {
                private int index;
                private String colName;
                private String sqlType;
                private TableImpl table;
                ColumnImpl() {

                        }

                @Override
                public int getIndex() {
                        return this.index;
                        }

                @Override
                public String getName() {
                        return this.colName;
                        }

                @Override
                public int hashCode() {
                        return colName.hashCode();
                        }

                @Override
                public Table getTable() {
                        return this.table;
                        }
                @Override
                public String toString() {
                        return getName()+" "+sqlType;
                        }
                }
        private class TableImpl implements Table
                {
                private String tableName;
                private int contig_col = -1;
                private int start_col = -1;
                private int end_col = -1;

                final List<Column> columns=new ArrayList<Column>();
                final Map<String,Column> name2column = new HashMap<String,Column>();
                void findLocatableColumns() {
                        Column c = name2column.get("chrom");
                        if(c!=null) contig_col = c.getIndex();
                        c = name2column.get("txStart");
                        if(c!=null)  start_col = c.getIndex();
                        if(start_col<0) {
                                c = name2column.get("chromStart");
                                if(c!=null)  start_col = c.getIndex();
                                }
                        if(start_col<0) {
                                c = name2column.get("cdsStart");
                                if(c!=null)  start_col = c.getIndex();
                                }
                        c = name2column.get("txEnd");
                        if(c!=null)  end_col = c.getIndex();
                        if(end_col<0) {
                                c = name2column.get("cdsEnd");
                                if(c!=null)  end_col = c.getIndex();
                                }
                        if(end_col<0) {
                                c = name2column.get("chromEnd");
                                if(c!=null)  end_col = c.getIndex();
                                }
                        }

                void add(final ColumnImpl col) {
                        if(col==null) return;
                        col.index = this.columns.size();
                        col.table = this;
                        name2column.put(col.getName(),col);
                        this.columns.add(col);
                        }

                /**  implements Locatable if UCSC columns for chrom/start/end were found */
                @Override
                public boolean isLocatable() {
                        return contig_col>=0 && start_col>=0 && end_col>=0;
                        }


                @Override
                public int hashCode() {
                        return tableName.hashCode();
                        }

                @Override
                public List<Column> getColumns() {
                        return this.columns;
                        }
                @Override
                public Optional<Column> getColumnByName(final String colName) {
                        return Optional.ofNullable(name2column.get(colName));
                        }

                @Override
                public Row toMap(final String[] tokens) {
                        return toMap(Arrays.asList(tokens));
                        }
                @Override
                public Row toMap(final List<String> tokens) {
                        return new RowImpl(this,tokens);
                        }

                @Override
                public String getName() {
                        return this.tableName;
                        }
                @Override
                public String toString() {
                        return "CREATE TABLE "+getName()+"(" + columns +")";
                        }
                }


    public static Table parseTable(final String uri) throws IOException {
        try(InputStream in = IOUtils.openURIForReading(uri)) {
                return parseTable(in);
                }
        }

    public static List<Table> parseTables(final String uri) throws IOException {
        try(InputStream in = IOUtils.openURIForReading(uri)) {
                return parseTables(in);
                }
        }


         public static Table parseTable(final InputStream in) throws IOException {
                final List<Table> tables = parseTables(in);
                if(tables.size()!=1) throw new IllegalStateException("expected only one table in schema but got "+tables.size());
                return tables.get(0);
                }

        public static List<Table> parseTables(InputStream in) throws IOException
                        {
                        try {
                                final SchemaParser parser = new SchemaParser(in);
                                parser.input();
                                return parser.tables;
                                }
                        catch(final Throwable err) {
                                throw new IOException(err);
                                }
                        }

  final private void input() throws ParseException {TableImpl t;
    label_1:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case DROP:
      case CREATE:
      case USE:
      case SET:
      case SEMICOLON:{
        ;
        break;
        }
      default:
        jj_la1[0] = jj_gen;
        break label_1;
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case DROP:{
        drop();
        break;
        }
      case SET:{
        set();
        break;
        }
      case CREATE:{
        create();
        break;
        }
      case USE:{
        use();
        break;
        }
      case SEMICOLON:{
        jj_consume_token(SEMICOLON);
        break;
        }
      default:
        jj_la1[1] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
    }
    jj_consume_token(0);
}

  final private void use() throws ParseException {
    jj_consume_token(USE);
    identifier();
    jj_consume_token(SEMICOLON);
}

  final private void create() throws ParseException {
    jj_consume_token(CREATE);
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case DATABASE:{
      createDatabase();
      break;
      }
    case TABLE:
    case TEMPORARY:{
      createTable();
      break;
      }
    default:
      jj_la1[2] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    jj_consume_token(SEMICOLON);
}

  final private void drop() throws ParseException {
    jj_consume_token(DROP);
    jj_consume_token(TABLE);
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IF:{
      jj_consume_token(IF);
      jj_consume_token(EXISTS);
      break;
      }
    default:
      jj_la1[3] = jj_gen;
      ;
    }
    identifier();
    jj_consume_token(SEMICOLON);
}

  final private void set() throws ParseException {
    jj_consume_token(SET);
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IDENTIFIER1:
    case IDENTIFIER2:{
      identifier();
      break;
      }
    case VARIABLE:{
      jj_consume_token(VARIABLE);
      break;
      }
    default:
      jj_la1[4] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    jj_consume_token(EQ);
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IDENTIFIER1:
    case IDENTIFIER2:{
      identifier();
      break;
      }
    case VARIABLE:{
      jj_consume_token(VARIABLE);
      break;
      }
    default:
      jj_la1[5] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    jj_consume_token(SEMICOLON);
}

  final private void createDatabase() throws ParseException {String s;
    jj_consume_token(DATABASE);
    s = identifier();
System.err.println("#############create database "+s);
}

  final private void createTable() throws ParseException {TableImpl table=new TableImpl(); String tableName; ColumnImpl col=null;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case TEMPORARY:{
      jj_consume_token(TEMPORARY);
      break;
      }
    default:
      jj_la1[6] = jj_gen;
      ;
    }
    jj_consume_token(TABLE);
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IF:{
      jj_consume_token(IF);
      jj_consume_token(NOT);
      jj_consume_token(EXISTS);
      break;
      }
    default:
      jj_la1[7] = jj_gen;
      ;
    }
    tableName = identifier();
    jj_consume_token(LPAR);
    col = component();
table.add(col);
    label_2:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case COMMA:{
        ;
        break;
        }
      default:
        jj_la1[8] = jj_gen;
        break label_2;
      }
      jj_consume_token(COMMA);
      col = component();
table.add(col);
    }
    jj_consume_token(RPAR);
    label_3:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case DEFAULTVALUE:
      case COLLATE:
      case CHARSET:
      case IDENTIFIER1:
      case IDENTIFIER2:
      case EQ:{
        ;
        break;
        }
      default:
        jj_la1[9] = jj_gen;
        break label_3;
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case IDENTIFIER1:
      case IDENTIFIER2:{
        identifier();
        break;
        }
      case EQ:{
        jj_consume_token(EQ);
        break;
        }
      case CHARSET:{
        jj_consume_token(CHARSET);
        break;
        }
      case COLLATE:{
        jj_consume_token(COLLATE);
        break;
        }
      case DEFAULTVALUE:{
        jj_consume_token(DEFAULTVALUE);
        break;
        }
      default:
        jj_la1[10] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
    }
table.tableName = tableName;
                table.findLocatableColumns();
                tables.add(table);
}

  final private ColumnImpl component() throws ParseException {ColumnImpl c=null;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case KEY:
    case UNIQUE:
    case PRIMARY:{
      key();
      break;
      }
    case IDENTIFIER1:
    case IDENTIFIER2:{
      c = column();
{if ("" != null) return c;}
      break;
      }
    default:
      jj_la1[11] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
{if ("" != null) return c;}
    throw new Error("Missing return statement in function");
}

  final private void key() throws ParseException {
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case UNIQUE:
    case PRIMARY:{
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case UNIQUE:{
        jj_consume_token(UNIQUE);
        break;
        }
      case PRIMARY:{
        jj_consume_token(PRIMARY);
        break;
        }
      default:
        jj_la1[12] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
      break;
      }
    default:
      jj_la1[13] = jj_gen;
      ;
    }
    jj_consume_token(KEY);
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IDENTIFIER1:
    case IDENTIFIER2:{
      identifier();
      break;
      }
    default:
      jj_la1[14] = jj_gen;
      ;
    }
    jj_consume_token(LPAR);
    identifier();
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case LPAR:{
      jj_consume_token(LPAR);
      integer();
      jj_consume_token(RPAR);
      break;
      }
    default:
      jj_la1[15] = jj_gen;
      ;
    }
    label_4:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case COMMA:{
        ;
        break;
        }
      default:
        jj_la1[16] = jj_gen;
        break label_4;
      }
      jj_consume_token(COMMA);
      identifier();
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case LPAR:{
        jj_consume_token(LPAR);
        integer();
        jj_consume_token(RPAR);
        break;
        }
      default:
        jj_la1[17] = jj_gen;
        ;
      }
    }
    jj_consume_token(RPAR);
}

  final private ColumnImpl column() throws ParseException {ColumnImpl c=new ColumnImpl();String s; String javaType=null;boolean nil=true;
    s = identifier();
    javaType = colType();
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case NOT:
    case NULL:{
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case NOT:{
        jj_consume_token(NOT);
nil=false;
        break;
        }
      default:
        jj_la1[18] = jj_gen;
        ;
      }
      jj_consume_token(NULL);
      break;
      }
    default:
      jj_la1[19] = jj_gen;
      ;
    }
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case DEFAULTVALUE:{
      jj_consume_token(DEFAULTVALUE);
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case NULL:{
        jj_consume_token(NULL);
        break;
        }
      case SIMPLE_QUOTE_LITERAL:
      case DOUBLE_QUOTE_LITERAL:{
        quoted();
        break;
        }
      case INT:{
        integer();
        break;
        }
      default:
        jj_la1[20] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
      break;
      }
    default:
      jj_la1[21] = jj_gen;
      ;
    }
c.colName=s;
          c.sqlType=javaType;
          {if ("" != null) return c;}
    throw new Error("Missing return statement in function");
}

  final private String colType() throws ParseException {String s; Token t;Set<String> set;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case VARCHAR:
    case CHAR:{
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case VARCHAR:{
        t = jj_consume_token(VARCHAR);
        break;
        }
      case CHAR:{
        t = jj_consume_token(CHAR);
        break;
        }
      default:
        jj_la1[22] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case LPAR:{
        dimension();
        break;
        }
      default:
        jj_la1[23] = jj_gen;
        ;
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case CHARACTER:
      case CHARSET:{
        charset();
        break;
        }
      default:
        jj_la1[24] = jj_gen;
        ;
      }
{if ("" != null) return t.image;}
      break;
      }
    case LONGBLOB:
    case BLOB:
    case MEDIUMBLOG:
    case TEXT:{
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case LONGBLOB:{
        t = jj_consume_token(LONGBLOB);
        break;
        }
      case MEDIUMBLOG:{
        t = jj_consume_token(MEDIUMBLOG);
        break;
        }
      case BLOB:{
        t = jj_consume_token(BLOB);
        break;
        }
      case TEXT:{
        t = jj_consume_token(TEXT);
        break;
        }
      default:
        jj_la1[25] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
{if ("" != null) return t.image;}
      break;
      }
    case INTEGER:{
      t = jj_consume_token(INTEGER);
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case LPAR:{
        dimension();
        break;
        }
      default:
        jj_la1[26] = jj_gen;
        ;
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case UNSIGNED:
      case SIGNED:{
        switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
        case SIGNED:{
          jj_consume_token(SIGNED);
          break;
          }
        case UNSIGNED:{
          jj_consume_token(UNSIGNED);
          break;
          }
        default:
          jj_la1[27] = jj_gen;
          jj_consume_token(-1);
          throw new ParseException();
        }
        break;
        }
      default:
        jj_la1[28] = jj_gen;
        ;
      }
{if ("" != null) return t.image;}
      break;
      }
    case SMALLINT:
    case TINYINT:{
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case SMALLINT:{
        t = jj_consume_token(SMALLINT);
        break;
        }
      case TINYINT:{
        t = jj_consume_token(TINYINT);
        break;
        }
      default:
        jj_la1[29] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case LPAR:{
        dimension();
        break;
        }
      default:
        jj_la1[30] = jj_gen;
        ;
      }
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case UNSIGNED:
      case SIGNED:{
        switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
        case SIGNED:{
          jj_consume_token(SIGNED);
          break;
          }
        case UNSIGNED:{
          jj_consume_token(UNSIGNED);
          break;
          }
        default:
          jj_la1[31] = jj_gen;
          jj_consume_token(-1);
          throw new ParseException();
        }
        break;
        }
      default:
        jj_la1[32] = jj_gen;
        ;
      }
{if ("" != null) return t.image;}
      break;
      }
    case DATETIME:{
      t = jj_consume_token(DATETIME);
{if ("" != null) return t.image;}
      break;
      }
    case ENUM:{
      jj_consume_token(ENUM);
      set = stringset();
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case CHARACTER:
      case CHARSET:{
        charset();
        break;
        }
      default:
        jj_la1[33] = jj_gen;
        ;
      }
{if ("" != null) return "ENUM("+String.join(",",set)+")";}
      break;
      }
    case SET:{
      jj_consume_token(SET);
      set = stringset();
{if ("" != null) return "ENUM("+String.join(",",set)+")";}
      break;
      }
    case FLOAT:{
      t = jj_consume_token(FLOAT);
{if ("" != null) return t.image;}
      break;
      }
    case DOUBLE:{
      t = jj_consume_token(DOUBLE);
{if ("" != null) return t.image;}
      break;
      }
    default:
      jj_la1[34] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    throw new Error("Missing return statement in function");
}

  final private Set<String> stringset() throws ParseException {Set<String> set=new HashSet<String>();String s;
    jj_consume_token(LPAR);
    s = quoted();
set.add(s);
    label_5:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case COMMA:{
        ;
        break;
        }
      default:
        jj_la1[35] = jj_gen;
        break label_5;
      }
      jj_consume_token(COMMA);
      s = quoted();
set.add(s);
    }
    jj_consume_token(RPAR);
{if ("" != null) return set;}
    throw new Error("Missing return statement in function");
}

  final private void charset() throws ParseException {
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case CHARACTER:{
      jj_consume_token(CHARACTER);
      jj_consume_token(SET);
      break;
      }
    case CHARSET:{
      jj_consume_token(CHARSET);
      break;
      }
    default:
      jj_la1[36] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IDENTIFIER1:
    case IDENTIFIER2:{
      identifier();
      break;
      }
    case SIMPLE_QUOTE_LITERAL:
    case DOUBLE_QUOTE_LITERAL:{
      quoted();
      break;
      }
    default:
      jj_la1[37] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case COLLATE:{
      jj_consume_token(COLLATE);
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case IDENTIFIER1:
      case IDENTIFIER2:{
        identifier();
        break;
        }
      case SIMPLE_QUOTE_LITERAL:
      case DOUBLE_QUOTE_LITERAL:{
        quoted();
        break;
        }
      default:
        jj_la1[38] = jj_gen;
        jj_consume_token(-1);
        throw new ParseException();
      }
      break;
      }
    default:
      jj_la1[39] = jj_gen;
      ;
    }
}

  final private void dimension() throws ParseException {
    jj_consume_token(LPAR);
    integer();
    jj_consume_token(RPAR);
}

  final private String identifier() throws ParseException {Token t;String s;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IDENTIFIER1:{
      t = jj_consume_token(IDENTIFIER1);
s=t.image;
      break;
      }
    case IDENTIFIER2:{
      t = jj_consume_token(IDENTIFIER2);
s=t.image.substring(1,t.image.length()-1);
      break;
      }
    default:
      jj_la1[40] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
{if ("" != null) return s;}
    throw new Error("Missing return statement in function");
}

  final private int integer() throws ParseException {Token t;
    t = jj_consume_token(INT);
{if ("" != null) return Integer.parseInt(t.image);}
    throw new Error("Missing return statement in function");
}

  final private String quoted() throws ParseException {Token t;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case SIMPLE_QUOTE_LITERAL:{
      t = jj_consume_token(SIMPLE_QUOTE_LITERAL);
{if ("" != null) return t.image;}
      break;
      }
    case DOUBLE_QUOTE_LITERAL:{
      t = jj_consume_token(DOUBLE_QUOTE_LITERAL);
{if ("" != null) return t.image;}
      break;
      }
    default:
      jj_la1[41] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    throw new Error("Missing return statement in function");
}

  /** Generated Token Manager. */
  public SchemaParserTokenManager token_source;
  SimpleCharStream jj_input_stream;
  /** Current token. */
  public Token token;
  /** Next token. */
  public Token jj_nt;
  private int jj_ntk;
  private int jj_gen;
  final private int[] jj_la1 = new int[42];
  static private int[] jj_la1_0;
  static private int[] jj_la1_1;
  static {
	   jj_la1_init_0();
	   jj_la1_init_1();
	}
	private static void jj_la1_init_0() {
	   jj_la1_0 = new int[] {0x80e0,0x80e0,0x700,0x800,0x0,0x0,0x200,0x800,0x0,0x100000,0x100000,0x20000,0x0,0x0,0x0,0x0,0x0,0x0,0x1000,0x3000,0x2000,0x100000,0xc0000,0x0,0x0,0x50600000,0x0,0x80000000,0x80000000,0x21000000,0x0,0x80000000,0x80000000,0x0,0x7fed8000,0x0,0x0,0x0,0x0,0x0,0x0,0x0,};
	}
	private static void jj_la1_init_1() {
	   jj_la1_1 = new int[] {0x2000,0x2000,0x0,0x0,0xe00,0xe00,0x0,0x0,0x10000,0x1630,0x1630,0x606,0x6,0x6,0x600,0x4000,0x10000,0x4000,0x0,0x0,0xc0100,0x0,0x0,0x4000,0x28,0x0,0x4000,0x1,0x1,0x0,0x4000,0x1,0x1,0x28,0x0,0x10000,0x28,0xc0600,0xc0600,0x10,0x600,0xc0000,};
	}

  /** Constructor with InputStream. */
  public SchemaParser(java.io.InputStream stream) {
	  this(stream, null);
  }
  /** Constructor with InputStream and supplied encoding */
  public SchemaParser(java.io.InputStream stream, String encoding) {
	 try { jj_input_stream = new SimpleCharStream(stream, encoding, 1, 1); } catch(java.io.UnsupportedEncodingException e) { throw new RuntimeException(e); }
	 token_source = new SchemaParserTokenManager(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 42; i++) jj_la1[i] = -1;
  }

  /** Reinitialise. */
  public void ReInit(java.io.InputStream stream) {
	  ReInit(stream, null);
  }
  /** Reinitialise. */
  public void ReInit(java.io.InputStream stream, String encoding) {
	 try { jj_input_stream.ReInit(stream, encoding, 1, 1); } catch(java.io.UnsupportedEncodingException e) { throw new RuntimeException(e); }
	 token_source.ReInit(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 42; i++) jj_la1[i] = -1;
  }

  /** Constructor. */
  public SchemaParser(java.io.Reader stream) {
	 jj_input_stream = new SimpleCharStream(stream, 1, 1);
	 token_source = new SchemaParserTokenManager(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 42; i++) jj_la1[i] = -1;
  }

  /** Reinitialise. */
  public void ReInit(java.io.Reader stream) {
	if (jj_input_stream == null) {
	   jj_input_stream = new SimpleCharStream(stream, 1, 1);
	} else {
	   jj_input_stream.ReInit(stream, 1, 1);
	}
	if (token_source == null) {
 token_source = new SchemaParserTokenManager(jj_input_stream);
	}

	 token_source.ReInit(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 42; i++) jj_la1[i] = -1;
  }

  /** Constructor with generated Token Manager. */
  public SchemaParser(SchemaParserTokenManager tm) {
	 token_source = tm;
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 42; i++) jj_la1[i] = -1;
  }

  /** Reinitialise. */
  public void ReInit(SchemaParserTokenManager tm) {
	 token_source = tm;
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 42; i++) jj_la1[i] = -1;
  }

  private Token jj_consume_token(int kind) throws ParseException {
	 Token oldToken;
	 if ((oldToken = token).next != null) token = token.next;
	 else token = token.next = token_source.getNextToken();
	 jj_ntk = -1;
	 if (token.kind == kind) {
	   jj_gen++;
	   return token;
	 }
	 token = oldToken;
	 jj_kind = kind;
	 throw generateParseException();
  }


/** Get the next Token. */
  final public Token getNextToken() {
	 if (token.next != null) token = token.next;
	 else token = token.next = token_source.getNextToken();
	 jj_ntk = -1;
	 jj_gen++;
	 return token;
  }

/** Get the specific Token. */
  final public Token getToken(int index) {
	 Token t = token;
	 for (int i = 0; i < index; i++) {
	   if (t.next != null) t = t.next;
	   else t = t.next = token_source.getNextToken();
	 }
	 return t;
  }

  private int jj_ntk_f() {
	 if ((jj_nt=token.next) == null)
	   return (jj_ntk = (token.next=token_source.getNextToken()).kind);
	 else
	   return (jj_ntk = jj_nt.kind);
  }

  private java.util.List<int[]> jj_expentries = new java.util.ArrayList<int[]>();
  private int[] jj_expentry;
  private int jj_kind = -1;

  /** Generate ParseException. */
  public ParseException generateParseException() {
	 jj_expentries.clear();
	 boolean[] la1tokens = new boolean[52];
	 if (jj_kind >= 0) {
	   la1tokens[jj_kind] = true;
	   jj_kind = -1;
	 }
	 for (int i = 0; i < 42; i++) {
	   if (jj_la1[i] == jj_gen) {
		 for (int j = 0; j < 32; j++) {
		   if ((jj_la1_0[i] & (1<<j)) != 0) {
			 la1tokens[j] = true;
		   }
		   if ((jj_la1_1[i] & (1<<j)) != 0) {
			 la1tokens[32+j] = true;
		   }
		 }
	   }
	 }
	 for (int i = 0; i < 52; i++) {
	   if (la1tokens[i]) {
		 jj_expentry = new int[1];
		 jj_expentry[0] = i;
		 jj_expentries.add(jj_expentry);
	   }
	 }
	 int[][] exptokseq = new int[jj_expentries.size()][];
	 for (int i = 0; i < jj_expentries.size(); i++) {
	   exptokseq[i] = jj_expentries.get(i);
	 }
	 return new ParseException(token, exptokseq, tokenImage);
  }

  private boolean trace_enabled;

/** Trace enabled. */
  final public boolean trace_enabled() {
	 return trace_enabled;
  }

  /** Enable tracing. */
  final public void enable_tracing() {
  }

  /** Disable tracing. */
  final public void disable_tracing() {
  }

        }
