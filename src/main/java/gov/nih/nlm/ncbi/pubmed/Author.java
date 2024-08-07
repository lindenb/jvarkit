//
// This file was generated by the Eclipse Implementation of JAXB, v4.0.5 
// See https://eclipse-ee4j.github.io/jaxb-ri 
// Any modifications to this file will be lost upon recompilation of the source schema. 
//


package gov.nih.nlm.ncbi.pubmed;

import java.util.ArrayList;
import java.util.List;
import jakarta.xml.bind.annotation.XmlAccessType;
import jakarta.xml.bind.annotation.XmlAccessorType;
import jakarta.xml.bind.annotation.XmlAttribute;
import jakarta.xml.bind.annotation.XmlElement;
import jakarta.xml.bind.annotation.XmlElements;
import jakarta.xml.bind.annotation.XmlRootElement;
import jakarta.xml.bind.annotation.XmlType;
import jakarta.xml.bind.annotation.adapters.CollapsedStringAdapter;
import jakarta.xml.bind.annotation.adapters.XmlJavaTypeAdapter;


/**
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName"
})
@XmlRootElement(name = "Author")
public class Author {

    @XmlAttribute(name = "ValidYN")
    @XmlJavaTypeAdapter(CollapsedStringAdapter.class)
    protected String validYN;
    @XmlElements({
        @XmlElement(name = "LastName", required = true, type = LastName.class),
        @XmlElement(name = "ForeName", required = true, type = ForeName.class),
        @XmlElement(name = "Initials", required = true, type = Initials.class),
        @XmlElement(name = "Suffix", required = true, type = Suffix.class),
        @XmlElement(name = "NameID", required = true, type = NameID.class),
        @XmlElement(name = "CollectiveName", required = true, type = CollectiveName.class)
    })
    protected List<java.lang.Object> lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName;

    /**
     * Gets the value of the validYN property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getValidYN() {
        if (validYN == null) {
            return "Y";
        } else {
            return validYN;
        }
    }

    /**
     * Sets the value of the validYN property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setValidYN(String value) {
        this.validYN = value;
    }

    /**
     * Gets the value of the lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName property.
     * 
     * <p>This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName property.</p>
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * </p>
     * <pre>
     * getLastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link CollectiveName }
     * {@link ForeName }
     * {@link Initials }
     * {@link LastName }
     * {@link NameID }
     * {@link Suffix }
     * </p>
     * 
     * 
     * @return
     *     The value of the lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName property.
     */
    public List<java.lang.Object> getLastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName() {
        if (lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName == null) {
            lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName = new ArrayList<>();
        }
        return this.lastNameOrForeNameOrInitialsOrSuffixOrNameIDOrCollectiveName;
    }

}
