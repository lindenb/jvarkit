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
import jakarta.xml.bind.annotation.XmlElement;
import jakarta.xml.bind.annotation.XmlRootElement;
import jakarta.xml.bind.annotation.XmlType;


/**
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "personalNameSubject"
})
@XmlRootElement(name = "PersonalNameSubjectList")
public class PersonalNameSubjectList {

    @XmlElement(name = "PersonalNameSubject", required = true)
    protected List<PersonalNameSubject> personalNameSubject;

    /**
     * Gets the value of the personalNameSubject property.
     * 
     * <p>This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the personalNameSubject property.</p>
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * </p>
     * <pre>
     * getPersonalNameSubject().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link PersonalNameSubject }
     * </p>
     * 
     * 
     * @return
     *     The value of the personalNameSubject property.
     */
    public List<PersonalNameSubject> getPersonalNameSubject() {
        if (personalNameSubject == null) {
            personalNameSubject = new ArrayList<>();
        }
        return this.personalNameSubject;
    }

}