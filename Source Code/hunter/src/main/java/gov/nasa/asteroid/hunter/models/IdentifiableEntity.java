/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.MappedSuperclass;

/**
 * Represents an identifiable entity.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@MappedSuperclass
public abstract class IdentifiableEntity {

    /**
     * Represents the ID.
     */
    @Id
    @GeneratedValue (strategy = GenerationType.IDENTITY)
    private long id;

    /**
     * <p>
     * The default constructor for class <code>IdentifiableEntity</code>.
     * </p>
     */
    protected IdentifiableEntity() {
        // does nothing
    }

    /**
     * <p>
     * Gets the id.
     * </p>
     * @return the id.
     */
    public long getId() {
        return id;
    }

    /**
     * <p>
     * Sets the id.
     * </p>
     * @param id the id.
     */
    public void setId(long id) {
        this.id = id;
    }
}
