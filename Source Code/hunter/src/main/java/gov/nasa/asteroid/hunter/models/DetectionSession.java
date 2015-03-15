/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Table;

/**
 * Represents a detection session.
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Entity
@Table(name = "detection_session")
public class DetectionSession extends IdentifiableEntity {
    /**
     * Represents the detection session progress.
     *
     * Required.
     */
    @Column(name = "progress", nullable = false)
    private float progress;

    /**
     * <p>
     * The default constructor for class <code>DetectionSession</code>.
     * </p>
     */
    public DetectionSession() {
        // does nothing
    }

    /**
     * <p>
     * Gets the detection session progress.
     * </p>
     * @return the detection session progress.
     */
    public float getProgress() {
        return progress;
    }

    /**
     * <p>
     * Sets the detection session progress.
     * </p>
     * @param progress the progress of the session.
     */
    public void setProgress(float progress) {
        this.progress = progress;
    }
}

