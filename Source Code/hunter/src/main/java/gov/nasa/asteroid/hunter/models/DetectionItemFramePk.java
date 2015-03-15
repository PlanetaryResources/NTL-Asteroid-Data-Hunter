/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import java.io.Serializable;

import javax.persistence.Column;
import javax.persistence.Embeddable;

/**
 * <p>
 * This class is the composite primary key for the DetectionItemFrame entity.
 * </p>
 *
 * <p>
 * <strong>Thread-Safety:</strong> This class is not thread safe since it is
 * mutable.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Embeddable
public class DetectionItemFramePk implements Serializable {

    /**
     * <p>
     * Represents the detectionItemId in the composite primary key.
     * </p>
     * <p>
     * Required.
     * </p>
     */
    @Column(name = "detectionItemId", nullable = false)
    private long detectionItemId;

    /**
     * <p>
     * Represents the frame in the composite primary key.
     * </p>
     * <p>
     * Required.
     * </p>
     */
    @Column(name = "frame", nullable = false)
    private int frame;

    /**
     * <p>
     * The default constructor for class <code>DetectionItemFramePk</code>.
     * </p>
     */
    public DetectionItemFramePk() {
        // does nothing
    }

    /**
     * <p>
     * Gets the frame in the composite primary key.
     * </p>
     * @return the frame in the composite primary key.
     */
    public int getFrame() {
        return frame;
    }

    /**
     * <p>
     * Sets the frame in the composite primary key.
     * </p>
     * @param frame the frame in the composite primary key.
     */
    public void setFrame(int frame) {
        this.frame = frame;
    }

    /**
     * <p>
     * Gets the detectionItemId in the composite primary key.
     * </p>
     * @return the detectionItemId in the composite primary key.
     */
    public long getDetectionItemId() {
        return detectionItemId;
    }

    /**
     * <p>
     * Sets the detectionItemId in the composite primary key.
     * </p>
     * @param detectionItemId the detectionItemId in the composite primary key.
     */
    public void setDetectionItemId(long detectionItemId) {
        this.detectionItemId = detectionItemId;
    }

    /**
     * <p>
     * The hash code of the entity.
     * </p>
     * @return the hash code of the instance
     */
    @Override
    public int hashCode() {
        return (int) (detectionItemId * 31 + frame);
    }

    /**
     * <p>
     * Checks if the entities are equal.
     * </p>
     * @param obj the object to compare.
     * @return true if the 2 objects are equal, otherwise false.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof DetectionItemFramePk)) {
            return false;
        }
        DetectionItemFramePk theObj = (DetectionItemFramePk) obj;
        return theObj.detectionItemId == this.detectionItemId
                && theObj.frame == this.frame;
    }

}
