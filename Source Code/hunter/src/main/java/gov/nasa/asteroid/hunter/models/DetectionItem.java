/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.models;

import java.util.Date;
import java.util.List;

import javax.persistence.Basic;
import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.OneToMany;
import javax.persistence.OrderBy;
import javax.persistence.Table;

/**
 * Represents a detection item.
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
@Table(name = "detection_item")
public class DetectionItem extends IdentifiableEntity {


    /**
     * Represents the timestamp.
     *
     * Required.
     */
    @Basic
    @Column(nullable = false)
    private Date timestamp;

    /**
     * Indicates whether the item is a NEO.
     *
     * Required.
     */
    @Basic
    @Column(nullable = false)
    private boolean neo;

    /**
     * Indicates whether the item is submitted.
     *
     * Required.
     */
    @Basic
    @Column(nullable = false)
    private boolean submitted;

    /**
     * Represents the frames data.
     *
     * Required.
     */
    @OneToMany (mappedBy = "detectionItem", cascade = CascadeType.ALL, fetch = FetchType.EAGER)
    @OrderBy("detectionItemFramePk.frame")
    private List<DetectionItemFrame> frames;

    /**
     * Indicates whether the item is a known by MPC.
     *
     * Required.
     */
    @Column(name = "known_by_mpc", nullable = false)
    private boolean knownByMPC;

    @Column(name = "significance", nullable = false)
    private double significance;
    
    @Column(name = "observatory_code", nullable = false)
    private String observatoryCode;
    
    @Column(name = "image_width", nullable = false)
    private double imageWidth;
    
    @Column(name = "image_height", nullable = false)
    private double imageHeight;
    
    @Column(name = "image_id", nullable = false)
    private long imageId;
    
    /**
     * <p>
     * The default constructor for class <code>DetectionItem</code>.
     * </p>
     */
    public DetectionItem() {
        // does nothing
    }

    /**
     * <p>
     * Gets the timestamp.
     * </p>
     * @return the timestamp.
     */
    public Date getTimestamp() {
        return timestamp;
    }

    /**
     * <p>
     * Sets the timestamp.
     * </p>
     * @param timestamp the timestamp.
     */
    public void setTimestamp(Date timestamp) {
        this.timestamp = timestamp;
    }

    /**
     * <p>
     * Gets whether the item is a NEO.
     * </p>
     * @return whether the item is a NEO.
     */
    public boolean isNeo() {
        return neo;
    }

    /**
     * <p>
     * Sets whether the item is a NEO.
     * </p>
     * @param neo whether the item is a NEO.
     */
    public void setNeo(boolean neo) {
        this.neo = neo;
    }

    /**
     * <p>
     * Gets whether the item is submitted.
     * </p>
     * @return whether the item is submitted.
     */
    public boolean isSubmitted() {
        return submitted;
    }

    /**
     * <p>
     * Sets whether the item is submitted.
     * </p>
     * @param submitted whether the item is submitted.
     */
    public void setSubmitted(boolean submitted) {
        this.submitted = submitted;
    }

    /**
     * <p>
     * Gets the frames data.
     * </p>
     * @return the frames data.
     */
    public List<DetectionItemFrame> getFrames() {
        return frames;
    }

    /**
     * <p>
     * Sets the frames data.
     * </p>
     * @param frames the frames data.
     */
    public void setFrames(List<DetectionItemFrame> frames) {
        this.frames = frames;
    }

    /**
     * <p>
     * Gets whether the item is a known by MPC.
     * </p>
     * @return whether the item is a known by MPC.
     */
    public boolean isKnownByMPC() {
        return knownByMPC;
    }

    /**
     * <p>
     * Sets whether the item is a known by MPC.
     * </p>
     * @param knownByMPC whether the item is a known by MPC.
     */
    public void setKnownByMPC(boolean knownByMPC) {
        this.knownByMPC = knownByMPC;
    }

    public double getSignificance() {
        return significance;
    }

    public void setSignificance(double significance) {
        this.significance = significance;
    }

    public String getObservatoryCode() {
        return observatoryCode;
    }

    public void setObservatoryCode(String observatoryCode) {
        this.observatoryCode = observatoryCode;
    }

    public double getImageWidth() {
        return imageWidth;
    }

    public void setImageWidth(double imageWidth) {
        this.imageWidth = imageWidth;
    }

    public double getImageHeight() {
        return imageHeight;
    }

    public void setImageHeight(double imageHeight) {
        this.imageHeight = imageHeight;
    }

    public long getImageId() {
        return imageId;
    }

    public void setImageId(long imageId) {
        this.imageId = imageId;
    }
}

