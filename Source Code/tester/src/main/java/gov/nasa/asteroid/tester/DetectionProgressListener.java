/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.tester;

/**
 * <p>
 * This interface defines the contract for a detection progress listener.
 * </p>
 * <p>
 * <strong>Thread Safety:</strong> Implementation does not have to be thread
 * safe.
 * </p>
 * 
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public interface DetectionProgressListener {
    
    /**
     * <p>
     * Actions on progress updated.
     * </p>
     * @param progress the progress.
     */
    public void updateProgress(float progress);
}

