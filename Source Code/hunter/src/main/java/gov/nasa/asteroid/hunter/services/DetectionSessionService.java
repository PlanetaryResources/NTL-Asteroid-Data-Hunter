/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services;

import gov.nasa.asteroid.hunter.models.DetectionSession;

/**
 * <p>
 * This service provides methods to manage the detection sessions.
 * </p>
 * <p>
 * <b>Thread Safety:</b> Implementations must be effectively thread safe. Refer
 * to ADS 1.3.4 for general assumptions and notes on thread safety.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public interface DetectionSessionService {

    /**
     * <p>
     * This method is used to get detection session.
     * </p>
     *
     * @param id the ID
     *
     * @return the detection session, null will be returned if there's no such entity.
     *
     * @throws IllegalArgumentException if id is not positive
     * @throws ServiceException if any other error occurred during the operation
     */
    public DetectionSession getDetectionSession(long id) throws ServiceException;

    /**
     * <p>
     * This method is used to create a new detection session.
     * </p>
     *
     * @return the detection session.
     *
     * @throws ServiceException if any other error occurred during the operation.
     */
    public DetectionSession newDetectionSession() throws ServiceException;

    /**
     * <p>
     * This method is used to update the progress of the session
     * </p>
     *
     * @param id the session id
     * @param progress the new progress
     *
     * @throws ServiceException if any other error occurred during the operation.
     * @throws EntityNotFoundException if the session of the id does not exists.
     * @throws IllegaArgumentException if the id or the progress is invalid.
     */
    public void updateDetectionSessionProgress(long id, float progress) throws ServiceException;

}
