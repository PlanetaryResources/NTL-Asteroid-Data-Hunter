/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.Helper;
import gov.nasa.asteroid.hunter.models.DetectionSession;
import gov.nasa.asteroid.hunter.services.DetectionSessionService;
import gov.nasa.asteroid.hunter.services.EntityNotFoundException;
import gov.nasa.asteroid.hunter.services.ServiceException;

import org.apache.log4j.Logger;
import org.springframework.transaction.annotation.Propagation;
import org.springframework.transaction.annotation.Transactional;

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
public class DetectionSessionServiceImpl  extends BasePersistenceService implements DetectionSessionService{
    /**
     * <p>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = DetectionSessionServiceImpl.class.getName();

    /**
     * <p>
     * Creates the instance of DetectionSessionServiceImpl.
     * </p>
     */
    public DetectionSessionServiceImpl() {
        // does nothing
    }
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
    @Override
    public DetectionSession getDetectionSession(long id) throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".getDetectionSession(long)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });

        // check the parameters
        Helper.checkPositive(logger, signature, "id", id);

        try {
            DetectionSession result = getEntityManager().find(DetectionSession.class, id);

            // log the exit
            LoggingHelper.logExit(logger, signature, new Object[] { result });

            return result;
        } catch (IllegalArgumentException e) {
            throw LoggingHelper.logException(logger, signature, new ServiceException("Failed to get the result.", e));
        }
    }

    /**
     * <p>
     * This method is used to create a new detection session.
     * </p>
     *
     * @return the detection session.
     *
     * @throws ServiceException if any other error occurred during the operation.
     */
    @Transactional
    @Override
    public DetectionSession newDetectionSession() throws ServiceException {
        // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".newDetectionSession()";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);

        // create the session
        DetectionSession session = new DetectionSession();

        // persist it
        persist(logger, signature, session);

        // log the exit
        LoggingHelper.logExit(logger, signature, new Object[] { session });

        return session;
    }

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
    @Transactional(propagation = Propagation.REQUIRES_NEW)
    @Override
    public void updateDetectionSessionProgress(long id, float progress) throws ServiceException {
     // prepare for logging
        Logger logger = getLogger();
        final String signature = CLASS_NAME + ".updateDetectionSessionProgress(long, progress)";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, new String[] { "id" }, new Object[] { id });

        // check the parameters
        Helper.checkPositive(logger, signature, "id", id);
        if (Double.compare(progress, 0) < 0 || Double.compare(progress, 1) > 0) {
            throw LoggingHelper.logException(
                    logger, signature, new IllegalArgumentException("The progress is invalid."));
        }
        try {
            DetectionSession session = getDetectionSession(id);
            if (session == null) {
                throw LoggingHelper.logException(logger, signature, new EntityNotFoundException(
                        "Detection Session with id" + id + " does not exist."));
            }
            // mark it
            session.setProgress(progress);

            // update the item.
            merge(logger, signature, session);

            // log the exit
            LoggingHelper.logExit(logger, signature, null);

        } catch (ServiceException e) {
            // log and re-throw
            throw LoggingHelper.logException(logger, signature, e);
        }
    }
}
