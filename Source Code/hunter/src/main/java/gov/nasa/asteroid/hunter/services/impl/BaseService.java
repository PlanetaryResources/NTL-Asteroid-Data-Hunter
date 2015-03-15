/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services.impl;

import gov.nasa.asteroid.hunter.LoggingHelper;

import javax.annotation.PostConstruct;

import org.apache.log4j.Logger;


/**
 * <p>
 * Base class of all service implementations. It holds a Logger.
 * </p>
 * <p>
 * <b>Thread Safety:</b> This class is effectively thread safe (injected
 * configurations are not considered as thread safety factor).
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public abstract class BaseService {
    /**
     * <p>
     * Represents the name of the class for logging.
     * </p>
     */
    private static final String CLASS_NAME = BaseService.class.getName();

    /**
     * <p>
     * Represents the Logger used to perform logging.
     * </p>
     *
     * <p>
     * Optional. If it is not configured, then no logging will be done in the service.
     * </p>
     */
    private Logger logger;



    /**
     * <p>
     * This is the default constructor for <code>BaseService</code>.
     * </p>
     */
    protected BaseService() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     *
     * Note, since there are no required fields in the base class, this function
     * always is successful.
     * </p>
     *
     * @throws ConfigurationException
     *             if any required field is not initialized properly.
     */
    @PostConstruct
    protected void checkConfiguration() {
        // just log the entrance and exit
        // prepare for logging
        final String signature = CLASS_NAME + ".checkConfiguration()";

        // log the entrance
        LoggingHelper.logEntrance(logger, signature, null, null);
        // log the exit
        LoggingHelper.logExit(logger, signature, null);
    }


    /**
     * <p>
     * Gets the logger for logging, can be null.
     * </p>
     * @return the logger for logging, can be null.
     */
    protected Logger getLogger() {
        return logger;
    }
    /**
     * <p>
     * Sets the logger.
     * </p>
     * @param logger the logger.
     */
    public void setLogger(Logger logger) {
        this.logger = logger;
    }
}

