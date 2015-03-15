/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.web.controllers;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;

/**
 * <p>
 * Base class for controller implementations.
 * </p>
 * <p>
 * <strong>Thread Safety: </strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * </p>
 * 
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
@Controller
public abstract class BaseController {
    
    /**
     * <p>
     * Represents the Logger used to perform logging.
     * </p>
     * <p>
     * Optional. If it is not configured, then no logging will be done in the
     * service.
     * </p>
     */
    @Autowired
    private Logger logger;

    /**
     * <p>
     * This is the default controller for base controller.
     * </p>
     */
    protected BaseController() {
        // does nothing
    }

    /**
     * <p>
     * Check if all required fields are initialized properly.
     * </p>
     * 
     * @throws ConfigurationException
     *             if any required field is not initialized properly.
     */
    protected void checkConfiguration() {
        // does nothing
    }

    /**
     * <p>
     * Gets the logger.
     * </p>
     * 
     * @return the logger.
     */
    protected Logger getLogger() {
        return logger;
    }
}
