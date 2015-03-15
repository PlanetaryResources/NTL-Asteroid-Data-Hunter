/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter;


/**
 * <p>
 * This exception will be thrown to indicate any configuration error.
 * </p>
 *
 * <p>
 * <strong>Thread Safety: </strong> Thread Safety is not an issue for Exceptions
 * since exceptions is used only in one thread.
 * </p>
 *
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class ConfigurationException extends RuntimeException {

    /**
     * <p>
     * This is the constructor of class <code>ConfigurationException</code> with message.
     * </p>
     * @param message the error message
     */
    public ConfigurationException(String message) {
        super(message);
    }

    /**
     * <p>
     * This is the constructor of class <code>ConfigurationException</code> with message
     * and error cause.
     * </p>
     * @param message the error message
     * @param cause the error cause.
     */
    public ConfigurationException(String message, Throwable cause) {
        super(message, cause);
    }
}

