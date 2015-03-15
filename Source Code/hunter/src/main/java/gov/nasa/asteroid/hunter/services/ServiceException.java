/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter.services;

/**
 * <p>
 * This exception will be thrown by services to indicate that there are some errors during the operations.
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
public class ServiceException extends Exception {
    /**
     * <p>
     * This is the constructor of class <code>ServiceException</code> with message.
     * </p>
     * @param message the error message
     */
    public ServiceException(String message) {
        super(message);
    }

    /**
     * <p>
     * This is the constructor of class <code>ServiceException</code> with message
     * and error cause.
     * </p>
     * @param message the error message
     * @param cause the error cause.
     */
    public ServiceException(String message, Throwable cause) {
        super(message, cause);
    }
}

