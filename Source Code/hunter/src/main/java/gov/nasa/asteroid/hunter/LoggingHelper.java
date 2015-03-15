/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package gov.nasa.asteroid.hunter;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import org.apache.log4j.Logger;


/**
 * <p>
 * This is the helper class of this assembly.
 * </p>
 *
 * <strong>Thread Safety: </strong> this class is immutable and thread-safe.
 *
 * @author TCSASSEMBLER
 * @version 1.0
 */
public class LoggingHelper {

    /**
     * <p>
     * Represents the entrance message.
     * </p>
     */
    private static final String MESSAGE_ENTRANCE = "Entering method %1$s.";

    /**
     * <p>
     * Represents the exit message.
     * </p>
     */
    private static final String MESSAGE_EXIT = "Exiting method %1$s.";

    /**
     * <p>
     * Represents the error message.
     * </p>
     */
    private static final String MESSAGE_ERROR = "Error in method %1$s. Details:";

    /**
     * <p>
     * Private empty constructor to prevent instantiating this class.
     * </p>
     */
    private LoggingHelper() {
        // does nothing
    }

    /**
     * <p>
     * Logs for entrance into public methods at <code>DEBUG</code> level.
     * </p>
     *
     * @param logger
     *            the log service.
     * @param signature
     *            the signature.
     * @param paramNames
     *            the names of parameters to log (not Null).
     * @param params
     *            the values of parameters to log (not Null).
     */
    public static void logEntrance(Logger logger, String signature, String[] paramNames, Object[] params) {
        if (logger == null) {
            // No logging
            return;
        }

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(MESSAGE_ENTRANCE, signature));
        sb.append(toString(paramNames, params));
        logger.debug(sb.toString());
    }

    /**
     * <p>
     * Logs for exit from public methods at <code>DEBUG</code> level.
     * </p>
     *
     * @param logger
     *            the log service.
     * @param signature
     *            the signature of the method to be logged.
     * @param value
     *            the return value to log.
     */
    public static void logExit(Logger logger, String signature, Object[] value) {
        if (logger == null) {
            // No logging
            return;
        }

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(MESSAGE_EXIT, signature));
        if (value != null) {
            sb.append(" Output parameter : " + toString(value[0]));
        }
        logger.debug(sb.toString());
    }

    /**
     * <p>
     * Logs the given exception and message at <code>ERROR</code> level.
     * </p>
     *
     * @param <T>
     *            the exception type.
     * @param logger
     *            the log service.
     * @param signature
     *            the signature of the method to log.
     * @param e
     *            the exception to log.
     *
     * @return the passed in exception.
     */
    public static <T extends Throwable> T logException(Logger logger, String signature, T e) {
        if (logger != null) {
            StringWriter sw = new StringWriter();
            e.printStackTrace(new PrintWriter(sw));
            logger.error(String.format(MESSAGE_ERROR, signature + " " + sw.toString()));
        }
        return e;
    }

    /**
     * <p>
     * Converts the parameters to string.
     * </p>
     *
     * @param paramNames
     *            the names of parameters.
     * @param params
     *            the values of parameters.
     * @return the string
     */
    static String toString(String[] paramNames, Object[] params) {
        StringBuffer sb = new StringBuffer(" Input parameters: {");
        if (params != null) {
            for (int i = 0; i < params.length; i++) {
                if (i > 0) {
                    sb.append(", ");
                }
                sb.append(paramNames[i]).append(":").append(toString(params[i]));
            }
        }
        sb.append("}.");
        return sb.toString();
    }

    /**
     * <p>
     * Converts the object to string.
     * </p>
     *
     * @param obj
     *            the object
     *
     * @return the string representation of the object.
     */
    static String toString(Object obj) {
        if (obj instanceof List<?>) {
            StringBuffer str = new StringBuffer("[");
            for (Object item : (List<?>) obj) {
                if (!str.toString().equals("[")) {
                    str.append(",");
                }
                str.append(toString(item));
            }
            str.append("]");
            return str.toString();
        } else {
            return String.valueOf(obj);
        }
    }

}
