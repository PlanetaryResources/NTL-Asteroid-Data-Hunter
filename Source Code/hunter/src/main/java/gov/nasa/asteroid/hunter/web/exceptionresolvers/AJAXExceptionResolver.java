package gov.nasa.asteroid.hunter.web.exceptionresolvers;

import gov.nasa.asteroid.hunter.LoggingHelper;
import gov.nasa.asteroid.hunter.services.EntityNotFoundException;

import java.io.IOException;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.servlet.ModelAndView;
import org.springframework.web.servlet.handler.SimpleMappingExceptionResolver;

/**
 * <p>
 * Spring Resolver for AJAX.
 * </p>
 * 
 * <p>
 * <strong>Thread Safety:</strong> This class is effectively thread safe
 * (injected configurations are not considered as thread safety factor).
 * </p>
 * 
 * @author albertwang, TCSASSEMBLER
 * @version 1.0
 */
public class AJAXExceptionResolver extends SimpleMappingExceptionResolver {
    
    /**
     * <p>
     * Represents the class name.
     * </p>
     */
    private static final String CLASS_NAME = AJAXExceptionResolver.class.getName();
    
    /**
     * Represents the Logger used to perform logging.
     * 
     * Optional. If it is not configured, then no logging will be done in the service.
     */
    @Autowired
    private Logger logger;

    /**
     * <p>
     * The default constructor of this class.
     * </p>
     */
    public AJAXExceptionResolver() {
        // does nothing
    }

    /**
     * Resolves the exception and sets proper http code.
     * 
     * @param request
     *            the http servlet request
     * @param response
     *            the http servlet response
     * @param handler
     *            the handler object
     * @param exception
     *            the thrown exception
     * 
     * @return Resulting ModelAndView to show to user
     */
    public ModelAndView resolveException(HttpServletRequest request, HttpServletResponse response, Object handler,
            Exception exception) {
        String signature = CLASS_NAME + "#resolveException(HttpServletRequest request, HttpServletResponse response,"
                + " Object handler, Exception exception)";

        LoggingHelper.logEntrance(logger, signature, new String[] { "request", "response", "handler", "exception" },
                new Object[] { request, response, handler, exception });

        ModelAndView result;
        // Perform the check only for AJAX requests
        if ("XMLHttpRequest".equals(request.getHeader("X-Requested-With"))) {
            try {
                String message = exception.getMessage();
                if (exception instanceof IllegalArgumentException) {
                    response.sendError(HttpServletResponse.SC_BAD_REQUEST, exception.getMessage());
                } else if (exception instanceof EntityNotFoundException) {
                    response.sendError(HttpServletResponse.SC_NOT_FOUND, exception.getMessage());
                } else {
                    response.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR, exception.getMessage());
                }
                // Write exception message to response
                response.getWriter().print(message);
            } catch (IOException e) {
                // Log and ignore the exception
                LoggingHelper.logException(logger, signature, e);
            }
            result = new ModelAndView();
        } else {
            return super.resolveException(request, response, handler, exception);
        }

        LoggingHelper.logExit(logger, signature, new Object[] { result });
        return result;
    }

}

