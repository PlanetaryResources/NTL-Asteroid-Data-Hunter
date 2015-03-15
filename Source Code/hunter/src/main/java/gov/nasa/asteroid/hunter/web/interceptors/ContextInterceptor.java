/*
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */

package gov.nasa.asteroid.hunter.web.interceptors;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.web.servlet.ModelAndView;
import org.springframework.web.servlet.handler.HandlerInterceptorAdapter;

/**
 * <p>
 * The context intercepter to inject common context.
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
public class ContextInterceptor extends HandlerInterceptorAdapter {
    /**
     * <p>
     * The constructor.
     * </p>
     */
    public ContextInterceptor() {
        // Empty Constructor
    }

    /**
     * <p>
     * Add more attributes to the context
     * </p>
     * 
     * @param request
     *            the request
     * @param response
     *            the response
     * @param handler
     *            the handler
     * @param modelAndView
     *            model and view
     */
    @Override
    public void postHandle(HttpServletRequest request, HttpServletResponse response, Object handler,
            ModelAndView modelAndView) throws Exception {

        if (modelAndView != null) {
            modelAndView.addObject("viewName", modelAndView.getViewName());
        }
        super.postHandle(request, response, handler, modelAndView);
    }

}