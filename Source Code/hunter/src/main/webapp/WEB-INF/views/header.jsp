<%--
  - Author: TCSASSEMBLER
  - Version: 1.0
  - Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
  -
  - Description: The header of the jsp pages
--%>
<%@ taglib prefix="c" uri="http://java.sun.com/jsp/jstl/core" %>
<%@ taglib prefix="fn" uri="http://java.sun.com/jsp/jstl/functions"%>

<!DOCTYPE html>

<html>
<head>

        
        <!-- title -->
        <title>Asteroid Data Hunter</title>
        
        <!-- metatags -->
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />

		<link rel="shortcut icon" href="<c:url value="/i/favicon.ico"/>">
		
        <!-- CSS -->
        <link href='//fonts.googleapis.com/css?family=Open+Sans:400,300,700' rel='stylesheet' type='text/css'>
        <link rel="stylesheet" href="<c:url value="/css/style.css"/>" media="all" />
        <!--[if IE 9]>
        <link rel="stylesheet" href="<c:url value="/css/style-ie9.css"/>" media="all" />
        <![endif]-->

        <link rel="stylesheet" href="<c:url value="/css/anytime.5.0.5.min.css"/>" media="all" />
        <link rel="stylesheet" href="<c:url value="/css/jcarousel.css"/>" media="all" />
        <link rel="stylesheet" href="<c:url value="/css/jquery-ui.css"/>" media="all" />
        
        <!-- JS -->
        <script type="text/javascript" src="<c:url value="/js/jquery-1.11.1.min.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery-dateFormat.min.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery.jcarousel.min.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery.vintage.min.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/anytime.5.0.5.min.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery-ui.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery.cookie.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery.elevatezoom.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/jquery.visible.min.js"/>"></script>
        <script type="text/javascript" src="<c:url value="/js/comm.js"/>"></script>

        <!-- js for each page -->
        <script type="text/javascript" src=""></script>
        
        <c:if test="${viewName == 'viewDashboard' || viewName == 'viewDetectionItem'}">
            <script type="text/javascript" src="<c:url value="/js/dashboard.js"/>"></script>
        </c:if>
        <c:if test="${viewName == 'listDetectionItems'}">
            <script type="text/javascript" src="<c:url value="/js/search.js"/>"></script>
        </c:if>
        <c:if test="${viewName == 'listNEOs'}">
            <script type="text/javascript" src="<c:url value="/js/neos.js"/>"></script>
        </c:if>
        <c:if test="${viewName == 'listHelpItems' || viewName == 'viewHelpItem'}">
            <script type="text/javascript" src="<c:url value="/js/topic.js"/>"></script>
        </c:if>
</head>

<body>
<div class="all-content">
<div class="center-content">
    <div class="header-content">
            <div class="logo-section">
                <a><img src="../i/logo.png" /></a>
                <div class="header-title">Asteroid Data Hunter</div>
            </div>
            <div class="header-tabs">
                <div class="v-seperator"></div>
                <div class="header-tab dashboard-tab ${viewName == 'viewDashboard' ? 'selected' : ''}" data-page="dashboard"><img src="../i/dashboard-icon.png"/>DASHBOARD</div>
                <div class="v-seperator"></div>
                <div class="header-tab search-tab ${viewName == 'listDetectionItems' || viewName == 'viewDetectionItem' ? 'selected' : ''}" data-page="detectionItems"><img src="../i/search-icon.png"/>ID List</div>
                <div class="v-seperator"></div>
                <div class="header-tab sky-tab ${viewName == 'listNEOs' ? 'selected' : ''}" data-page="neos"><img src="../i/sky-icon.png"/>SKY</div>
                <div class="v-seperator"></div>
                <div class="hidden header-tab right-seperator helper-tab ${(viewName == 'listHelpItems' || viewName == 'viewHelpItem') ? 'selected' : ''}" data-page="helpItems"><img src="../i/Help-icon.png"/>HELP</div>
                <div class="hidden v-seperator"></div>
            </div>
        
    </div>
    <div class="seperator">
    </div>
    <div class="main-content">
