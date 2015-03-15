<%--
  - Author: TCSASSEMBLER
  - Version: 1.0
  - Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
  -
  - Description: The help items page.
--%>
<%@taglib uri="http://www.springframework.org/tags/form" prefix="form" %>

<%@ page import="com.fasterxml.jackson.databind.ObjectMapper" %>

<%@include file="header.jsp" %>

<script>
    var gTopics = <% out.print(new ObjectMapper().writeValueAsString(request.getAttribute("helpTopics"))); %>
    var gHelpItem = <% out.print(new ObjectMapper().writeValueAsString(request.getAttribute("helpItem"))); %>
    var gViewName = "${viewName}"
</script>
        <div class="page-content help-page">
            <div class="round-box help-page-content">
                <div class="title-panel">
                    <img src='../i/Help-icon.png' /><span class="title-label">HELP</span>
                    <div class="search-component">
                        SEARCH:
                        <div class="search-box">
                            <input type="text" />
                            <img class="help-search-btn" src = "../i/search-blue.png" />
                        </div>
                    </div>
                </div>
                <div class="seperator"></div>
                <div class="topics-and-items">
                    <div class="topics-panel">
                        <div class="topics-panel-title blue-text">Help Topics</div>
                        <div class="seperator"></div>
                        <div class="topics-list">
                            <div class="topic-item selected">
                                <span class="topic-name">Getting Started</span>
                            </div>
                            <div class="seperator"></div>
                            <div class="topic-item">
                                <span class="topic-name">Application</span>
                            </div>
                            <div class="seperator"></div>
                        </div>
                    </div>
                    <div class="items-panel content-panel">
                        <div class="item-topic-title blue-text">Getting Started</div>
                        <div class="item-list asteroid-table-container">
                            <table>
                                <thead>
                                    <tr class="hidden">
                                        <th>Item Name</th>
                                    </tr>
                                </thead>
                                <tbody>
        
                                </tbody>
                            </table>
                            <div class="page-section"></div>
                        </div>
                    </div>
                
                    <div class="item-content-panel content-panel hidden">
                        <div class="item-topic-title blue-text">Getting Started</div>
                        <div class="item-content"></div>
                    </div>
                    
                    <div class="search-result-panel content-panel hidden">
                        <div class="item-topic-title blue-text">Search Results For &quot;<span class="search-term-title">abc</span>&quot;</div>
                        <div class="items-result-table asteroid-table-container">
                            <table>
                                <thead>
                                    <tr>
                                        <th>TOPIC</th>
                                        <th>HELP ITEM</th>
                                    </tr>
                                </thead>
                                <tbody>
        
                                </tbody>
                            </table>
                            <div class="page-section"></div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
<%@include file="footer.jsp" %>