<%--
  - Author: TCSASSEMBLER
  - Version: 1.0
  - Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
  -
  - Description: The detection item details page.
--%>
<%@ page import="com.fasterxml.jackson.databind.ObjectMapper" %>

<%@include file="header.jsp" %>

<script>
    var latestItem = <% out.print(new ObjectMapper().writeValueAsString(request.getAttribute("detectionItem"))); %>
    var gSearchResult = null
    var observationsMailAddress =  "${observationsMailAddress}"
    var observationsMailSubject = "${observationsMailSubject}"
    var observationsMailBody = "${observationsMailBody}"
</script>
        <div class="page-content dashboard-page">
            <div class="round-box latest-result">
                <div class="latest-result-title"><img src="../i/latest-icon.png"/>SEARCH DETAILS<div class="constrast-slider">CONTRAST&nbsp;&nbsp;<div class="slider" ></div></div></div>
                <div class="latest-result-content">

                    <div class="result-img box">     
                    	<c:if test="${detectionItem != null}">               
                        <div class="img-container">
                        	<img class="bg-img" data-frame-index="0" src="/frameimage/${detectionItem.id}/frame/0/image" alt="">
                        	<img class="red-circle" src="../i/red-circle.png" />
                        </div>  
                        </c:if>  
                    </div>

                    <div class="result-info box">
                        <div class="info-item">
                            <div class="info-item-name">NEO?</div>
                            <div class="info-item-value blue-text">Asteroid</div>
                            <div class="seperator"></div>
                        </div>
                        <div class="info-item">
                            <div class="info-item-name">NEO?</div>
                            <div class="info-item-value blue-text">Asteroid</div>
                            <div class="seperator"></div>
                        </div>
                    </div>

                    <div class="info-item-template hidden">
                        <div class="info-item-name">NEO?</div>
                        <div class="info-item-value blue-text">Asteroid</div>
                        <div class="seperator"></div>
                    </div>

                    <div class="result-controls">
                    <a href="javascript:void(0)" class="play-pause-btn" title="Will play all 4 images in order"><img class="play-btn" src="../i/play-btn.png"/><img class="pause-btn hidden" src="../i/pause-btn.png"/></a>
                        <div class="frame-control">
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="0">1</a>
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="1">2</a>
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="2">3</a>
                            <a href="javascript:void(0);" class="frame-index-btn" data-frame-index="3">4</a>
                        </div>
                        <a href="javascript:void(0)" class="invert-btn" title="Will invert the image"><img src="../i/invert-btn.png"/></a>
                        <a href="javascript:void(0)" class="flicker-btn hidden"><img src="../i/circle-btn.png"/></a>
                    </div>
                </div>
                <div class="latest-result-no-content hidden">
                    <p>No Recent Results</p>
                </div>
            </div>
            <div class="round-box submit-panel">
                <div class="submit-section-title"><img src="../i/submit-mpc.png"/>SUBMIT TO MPC</div>
                <div class="submit-section-note red-text">Are you really sure to submit this?</div>
                <div class="submit-item round-box">
                    <div class="submit-panel-title">OBSERVATIONS</div>
                    <div class="submit-panel-content">If a signifigant detection has been made, and does NOT exist at the Minor Planet Center (refer to the &quot;Exists at MPC&quot; field to the left), then please use the button below to submit.</div>
                    <div class="send-email-btn"><a href="javascript:;"><img src="../i/email.png"/>${observationsEmailAddress}</a></div>
                </div>
                <div class="submit-item round-box hidden">
                    <div class="submit-panel-title">NEW COMET REPORTS</div>
                    <div class="submit-panel-content">Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed id justo venenatis, congue dui vel, placerat elit. Maecenas tempor, eros rutrum elementum vestibulum, sapien orci mollis justo, quis auctor elit ante et quam. Nulla a pretium quam, mattis sagittis nisi. Vivamus pharetra.</div>
                    <div class="send-email-btn"><a href="mailto:${newCometReportMailTo}"><img src="../i/email.png"/>${newCometReportsEmailAddress}</a></div>
                </div>
            </div>
        </div>

<%@include file="footer.jsp" %>