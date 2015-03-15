/**
 * The javascript using in dash board page.
 * @author TCSASSEMBLER
 * @version 1.0
 */

$(document).ready(function() {

	    /* file upload */
    $(document).on('change', '#fileUpload', function (e) {
        var files = $("#fileUpload")[0].files;
        var text = ""
		for (var i = 0; i < files.length; i++) {
			if (text != "") {
				text += ", "
				
			}
			text += files[i].name
		}
    
        $('#fileUploadTxt').val(text);
        e.stopPropagation();
    });


    var uploader = new HunterUploader()
    uploader.progressCallback = function(session) {
        if (currentProgress + 1e-10 >= 1) {
            return
        }
        $.cookie("currentUploadSessionID", session.id);
        $.cookie("currentUploadProgress", session.progress);
    	setProgressView(session.progress)
    }
    var currentProgress = 0
    function setProgressView(progress) {
        if (currentProgress + 1e-10 >= 1) {
            return
        }
        currentProgress = progress
        $('.progress-bar').animate({
            'width': (progress * 100) + '%'
        }, 300, function() {
            var offset = $('.progress-bar').offset()
            $('.progress-tip').text(parseInt(progress * 100) + '%')
            $('.progress-tip').css({top:offset.top - 30, left:offset.left + $('.progress-bar').width() - 25})
            if (progress + 1e-10 >= 1) {
                console.log('cookie removed.')
                $.removeCookie("currentUploadSessionID")
                $.removeCookie("currentUploadProgress")
            	window.setTimeout(function() {
                    alert("Hunting completed successfully!")
                    $('#fileUploadTxt').val('')
                    $('.upload-section').removeClass("hidden")
                    $('.progress-section').addClass("hidden")
                    window.location.reload()
                }, 2000)
            }
        })
    }
    
    uploader.finishCallback = function() {
        setProgressView(1.0)
    }
    
    uploader.huntingStartCallback = function() {
        $('.upload-section').addClass("hidden")
        $('.progress-section').removeClass("hidden")
        currentProgress = 0
        setProgressView(0)
    }

    uploader.errorCallback = function(statusCode) {
        $.removeCookie("currentUploadSessionID")
        $.removeCookie("currentUploadProgress")
        if (statusCode == 409) {
            var r = confirm("These images have already been processed. Are you sure you want to process the images again?")
            if (r) {
                $('.upload-section').removeClass("hidden")
                $('.progress-section').addClass("hidden")
                uploader.doHunting(true)
            } else {
                $('.upload-section').removeClass("hidden")
                $('.progress-section').addClass("hidden")
            }
        } else {
            alert('Hunting Failed due to the backend issue.')
            $('.upload-section').removeClass("hidden")
            $('.progress-section').addClass("hidden")
        }
    }

    $(document).on('click', '.start-hunt-btn', function() {
        // create the session first
        uploader.doHunting()
    })

    $(document).on('click', '.frame-index-btn', function() {
        var frameIndex = $(this).data('frame-index')

        render.gotoFrame(frameIndex)
    })

    // check if already hunting
    var progressInCookie = $.cookie('currentUploadProgress')
    var sessionIDInCookie = $.cookie('currentUploadSessionID')
    if (progressInCookie && sessionIDInCookie) {
        uploader.huntingStartCallback()
        currentProgress = progressInCookie
        uploader.setSessionID(sessionIDInCookie)
        uploader.setHunting()
        setProgressView(currentProgress)
    }

    var huntingRequest = null

    $(window).bind('beforeunload', function(){
          if (uploader.isHunting() && huntingRequest) {
             huntingRequest.abort()
          }
    });

    var render = new PageRenderer()
    render.latestItem = latestItem
    render.render()

    function HunterUploader() {
        var sessionId = -1
        var hunting = false
        this.progressCallback = null
        this.finishCallback = null
        this.errorCallback = null
        this.huntingStartCallback = null
        this.isHunting = function() {
            return hunting
        }
        this.setHunting = function() {
            hunting = true
        }
        this.setSessionID = function(theSessionID) {
            sessionId = theSessionID
        }
        this.uploadFiles = function(isForced) {
            var that = this
            var formData = new FormData();
            // check the observatory code
            var observatoryCode = $('.observatory-code-input').val()
            if (observatoryCode === null || observatoryCode === "") {
                alert("You must input the Observatory Code")
                return
            }
            formData.append("observatoryCode", observatoryCode)

            var fileUpload = document.getElementById("fileUpload");
            if (fileUpload.files.length != 4) {
                alert("You must select exactly 4 files")
                return
            }

            // check file extension
            for (var i = 0; i < fileUpload.files.length; i++) {
                var uploadFile = fileUpload.files[i]
                var fileName = ''
                if ('name' in uploadFile) {
                    fileName = uploadFile.name
                } else {
                    fileName = uploadFile.fileName
                }
                if (fileName.indexOf(".H", fileName.length - ".H".length) === -1
                    && fileName.indexOf(".arch", fileName.length - ".arch".length) === -1
                    && fileName.indexOf(".fits", fileName.length - ".fits".length) === -1
                    && fileName.indexOf(".fit", fileName.length - ".fit".length) === -1) {
                    alert("Unsupported File type")
                    return
                }
            }
            if (sessionId < 0) {
                alert("No session created")
                return
            }
            for (var i = 0; i < fileUpload.files.length; i++) {
                formData.append("files", fileUpload.files[i]);
            }
            hunting = true
            console.log(that.huntingStartCallback)

            if (that.huntingStartCallback) {
                that.huntingStartCallback()
            }
            huntingRequest = $.ajax({
                global:false,
                url : "/detectAsteroids?id=" + sessionId + '&forced=' + (isForced === true ? "true" : "false"),
                type : "POST",
                data : formData,
                cache : false,
                contentType : false,
                processData : false,
                success : function() {
                    hunting = false
                    console.log("Upload success")
                    if (that.finishCallback) {
                        that.finishCallback()
                    }
                },
                error: function(jqXHR, textStatus, errorThrown) {
                    hunting = false
                    console.log("The status: " + textStatus)
                    if (textStatus == "abort") {
                        console.log("Request aborted..")
                        return
                    }

                    if (that.errorCallback) {
                        that.errorCallback(jqXHR.status)
                    }
                    console.log("Failed to upload the files")
                }

            })
        }

        this.doHunting = function(isForced) {
            var that = this
            $.ajax({
                url: "/detectionSessions",
                type: "POST",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success: function(session) {
                    sessionId = session.id
                    console.log("Session created")
                    that.uploadFiles(isForced === true)
                },
                error: function() {
                     alert("Failed to create the session")
                }
            })
        }

        // create the progress checker
        var that = this
        setInterval(function(){
            if (!hunting) {
                return
            }
            if (sessionId < 0) {
                return
            }
            $.ajax({
                url: "/detectionSessions/" + sessionId + '?t=' + Math.random(),
                type: "GET",
                cache: 'false',
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success: function(session) {
                    if (that.progressCallback) {
                        that.progressCallback(session)
                    }
                },
                error: function() {
                    console.log("Failed to get the hunting progress")
                }
            })
        }, 1000)


    }

    function PageRenderer() {
        this.latestItem = null
        this.searchResult = null
        var that = this
        var currentFrame = -1
        var jcarousel = new JcarouselManager($('.jcarousel'))
        this.render =  function() {
            renderLatest()
            if (gSearchResult != null) {
            	renderSearchResult()
            }
        }
        this.gotoFrame = function(index) {
            that.showFrame(index)
        }

        this.showFrame = function(frameIndex) {
        	$('.frame-control a').each(function() {
                $(this).removeClass("selected")
                var index = $(this).data('frame-index')
                if (index == frameIndex) {
                    $(this).addClass("selected")
                }
            })
            var detectItem = that.latestItem
            if (detectItem == null) {
            	return;
            }
            currentFrame = frameIndex
            var frame = this.latestItem.frames[frameIndex]
            $(".result-info").html('')
            var items = []
            items.push({name: "NEO?", value: detectItem.neo ? "Yes" : "No"})
            items.push({name: "RIGHT ACSENSION", value: formatHMS(frame.rightAscension)})
            items.push({name: "DECLINATION", value: formatDMS(frame.declination)})
            items.push({name: "ROUGH MAGNITUDE", value: frame.roughMagnitude.toFixed(2)})
            items.push({name: "IMAGE CREATE DATE / TIME", value: formatDate(detectItem.timestamp)})
            items.push({name: "OBSERVATORY CODE", value: detectItem.observatoryCode})

            items.push({name: "EXISTS AT MPC?", value: detectItem.submitted ? "Yes" : "No Match"})

            for (var i = 0; i < items.length; i++) {
                var item = items[i]
                var itemContent = $(".info-item-template").clone()
                itemContent.removeClass("hidden")
                itemContent.removeClass("info-item-template")
                itemContent.addClass("info-item")

                itemContent.find(".info-item-name").text(item.name)
                itemContent.find(".info-item-value").text(item.value)
                if (i == items.length - 1) {
                    itemContent.find(".seperator").addClass('hidden')
                    var textColorClass = detectItem.submitted ? "green-text" : "red-text"
                    itemContent.find(".info-item-value").removeClass("blue-text").addClass(textColorClass)
                }
                $('.result-info').append(itemContent)
            }
            // show the red circle
            var imageOriginalWidth = detectItem.imageWidth;
            var imageOriginalHeight = detectItem.imageHeight;
            var imageShowWidth = $('.result-img .bg-img').width();
            var imageShowHeight = $('.result-img .bg-img').height();
            var x = frame.x;
            var y = frame.y;
            var circleWidth = 24;
            var circleHeight = 24;
            var top = y * imageShowHeight / imageOriginalHeight - circleWidth / 2;
            var left = x * imageShowWidth / imageOriginalWidth - circleHeight / 2;
            $('.red-circle').css({top: top + "px", left:left + "px"})
            var elevateZoom = $('.result-img .bg-img').data('elevateZoom')
            if (elevateZoom) {
                console.log("here")
                elevateZoom.options.circle = {
                    bgWidth : imageOriginalWidth,
                    bgHeight: imageOriginalHeight,
                    widthPercent : x / imageOriginalWidth,
                    heightPercent: y / imageOriginalHeight
                }
                elevateZoom.resetRedCircle()
            }

        }
        function renderLatest() {
             if (that.latestItem == null) {
                 // show no recent
                 $(".latest-result-no-content").removeClass("hidden")
                 $(".latest-result-content").addClass("hidden")
                 $(".constrast-slider").addClass("hidden")
                 return
             }
            $(".latest-result-content").removeClass("hidden")
            $(".latest-result-no-content").addClass("hidden")
            $(".constrast-slider").removeClass("hidden")

            $('.result-img .bg-img').each(function() {
                removeElevateZoom($(this))
            })

            addElevateZoom($(".result-img .bg-img"))
            that.showFrame(0)

        }

        function renderSearchResult() {
        	
            var searchHistoryTable = new AsteroidTable($('.asteroid-table-container'))
            searchHistoryTable.renderCell = function(item, columnIndex, cell) {
                var text = ""
                if (columnIndex == 0) {
                    text = formatDate(item.timestamp)
                } else if (columnIndex == 2) {
                    text = formatHMS(item.frames[0].rightAscension) + '\n' + formatDMS(item.frames[0].declination)
                } else if (columnIndex == 1) {
                    text =  item.observatoryCode
                }
                cell.text(text)
                cell.closest('tr').addClass('clickable')
                cell.closest('tr').data('item-id', item.id)
                cell.closest('tr').data('detection-item', item)
            }
            searchHistoryTable.renderTable(gSearchResult.values)
        }

        $(document).on('click', '.asteroid-table-container tbody tr', function() {
            window.location.href = '/detectionItems/' + $(this).data('item-id')
        })

        $(document).on('mouseover', '.asteroid-table-container tbody tr', function() {
            that.latestItem = $(this).data('detection-item')
            $('.result-img .bg-img').each(function() {
                if (that.latestItem.imageId != $(this).data('imageId')) {
                    // update the image
                    $(this).attr('src', '/frameimage/' + that.latestItem.id + '/frame/' + 0 + '/image')
                    $(this).data('imageId', that.latestItem.imageId)
                }
            })
            renderLatest()
        })

        var mouseEvent = null
        $(document).on('mousemove', function(event) {
            if (!event.pageX) {
                return
            }
            mouseEvent = $.extend($.Event(event.type /* click, mousedown, touchstart, ect. ect. */), {
                which: 1,
                clientX: event.clientX,
                clientY: event.clientY,
                pageX: event.pageX,
                pageY: event.pageY,
                screenX: event.screenX,
                screenY: event.screenY
            });
        })

        jcarousel.init({
            animation: {
                duration: 0 // make changing image immediately
            }
        })
        jcarousel.onTargetIn = function(index) {
            $('.jcarousel li img').each(function() {
                removeElevateZoom($(this))
            })
            that.showFrame(index)
            var item = $(".jcarousel li img:eq(" + index + ")")
            addElevateZoom(item)
            
            setTimeout(function() {
                $('.zoomContainer').trigger(mouseEvent)
                $('.zoomLens').trigger(mouseEvent)
            }, 100)

            // change the slider
           var contrastValue = item.data('contrastValue')
           if (!contrastValue) {
               contrastValue = 0
           }

            ignoreSlideChange = true
           $('.slider').slider('value', contrastValue);
            updateInvertButtonIcon()
        }

        function ImagePlayer() {
            var index = 0;
            var playInterval = null;
            this.autoPlay = function() {
                if (playInterval != null) {
                    return;
                }
                playInterval = setInterval(function() {
                    index++;
                    if (index > 3) {
                        index = 0;
                    }
                    that.showFrame(index);
                }, 2000)
            }
            this.pausePlay = function() {
                if (playInterval != null) {
                    clearInterval(playInterval)
                    playInterval = null
                }
            }
        }
        var player = new ImagePlayer()
        $(document).on('click', '.play-btn', function() {
            $(this).addClass("hidden")
            $(".pause-btn").removeClass("hidden")
            player.autoPlay()
        })

        $(document).on('click', '.pause-btn', function() {
            $(this).addClass("hidden")
            $(".play-btn").removeClass("hidden")
            player.pausePlay()
        })

        // effects

        $(document).on('click', '.invert-btn', function() {
            var item = $('.img-container .bg-img')
            var vjsAPI = item.data('vintageJS');
            if (item.hasClass('effected-on')) {
                item.removeClass('effected-on')


                vjsAPI.reset()
                updateInvertButtonIcon()
                // recover
                resetElevateZoom()
                return
            }
            item.addClass('effected-on')
            // Invert color
            var invertCurves = {
                r : [],
                g : [],
                b : []
            };

            for (var i = 0; i < 256; i++) {
                invertCurves.r.push(255-i);
                invertCurves.g.push(255-i);
                invertCurves.b.push(255-i);
            }

            vjsAPI.vintage({
                curves : invertCurves
            });

            updateInvertButtonIcon()
            resetElevateZoom()
        })


       $(document).on('click', '.flicker-btn', function() {
           jcarousel.showFlickerEffect()
        })


        function addElevateZoom(element) {
            removeElevateZoom(element)
            element.elevateZoom({
                zoomType	: "lens", lensShape : "round", lensSize : 200, scrollZoom : true
            });
        }
        function removeElevateZoom(element) {
            element.removeData('elevateZoom');
            $('.zoomContainer').remove();
        }

        function resetElevateZoom() {

            removeElevateZoom($(".img-container .bg-img"))
            addElevateZoom($(".img-container .bg-img"))

        }

        $(document).on('RESET_ELEVATE_EVENT', function() {
            resetElevateZoom()
        })

        function updateInvertButtonIcon() {
            var item = $('.img-container .bg-img')
            if (item.hasClass('effected-on')) {
                // hide the slider bar
                $('.constrast-slider').addClass("hidden")
                $('.invert-btn img').attr('src', '../i/invert-btn-on.png')
            } else {
                $('.constrast-slider').removeClass("hidden")
                $('.invert-btn img').attr('src', '../i/invert-btn.png')
            }
        }
    }

    function JcarouselManager(element) {

        var jcarouselElement = element
        var that = this
        var playInterval = null
        var currentIndex = 0
        this.onTargetIn = function(index) {
            // should be override
        }

        this.init = function(config) {
            jcarouselElement.jcarousel(config);
        }

        this.showFlickerEffect = function() {
            if (playInterval != null) {
                // already playing
                return
            }
            jcarouselElement.jcarousel('scroll', "0");
            playInterval = setInterval(function(){
                if (currentIndex == 3) {
                    // the last one
                    jcarouselElement.jcarousel('scroll', "0");
                    clearInterval(playInterval)
                    playInterval = null
                } else {
                    jcarouselElement.jcarousel('scroll', "+=1");
                }
            }, 1000)
        }

        this.autoPlay = function() {
            if (playInterval != null) {
                return
            }

            playInterval = setInterval(function(){
                if (currentIndex == 3) {
                    // the last one
                    jcarouselElement.jcarousel('scroll', "0");
                } else {
                    jcarouselElement.jcarousel('scroll', "+=1");
                }
            }, 1000)
        }

        this.pausePlay = function() {
            if (playInterval != null) {
                clearInterval(playInterval)
                playInterval = null
            }
        }
        this.gotoFrame = function (index) {
            jcarouselElement.jcarousel('scroll', index + "");
        }
        jcarouselElement.on('jcarousel:targetin', "li", function() {
            var index = jcarouselElement.find("li").index($(this))
            currentIndex = index
            that.onTargetIn(index)
        })
    }


    var sendGmail = function(opts){
        var str = 'http://mail.google.com/mail/?view=cm&fs=1'+
            '&to=' + opts.to +
            '&su=' + opts.subject +
            '&body=' + opts.message +
            '&ui=1';
        window.open(str, '_blank');
    }

    $('.send-email-btn a').on('click', function() {
        sendGmail({to: observationsMailAddress, subject: observationsMailSubject, message: observationsMailBody})
    })
})
