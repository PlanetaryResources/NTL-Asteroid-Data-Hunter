/**
 * The javascript using in the Topic page.
 * @author TCSASSEMBLER
 * @version 1.0
 */
$(document).ready(function() {
    showTopics()
    var lister = null

    if (gViewName == "listHelpItems") {
        showTopicItems(gTopics[0].id)
    } else if (gViewName == "viewHelpItem") {
         showHelpItem(gHelpItem)
    }



    var itemsTable = new AsteroidTable($('.item-list'))
    itemsTable.changePageCallBack = function(pageIndex) {
        lister.pageNumber = pageIndex
        lister.search()
    }
    itemsTable.renderCell = function(item, columnIndex, cell) {
        cell.text(item.title)
        cell.addClass("left-align")
        cell.closest('tr').addClass('clickable')
        cell.closest('tr').data('topic-item-id', item.id)
    }

    function TopicItemsLister(criteria) {
        this.successCallback = null
        this.pageSize = 9
        this.pageNumber = 1
        this.sortBy = 'id'
        this.sortType = 'ASC'
        var that = this
        this.search = function() {
            var options = $.extend(criteria, {
                pageSize: that.pageSize,
                pageNumber: that.pageNumber,
                sortBy: that.sortBy,
                sortType: that.sortType})
            $.ajax({
                url:'/search/helpItems?' + $.param(options),
                type: "GET",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success:function(result) {
                    if (that.successCallback) {
                        that.successCallback(result)
                    }
                    itemsTable.renderTable(result.values)
                    itemsTable.renderPaging(result.totalPages, result.pageNumber)
                },
                error: function() {
                	// we do not alert this, because this error will promp if we manually abort the requests.
                	console.log("Error failed in request the topics.");
                }
            })
        }
    }

    function showTopicItems(topicId) {
        lister = new TopicItemsLister({
            topicId: topicId
        })
        lister.successCallback = function(result) {
            // choose the selected topic
           highlightTopic(topicId)

            $('.content-panel').each(function() {
                $(this).addClass("hidden")
            })
            $(".items-panel").removeClass("hidden")

            $('.items-panel .item-topic-title').text(getTopic(topicId).name)
        }
        lister.search()
    }

    function highlightTopic(topicId) {
        // choose the selected topic
        $('.topics-list .topic-item').each(function() {
            $(this).removeClass("selected")
            if ($(this).data('topic-id') == topicId) {
                $(this).addClass("selected")
            }
        })
    }

    function getTopic(topicId) {
        for (var i = 0; i < gTopics.length; i++) {
            if (gTopics[i].id == topicId) {
                return gTopics[i]
            }
        }
        return null
    }

    function showHelpItem(item) {
        $('.content-panel').each(function() {
            $(this).addClass('hidden')
        })
        $('.item-content-panel').removeClass("hidden")
        $('.item-content-panel .item-content').html(item.content)
        $('.item-content-panel .item-topic-title').text(item.title)

        // highlight the item's topic
        highlightTopic(item.topic.id)

    }

    function showTopics() {
       var list = $('.topics-list')
        list.html('')
        for (var i = 0; i < gTopics.length; i++) {
            var topic = gTopics[i];
            var topicItem = $('<div></div>').addClass("topic-item").data('topic-id', topic.id)
            var topicContent = $('<span></span>').addClass("topic-name").text(topic.name)
            topicItem.append(topicContent)
            list.append(topicItem)
            var seperator = $('<div></div>').addClass('seperator')
            list.append(seperator)
        }
    }

    // for search
    var helpItemSearcher = null;
    var searchTable = new AsteroidTable($('.items-result-table'))
    searchTable.changePageCallBack = function(pageIndex) {
        helpItemSearcher.pageNumber = pageIndex
        helpItemSearcher.search()
    }
    searchTable.renderCell = function(item, columnIndex, cell) {
        var text = "";
        if (columnIndex == 0) {
            text = item.topic.name
        } else if (columnIndex == 1) {
            text = item.title
        }
        cell.text(text)
        cell.closest('tr').addClass('clickable')
        cell.closest('tr').data('topic-item-id', item.id)
    }

    function HelpItemSearcher(criteria) {
        this.successCallback = null
        this.pageSize = 9
        this.pageNumber = 1
        this.sortBy = 'id'
        this.sortType = 'ASC'
        var that = this
        this.search = function() {
            var options = $.extend(criteria, {
                pageSize: that.pageSize,
                pageNumber: that.pageNumber,
                sortBy: that.sortBy,
                sortType: that.sortType})
            $.ajax({
                url:'/search/helpItems?' + $.param(options),
                type: "GET",
                headers: {
                    'Accept': 'application/json'
                },
                dataType : 'json',
                success:function(result) {
                    if (that.successCallback) {
                        that.successCallback(result)
                    }
                    searchTable.renderTable(result.values)
                    searchTable.renderPaging(result.totalPages, result.pageNumber)
                },
                error: function() {
                	// we do not alert this, because this error will promp if we manually abort the requests.
                	console.log("Error failed in request the help items.");
                }
            })
        }
    }

    $(document).on('click', '.topic-item', function() {
        var topicId = $(this).data('topic-id')
        showTopicItems(topicId)
    })

    $(document).on('click', '.item-list tbody tr', function() {
        window.location.href = '/helpItems/' + $(this).data('topic-item-id')
    })
    $(document).on('click', '.search-result-panel tbody tr', function() {
        window.location.href = '/helpItems/' + $(this).data('topic-item-id')
    })

    $(document).on('click', '.help-search-btn', function() {
        var term = $(this).closest('.search-box').find('input').val()
        helpItemSearcher = new HelpItemSearcher({
            keyword : term
        })
        helpItemSearcher.successCallback = function(result) {
            $('.content-panel').each(function() {
                $(this).addClass("hidden")
            })
            
            $('.search-term-title').text(term)
            if (!term) {
            	$('.search-term-title').parent().addClass("hidden")
            } else {
            	$('.search-term-title').parent().removeClass("hidden")
            }
            $(".search-result-panel").removeClass("hidden")
            highlightTopic(-1)
        }
        helpItemSearcher.search()

    })
    
    $('.search-box').on("keypress", function(e) {
        /* ENTER PRESSED*/
        if (e.keyCode == 13) {
           $('.help-search-btn').trigger('click');
            return false;
        }
    })
})