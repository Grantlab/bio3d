actionButton3 <- function (inputId, label, icon = NULL, cl="btn btn-default action-button",...) {
  tags$button(id = inputId, type = "button", class = cl, list(icon, label), ...)
}

dialogBox <- function(id = 1, icon = "question-circle", title = NULL, style = NULL, ... ) {

  ## CSS for this stuff is in www/css/jquery-ui.css
  ido <- paste0("opener", id)
  idd <- paste0("dialog", id)

  script <- sprintf(
    '$(function() {$( "#%s" ).dialog({autoOpen: false, width: 600, modal: true, show: {effect: "fade", duration: 500}}); $( "#%s" ).click(function() { $( "#%s" ).dialog( "open" ); });  });',
    idd, ido, idd
    )

  tags$div(
    tags$div(id=ido, icon(icon), style = style),
    tags$div(id=idd, title = title, ...),

    tags$script(
      script
      )
    )
}

popoverQuestion <- function(id = 'popQues1', icon = 'question-circle', title = 'Title',
                      content = '', trigger = 'hover', style = "position: absolute; right: 25px; top: 5px;",
                      placement = 'right', data_toggle = "pop-box", ...) {
    script <- sprintf(
        'jQuery(function(){
         $("[data-toggle=\'%s\']")
         .popover({html: true,
         container: "body",
         content: "%s",
         delay: { "show": 50, "hide": 20 },
         placement: "auto %s"
         });
         $("#%s").click(function(e){
         $("[data-toggle=\'%s\']").popover("show");
         });
         $(document).click(function (e) {
         if (($(".popover").has(e.target).length == 0) && $("#%s").has(e.target).length == 0) {
         $("[data-toggle=\'%s\']").popover("hide");
         }
         });
         });
         ', data_toggle, content, placement, id, data_toggle, id, data_toggle
    )

    tags$div(
        tags$div(id=id, icon(icon), style = style
#                 ,"data-toggle" = "pop-box",
#                 "tabindex" = 0, "data-trigger" = "manual"
                 )
        ,tags$script(HTML( script ))
    )

}

popRadioButtons <- function(inputId, label, choices, selected = NULL, inline = FALSE,
                            title = 'Title', content = '' , placement = 'right',
                            data_toggle = 'pop-box', ...) {
        tags$div(id = paste0("pop-", inputId), "data-toggle" = data_toggle,
            "data-trigger" = "manual", "data-original-title" = title, ...,
            radioButtons(inputId, label, choices, selected, inline)
        )
}

popTextInput <- function(inputId, label, value = "",
                         title = 'Title', content = '', placement = 'right',
                         data_toggle = "pop-box", ...) {
    tags$div(id = paste0('pop-', inputId), "data-toggle" = data_toggle,
             "data-trigger" = "manual", "data-original-title" = title,
             "data-content" = content, ...,
             textInput(inputId, label, value)
    )
}

modalBox <- function(id, button_label = "", heading = NULL,
                     icon = NULL, content = NULL,
                     cl = "btn btn-primary btn-input action-button",
                     ...) {

  idm1 <- paste0("myModal", id)
  idm2 <- paste0("#", idm1)

  tags$div(
    tags$button(type="button", class=cl,
                "data-toggle"="modal", "data-target"=idm2,
                style = "float: right",
                list(icon("question"), button_label)
                ),

  tags$div(
    class="modal fade", id=idm1, tabindex="-1", role="dialog",
    "aria-labelledby"="myModalLabel",


    tags$div(
      class="modal-dialog", role="document",

      tags$div(
        class="modal-content",

        tags$div(class="modal-header",
                 tags$button(
                   type="button", class="close",
                   "data-dismiss"="modal",
                   "aria-label"="Close",
                   tags$span("aria-hidden"="true", "Close")
                   ),
                 h4(heading, class="modal-title", id="myModalLabel")
                 ),

        tags$div(
          class="modal-body",
          content
          )
        )
      )
    )
    )
}

