actionButton3 <- function (inputId, label, icon = NULL, cl="btn btn-default action-button",...) {
  tags$button(id = inputId, type = "button", class = cl, list(icon, label), ...)
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
         if($("input:radio[name=input_type]:checked").val() == "pdb"){
             $("#pop-pdbid").popover("enable");
         } else {
             console.log("conditions for sequence and multipdb will be added")
             $("#pop-pdbid").popover("disable");
         }
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
                     icon = 'question', content = NULL,
                     cl = "btn btn-primary btn-input action-button",
                     ...) {

  idm1 <- paste0("myModal", id)
  idm2 <- paste0("#", idm1)

  if(!is.null(icon))
    icon <- icon(icon)

  tags$div(
    tags$button(type="button", class=cl,
                "data-toggle"="modal", "data-target"=idm2,
                style = "float: right",
                list(icon, button_label)
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

