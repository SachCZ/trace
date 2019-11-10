#Disable messages from extrenaly build projects
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif ()
endfunction()