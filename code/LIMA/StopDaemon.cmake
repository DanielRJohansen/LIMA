# Stops the daemon of the lima executable about to be replaced. Lets its running commands finish first.
if(EXISTS "${LIMA_EXE}")
    execute_process(COMMAND "${LIMA_EXE}" daemon stop OUTPUT_QUIET ERROR_QUIET)
endif()
