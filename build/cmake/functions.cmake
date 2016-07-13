#/////////////////////// Utility functions to be used by CMake project

function(ncbi_copy_exec_to_old_outdir target_name)
    if (NGS_TOOLS_OUTDIR_ENABLED)
        if (UNIX)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${target_name}> ${NGS_TOOLS_OUTDIR_PREFIX}/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/bin/$<TARGET_FILE_NAME:${target_name}>
            )
        elseif (WIN32)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${target_name}> ${NGS_TOOLS_OUTDIR_PREFIX}/${OS}/${PLATFORM_TOOLSET}/${WIN_PLATFORM}/$<CONFIGURATION>/bin/$<TARGET_FILE_NAME:${target_name}>
            )
        endif()
    endif()
endfunction()
