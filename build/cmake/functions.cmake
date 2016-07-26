#/////////////////////// Utility functions to be used by CMake project

function(ncbi_copy_exec_to_old_outdir target_name)
    if (NGS_TOOLS_OUTDIR_ENABLED)
        if (UNIX)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory ${NGS_TOOLS_OUTDIR_PREFIX}/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/bin/
                # copy both TARGET_FILE (which is target executable) $NGS_TOOLS_OUTDIR
                COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${target_name}> ${NGS_TOOLS_OUTDIR_PREFIX}/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/bin/$<TARGET_FILE_NAME:${target_name}>
            )
        elseif (WIN32)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory ${NGS_TOOLS_OUTDIR_PREFIX}/${OS}/${PLATFORM_TOOLSET}/${WIN_PLATFORM}/$<CONFIGURATION>/bin/
                # copy both TARGET_FILE (which is exe file) and ${target_name}.pdb to $NGS_TOOLS_OUTDIR (pdb file name is specified via "${target_name}.pdb" since CMake 2.8.12 does not support variable TARGET_PDB_FILE)
                COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${target_name}> $<$<CONFIG:Debug>:$<TARGET_FILE_DIR:${target_name}>/${target_name}.pdb> ${NGS_TOOLS_OUTDIR_PREFIX}/${OS}/${PLATFORM_TOOLSET}/${WIN_PLATFORM}/$<CONFIGURATION>/bin/
            )
        endif()
    endif()
endfunction()
