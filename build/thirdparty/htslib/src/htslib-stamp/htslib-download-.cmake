

set(command "/usr/local/bin/cmake;-P;/storage/mlamkin/projects/chips/build/thirdparty/htslib/tmp/htslib-gitclone.cmake")
execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_FILE "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-download-out.log"
  ERROR_FILE "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-download-err.log"
  )
if(result)
  set(msg "Command failed: ${result}\n")
  foreach(arg IN LISTS command)
    set(msg "${msg} '${arg}'")
  endforeach()
  set(msg "${msg}\nSee also\n  /storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-download-*.log")
  message(FATAL_ERROR "${msg}")
else()
  set(msg "htslib download command succeeded.  See also /storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-download-*.log")
  message(STATUS "${msg}")
endif()
