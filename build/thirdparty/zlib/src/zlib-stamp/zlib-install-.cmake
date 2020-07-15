

set(command "${make};install")
execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_FILE "/storage/mlamkin/projects/chips/build/thirdparty/zlib/src/zlib-stamp/zlib-install-out.log"
  ERROR_FILE "/storage/mlamkin/projects/chips/build/thirdparty/zlib/src/zlib-stamp/zlib-install-err.log"
  )
if(result)
  set(msg "Command failed: ${result}\n")
  foreach(arg IN LISTS command)
    set(msg "${msg} '${arg}'")
  endforeach()
  set(msg "${msg}\nSee also\n  /storage/mlamkin/projects/chips/build/thirdparty/zlib/src/zlib-stamp/zlib-install-*.log")
  message(FATAL_ERROR "${msg}")
else()
  set(msg "zlib install command succeeded.  See also /storage/mlamkin/projects/chips/build/thirdparty/zlib/src/zlib-stamp/zlib-install-*.log")
  message(STATUS "${msg}")
endif()
