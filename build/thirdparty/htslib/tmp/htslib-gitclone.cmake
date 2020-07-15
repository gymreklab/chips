if("1.3.1" STREQUAL "")
  message(FATAL_ERROR "Tag for git checkout should not be empty.")
endif()

set(run 0)

if("/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-gitinfo.txt" IS_NEWER_THAN "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-gitclone-lastrun.txt")
  set(run 1)
endif()

if(NOT run)
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib'")
endif()

set(git_options)

# disable cert checking if explicitly told not to do it
set(tls_verify "")
if(NOT "x" STREQUAL "x" AND NOT tls_verify)
  list(APPEND git_options
    -c http.sslVerify=false)
endif()

set(git_clone_options)

set(git_shallow "")
if(git_shallow)
  list(APPEND git_clone_options --depth 1 --no-single-branch)
endif()

# try the clone 3 times incase there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git" ${git_options} clone ${git_clone_options} --origin "origin" "https://github.com/samtools/htslib.git" "htslib"
    WORKING_DIRECTORY "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/samtools/htslib.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git" ${git_options} checkout 1.3.1
  WORKING_DIRECTORY "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '1.3.1'")
endif()

execute_process(
  COMMAND "/usr/bin/git" ${git_options} submodule init 
  WORKING_DIRECTORY "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to init submodules in: '/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib'")
endif()

execute_process(
  COMMAND "/usr/bin/git" ${git_options} submodule update --recursive --init 
  WORKING_DIRECTORY "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-gitinfo.txt"
    "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-gitclone-lastrun.txt"
  WORKING_DIRECTORY "/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/storage/mlamkin/projects/chips/build/thirdparty/htslib/src/htslib-stamp/htslib-gitclone-lastrun.txt'")
endif()

