# 输入变量:
# 输入zlib_PATH, THIRDLIB_PATH, 
# 依据这两个输入路径搜索需要的头文件和库文件,
#
# 输出变量:
# 输出zlib_FOUND, zlib_INCLUDES, zlib_LIBRARIES.
# 并将zlib_INCLUDES追加到THIRDLIB_INCLUDES,
# 将zlib_LIBRARIES追加到THIRDLIB_LIBRARIES


# 定义搜索路径
set(zlib_search_PATHS ${zlib_PATH})
foreach (path ${THIRDLIB_PATH})
  list(APPEND zlib_search_PATHS ${path}/zlib)
endforeach()

# 这个模块由很多的库文件组成
if (NOT zlib_FIND_COMPONENTS)
  set(zlib_FIND_COMPONENTS z)
endif ()

# 设置搜索路径
set(_OLD_PREFIX_PATH ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${zlib_search_PATHS})
if (NOT THIRDLIB_FIND_ARGUMENTS)
  set(THIRDLIB_FIND_ARGUMENTS "PATHS ${zlib_search_PATHS} PATH_SUFFIXES lib NO_DEFAULT_PATH")
endif ()

# 搜索需要的头文件
find_file(zlib_INCLUDES NAMES include ${THIRDLIB_FIND_ARGUMENTS})
if (zlib_INCLUDES)
  # do nothing
elseif (zlib_FIND_REQUIRED)
  message(FATAL_ERROR "Not found 'include' directory in ${zlib_search_PATHS}")
endif ()

# 搜索需要的库文件
foreach (module ${zlib_FIND_COMPONENTS})
  find_library(zlib_LIBRARIES_${module} NAMES ${module} ${THIRDLIB_FIND_ARGUMENTS})
  if (zlib_LIBRARIES_${module})
    list(APPEND zlib_LIBRARIES ${zlib_LIBRARIES_${module}})
  else ()
    list(APPEND zlib_NO_LIBRARIES ${zlib_LIBRARIES_${module}})
    if (zlib_FIND_REQUIRED)
      message(FATAL_ERROR "Not found lib/'${module}' in ${zlib_search_PATHS}")
    endif ()
  endif ()
endforeach ()

# 恢复CMAKE_PREFIX_PATH
set(CMAKE_PREFIX_PATH ${_OLD_PREFIX_PATH})

# 输出信息
if (NOT zlib_FIND_QUIETLY)
  message(STATUS "Found zlib: ${zlib_INCLUDES}, ${zlib_LIBRARIES};${zlib_NO_LIBRARIES} in ${zlib_search_PATHS}")
endif ()

# 定义输出变量
if (zlib_INCLUDES AND zlib_LIBRARIES)
  set(zlib_FOUND TRUE)
  string(TOUPPER zlib _UPPER_CASE)
  set(HAVE_${_UPPER_CASE} 1)
  unset(_UPPER_CASE)
else ()
  set(zlib_FOUND FALSE)
endif ()

# 只是为了方便, 其实破坏了封装
# 追加包含变量到THIRDLIB变量里面
if (zlib_INCLUDES)
  list(APPEND THIRDLIB_INCLUDES ${zlib_INCLUDES})
endif ()
# 追加库文件变量到THIRDLIB变量里面, 只是为了方便, 其实破坏了封装
if (zlib_LIBRARIES)
  list(APPEND THIRDLIB_LIBRARIES ${zlib_LIBRARIES})
endif ()

