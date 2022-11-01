# 输入变量:
# 输入XXX_PATH, THIRDLIB_PATH, 
# 依据这两个输入路径搜索需要的头文件和库文件,
#
# 输出变量:
# 输出XXX_FOUND, silo_INCLUDES, silo_LIBRARIES.
# 并将XXX_INCLUDES追加到THIRDLIB_INCLUDES,
# 将XXX_LIBRARIES追加到THIRDLIB_LIBRARIES


# 定义搜索路径
set(silo_search_PATHS ${silo_PATH})
foreach (path ${THIRDLIB_PATH})
  list(APPEND silo_search_PATHS ${path}/silo)
endforeach()

# 这个模块由很多的库文件组成
# silo with HDF5 driver
if (NOT silo_FIND_COMPONENTS)
  set(silo_FIND_COMPONENTS siloh5)
endif ()

# 设置搜索路径
set(_OLD_PREFIX_PATH ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${silo_search_PATHS})
if (NOT THIRDLIB_FIND_ARGUMENTS)
  set(THIRDLIB_FIND_ARGUMENTS "PATHS ${silo_search_PATHS} PATH_SUFFIXES lib NO_DEFAULT_PATH")
endif ()

# 搜索需要的头文件
find_file(silo_INCLUDES NAMES include ${THIRDLIB_FIND_ARGUMENTS})
if (silo_INCLUDES)
  # do nothing
elseif (silo_FIND_REQUIRED)
  message(FATAL_ERROR "Not found 'include' directory in ${silo_search_PATHS}")
endif ()

# 搜索需要的库文件
foreach (module ${silo_FIND_COMPONENTS})
  find_library(silo_LIBRARIES_${module} NAMES ${module} ${THIRDLIB_FIND_ARGUMENTS})
  if (silo_LIBRARIES_${module})
    list(APPEND silo_LIBRARIES ${silo_LIBRARIES_${module}})
  else ()
    list(APPEND silo_NO_LIBRARIES ${silo_LIBRARIES_${module}})
    if (silo_FIND_REQUIRED)
      message(FATAL_ERROR "Not found lib/'${module}' in ${silo_search_PATHS}")
    endif ()
  endif ()
endforeach ()

# 恢复CMAKE_PREFIX_PATH
set(CMAKE_PREFIX_PATH ${_OLD_PREFIX_PATH})

# 输出信息
if (NOT silo_FIND_QUIETLY)
  message(STATUS "Found silo: ${silo_INCLUDES}, ${silo_LIBRARIES};${silo_NO_LIBRARIES} in ${silo_search_PATHS}")
endif ()

# 定义输出变量
if (silo_INCLUDES AND silo_LIBRARIES)
  set(silo_FOUND TRUE)
  string(TOUPPER silo _UPPER_CASE)
  set(HAVE_${_UPPER_CASE} 1)
  unset(_UPPER_CASE)
else ()
  set(silo_FOUND FALSE)
endif ()

# 只是为了方便, 其实破坏了封装
# 追加包含变量到THIRDLIB变量里面
if (silo_INCLUDES)
  list(APPEND THIRDLIB_INCLUDES ${silo_INCLUDES})
endif ()
# 追加库文件变量到THIRDLIB变量里面, 只是为了方便, 其实破坏了封装
if (silo_LIBRARIES)
  list(APPEND THIRDLIB_LIBRARIES ${silo_LIBRARIES})
endif ()

