#gcta CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(GCTA C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)


#set project source directory
set(PRJ_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CUSTORM_INCLUDE_SERACH_PATH FALSE CACHE PATH "")
set(CUSTORM_LIB_SEARCH_PATH FALSE CACHE PATH "")
set(BUILD_STATIC FALSE CACHE BOOL "")


#information of CPU architecture
message(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR: " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message("")


#information of operation system
message(STATUS ${CMAKE_HOST_SYSTEM})
if (UNIX)
    set(system "UNIX_like")
elseif(APPLE)
    set(system "macOS")
elseif(WIN32)
    set(system "windows")
else()
    message(FATAL_ERROR "Operation system not recognized, only support unix apple and windows.")
endif()
message(STATUS "OS: " ${system})
message("")


#information of compilation toolchain
message(STATUS "CMAKE_C_COMPILER:  " ${CMAKE_C_COMPILER} " " ${CMAKE_C_COMPILER_ID} ":" ${CMAKE_C_COMPILER_VERSION} " " ${CMAKE_C_STANDARD})
message(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} " " ${CMAKE_CXX_COMPILER_ID} ":" ${CMAKE_CXX_COMPILER_VERSION} " " ${CMAKE_CXX_STANDARD})
message(STATUS "CMAKE_AR: " ${CMAKE_AR})
message(STATUS "CMAKE_RANLIB: " ${CMAKE_RANLIB})
message(STATUS "CMAKE_COMMAND: " ${CMAKE_COMMAND})
message(STATUS "CMAKE_VERSION: " ${CMAKE_VERSION})
message(STATUS "CMAKE_GENERATOR: " ${CMAKE_GENERATOR})
message(STATUS "CMAKE_BUILD_TOOL: " ${CMAKE_BUILD_TOOL})

message(STATUS "CMAKE_SHARED_LIBRARY_PREFIX: " ${CMAKE_SHARED_LIBRARY_PREFIX})
message(STATUS "CMAKE_SHARED_LIBRARY_SUFFIX: " ${CMAKE_SHARED_LIBRARY_SUFFIX})
message(STATUS "CMAKE_STATIC_LIBRARY_PREFIX: " ${CMAKE_STATIC_LIBRARY_PREFIX})
message(STATUS "CMAKE_STATIC_LIBRARY_SUFFIX: " ${CMAKE_STATIC_LIBRARY_SUFFIX})

message(STATUS "CMAKE_FIND_LIBRARY_PREFIXES: " ${CMAKE_FIND_LIBRARY_PREFIXES})
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
message(STATUS "CMAKE_FIND_LIBRARY_SUFFIXES: " ${CMAKE_FIND_LIBRARY_SUFFIXES})

message(STATUS "CMAKE_PROJECT_NAME: " ${CMAKE_PROJECT_NAME})
message(STATUS "CMAKE_PROJECT_VERSION: " ${CMAKE_PROJECT_VERSION})
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message(STATUS "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message("")


#detect required packages
function(FIND_PKG pkg_path include_path lib_path include_dir bool_include lib_dir bool_lib head_name lib_name)
    if (pkg_path)
        if (bool_include)
            set(${include_dir} ${pkg_path}/include PARENT_SCOPE)
        endif ()
        if (bool_lib)
            set(${lib_dir} ${pkg_path}/lib PARENT_SCOPE)
        endif ()
    elseif (include_path OR lib_path)
        if (include_path AND bool_include)
            set(${include_dir} ${include_path} PARENT_SCOPE)
        endif()
        if (lib_path AND bool_lib)
            set(${lib_dir} ${lib_path} PARENT_SCOPE)
        endif()
    else ()
        if (bool_include)
            unset(include_path_find CACHE)
            find_path(include_path_find ${head_name} ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES} ${CUSTORM_INCLUDE_SERACH_PATH})
            set(${include_dir} ${include_path_find} PARENT_SCOPE)
        endif ()
        if (bool_lib)
            unset(lib_path_find CACHE)
            find_library(lib_path_find ${lib_name} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES} ${CUSTORM_LIB_SEARCH_PATH})
            get_filename_component(lib_path_find ${lib_path_find} DIRECTORY)
            set(${lib_dir} ${lib_path_find} PARENT_SCOPE)
        endif()
    endif ()
endfunction()

#detect MKL
set(mkl_path FALSE CACHE PATH "")
set(mkl_include_path FALSE CACHE PATH "")
set(mkl_lib_path FALSE CACHE PATH "")
set(mkl_include_dir FALSE)
set(mkl_lib_dir FALSE)
FIND_PKG(${mkl_path} ${mkl_include_path} ${mkl_lib_path} "mkl_include_dir" TRUE "mkl_lib_dir" TRUE "mkl.h" "mkl_core")
if (mkl_include_dir AND mkl_lib_dir)
    message(STATUS "Find mkl head file: " ${mkl_include_dir})
    message(STATUS "Find mkl library: " ${mkl_lib_dir})
else ()
    message(WARNING "Can not find MKL header or library")
endif ()

#detect openblas
set(openblas_path FALSE CACHE PATH "" )
set(openblas_include_path FALSE CACHE PATH "" )
set(openblas_lib_path FALSE CACHE PATH "" )
set(openblas_include_dir FALSE)
set(openblas_lib_dir FALSE)
FIND_PKG(${openblas_path} ${openblas_include_path} ${openblas_lib_path} "openblas_include_dir" TRUE "openblas_lib_dir" TRUE "cblas.h" "openblas")
if (openblas_include_dir AND openblas_lib_dir)
    message(STATUS "Find openblas head file: " ${openblas_include_dir})
    message(STATUS "Find openblas library: " ${openblas_lib_dir})
else ()
    message(STATUS "Can not find openblas header or library.")
endif ()

if ((mkl_include_dir AND mkl_lib_dir) OR (openblas_include_dir AND openblas_lib_dir))
else ()
    message(FATAL_ERROR "Can not found neither MKL or openblas.")
endif ()

#detect boost
set(boost_path FALSE CACHE PATH "")
set(boost_include_path FALSE CACHE PATH "")
set(boost_lib_path FALSE CACHE PATH "")
set(boost_include_dir FALSE)
set(boost_lib_dir FALSE)
FIND_PKG(${boost_path} ${boost_include_path} ${boost_lib_path} "boost_include_dir" TRUE "boost_lib_dir" FALSE "boost" FALSE)
if (boost_include_dir)
    message(STATUS "Find boost head file: " ${boost_include_dir})
else ()
    message(FATAL_ERROR "Can not find boost header.")
endif ()

#detect eigen
set(eigen_path FALSE CACHE PATH "")
set(eigen_include_path FALSE CACHE PATH "")
set(eigen_lib_path FALSE CACHE PATH "")
set(eigen_include_dir FALSE)
set(eigen_lib_dir FALSE)
FIND_PKG(${eigen_path} ${eigen_include_path} ${eigen_lib_path} "eigen_include_dir" TRUE "eigen_lib_dir" FALSE "eigen3" FALSE)
if (eigen_include_dir)
    message(STATUS "Find eigen head file: " ${eigen_include_dir})
else ()
    message(FATAL_ERROR "Can not find eigen.")
endif ()

#detect spectra
set(spectra_path FALSE CACHE PATH "")
set(spectra_include_path FALSE CACHE PATH "")
set(spectra_lib_path FALSE CACHE PATH "")
set(spectra_include_dir FALSE)
set(spectra_lib_dir FALSE)
FIND_PKG(${spectra_path} ${spectra_include_path} ${spectra_lib_path} "spectra_include_dir" TRUE "spectra_lib_dir" FALSE "Spectra" FALSE)
if (spectra_include_dir)
    message(STATUS "Find spectra head file: " ${spectra_include_dir})
else ()
    message(FATAL_ERROR "Can not find spectra head file.") 
endif ()

#detect zlib
set(zlib_path FALSE CACHE PATH "")
set(zlib_include_path FALSE CACHE PATH "")
set(zlib_lib_path FALSE CACHE PATH "")
set(zlib_include_dir FALSE)
set(zlib_lib_dir FALSE)
FIND_PKG(${zlib_path} ${zlib_include_path} ${zlib_lib_path} "zlib_include_dir" TRUE "zlib_lib_dir" TRUE "zlib.h" "z")
if (zlib_include_dir AND zlib_lib_dir)
    message(STATUS "Find zlib head file: " ${zlib_include_dir})
    message(STATUS "Find zlib library: " ${zlib_lib_dir})
else ()
    message(FATAL_ERROR "Can not find zlib header or library.")
endif ()

#detect zstd
set(zstd_path FALSE CACHE PATH "")
set(zstd_include_path FALSE CACHE PATH "")
set(zstd_lib_path FALSE CACHE PATH "")
set(zstd_include_dir FALSE)
set(zstd_lib_dir FALSE)
FIND_PKG(${zstd_path} ${zstd_include_path} ${zstd_lib_path} "zstd_include_dir" TRUE "zstd_lib_dir" TRUE "zstd.h" "zstd")
if (zstd_include_dir AND zstd_lib_dir)
    message(STATUS "Find zstd head file: " ${zstd_include_dir})
    message(STATUS "Find zstd library: " ${zstd_lib_dir})
else ()
    message(FATAL_ERROR "Can not find zstd header or librayr.")
endif ()

#detect gsl
set(gsl_path FALSE CACHE PATH "")
set(gsl_include_path FALSE CACHE PATH "")
set(gsl_lib_path FALSE CACHE PATH "")
set(gsl_include_dir FALSE)
set(gsl_lib_dir FALSE)
FIND_PKG(${gsl_path} ${gsl_include_path} ${gsl_lib_path} "gsl_include_dir" TRUE "gsl_lib_dir" TRUE "gsl" "gsl")
if (gsl_include_dir AND gsl_lib_dir)
    message(STATUS "Find gsl head files: " ${gsl_include_dir})
    message(STATUS "Find gsl library: " ${gsl_lib_dir})
else ()
    message(FATAL_ERROR "Can not find gsl header or library.")
endif ()


#get source files
file(GLOB SRC_FILES_CXX ${PRJ_SOURCE_DIR}/*.cpp)
file(GLOB SRC_FILES_CC ${PRJ_SOURCE_DIR}/*.c)
set(SRC_FILES ${SRC_FILES_CXX} ${SRC_FILES_CC})
message(${SRC_FILES})



#build
add_executable(gcta ${SRC_FILES})
add_subdirectory(submods)


if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
    #choose mkl or openblas
    if (mkl_include_dir AND mkl_lib_dir)
        set(blas_include_dir ${mkl_include_dir})
        set(blas_lib_dir ${mkl_lib_dir})
    else ()
        set(blas_include_dir ${openblas_include_dir})
        set(blas_lib_dir ${openblas_include_dir})
    endif ()

    if (system STREQUAL "UNIX_like")
        if (mkl_include_dir AND mkl_lib_dir)
            set(blas_link_item "mkl_core" "mkl_intel_lp64" "mkl_gnu_thread")
        else ()
            set(blas_link_item "")
        endif ()
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            #target_compile_definitions(gcta)
            target_compile_options(gcta PRIVATE "-fopenmp")
            target_include_directories(gcta PRIVATE ${blas_include_dir} PRIVATE ${boost_include_dir} 
                PRIVATE ${eigen_include_dir} PRIVATE ${spectra_include_dir} PRIVATE ${zlib_include_dir} 
                PRIVATE ${zstd_include_dir} PRIVATE ${gls_include_dir} PRIVATE ${CMAKE_SOURCE_DIR}) 
            target_link_directories(gcta PRIVATE ${blas_lib_dir} PRIVATE ${zlib_lib_dri} PRIVATE ${zstd_lib_dir} PRIVATE ${gsl_lib_dir})
            target_link_options(gcta PRIVATE "-static-libgcc" PRIVATE "-static-libgcc++")
            target_link_libraries(gcta Pgenlib -Wl,--start-group ${blas_link_item} -Wl,--end-group pthread z zstd gsl gslcblas)
        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        else ()
            message(WARNING "Under UNIX_list system, one GCC and Clang was supported.")
        endif ()
    elseif (system STREQUAL "macOS")

    elseif (system STREQUAL "windows")

    else ()
    endif ()
elseif (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")

else ()
    message(FATAL_ERROR "CPU architecture not supported.")
endif ()





#choose mkl or openblas, if mkl present then use it, and it is depend on system too.










