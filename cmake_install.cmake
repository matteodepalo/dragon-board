# Install script for directory: /Users/matteodepalo/Documents/Università/dragon-board

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    FILE(INSTALL DESTINATION "/Users/matteodepalo/Documents/Università/dragon-board/lib" TYPE SHARED_LIBRARY FILES "/Users/matteodepalo/Documents/Università/dragon-board/lib/Debug/libdragon-board.so")
    IF(EXISTS "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
        -id "libdragon-board.so"
        "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDIF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  IF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    FILE(INSTALL DESTINATION "/Users/matteodepalo/Documents/Università/dragon-board/lib" TYPE SHARED_LIBRARY FILES "/Users/matteodepalo/Documents/Università/dragon-board/lib/Release/libdragon-board.so")
    IF(EXISTS "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
        -id "libdragon-board.so"
        "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDIF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  IF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    FILE(INSTALL DESTINATION "/Users/matteodepalo/Documents/Università/dragon-board/lib" TYPE SHARED_LIBRARY FILES "/Users/matteodepalo/Documents/Università/dragon-board/lib/MinSizeRel/libdragon-board.so")
    IF(EXISTS "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
        -id "libdragon-board.so"
        "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDIF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  IF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    FILE(INSTALL DESTINATION "/Users/matteodepalo/Documents/Università/dragon-board/lib" TYPE SHARED_LIBRARY FILES "/Users/matteodepalo/Documents/Università/dragon-board/lib/RelWithDebInfo/libdragon-board.so")
    IF(EXISTS "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
        -id "libdragon-board.so"
        "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/Users/matteodepalo/Documents/Università/dragon-board/lib/libdragon-board.so")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDIF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/Users/matteodepalo/Documents/Università/dragon-board/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/Users/matteodepalo/Documents/Università/dragon-board/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
