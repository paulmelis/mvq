CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT(MVQ)

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/..")

ADD_EXECUTABLE(t_quat "t_quat.cpp")