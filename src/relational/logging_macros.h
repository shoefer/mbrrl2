#ifndef LOGGING_MACROS_H
#define LOGGING_MACROS_H

#ifdef ROS_LOGGING
#include <ros/ros.h>
#include <ros/console.h>

#define RLOG_DEBUG         ROS_DEBUG
#define RLOG_DEBUG_STR     ROS_DEBUG_STREAM
#define RLOG_INFO          ROS_INFO
#define RLOG_INFO_STR      ROS_INFO_STREAM
#define RLOG_WARN       ROS_WARN
#define RLOG_WARN_STR   ROS_WARN_STREAM
#define RLOG_ERROR         ROS_ERROR
#define RLOG_ERROR_STR     ROS_ERROR_STREAM

#define RLOG_DEBUG_NAMED   ROS_DEBUG_NAMED
#define RLOG_DEBUG_STR_NAMED   ROS_DEBUG_STREAM_NAMED
#define RLOG_INFO_NAMED   ROS_INFO_NAMED
#define RLOG_INFO_STR_NAMED   ROS_INFO_STREAM_NAMED
#define RLOG_WARN_NAMED   ROS_WARN_NAMED
#define RLOG_WARN_STR_NAMED   ROS_WARN_STREAM_NAMED
#define RLOG_ERROR_NAMED   ROS_ERROR_NAMED
#define RLOG_ERROR_STR_NAMED   ROS_ERROR_STREAM_NAMED

#else
#include <stdio.h>
#include <iostream>

// no debug output if not in debug mode:
#ifdef NDEBUG
#define RLOG_DEBUG(...)       (void)0
#define RLOG_DEBUG_STR(...)   (void)0
#define RLOG_DEBUG_NAMED(...)       (void)0
#define RLOG_DEBUG_STR_NAMED(...)   (void)0

#else
//#define RLOG_DEBUG(...)        fprintf(stderr, __VA_ARGS__), fflush(stderr)
//#define RLOG_DEBUG_STR(args)   std::cerr << args << std::endl
#define RLOG_DEBUG(...)      fprintf(stderr, "[DEBUG] "), fprintf(stdout, __VA_ARGS__), fflush(stdout)
#define RLOG_DEBUG_STR(args) std::cout << "[DEBUG] " << args << std::endl
#define RLOG_DEBUG_NAMED(...)   fprintf(stderr, "[DEBUG] "), fprintf(stdout, __VA_ARGS__), fflush(stderr)
#define RLOG_DEBUG_STR_NAMED(name, args)   std::cout << "[DEBUG] " << "[" << name << "] " << args << std::endl
#endif


#ifdef NINFO
#define RLOG_INFO(...)       (void)0
#define RLOG_INFO_STR(...)   (void)0
#define RLOG_INFO_NAMED(...)   (void)0
#define RLOG_INFO_STR_NAMED(name, args)   (void)0
#else
#define RLOG_INFO(...)         fprintf(stdout, __VA_ARGS__), fflush(stderr)
#define RLOG_INFO_STR(args)    std::cout << args << std::endl
#define RLOG_INFO_NAMED(...)   fprintf(stdout, __VA_ARGS__), fflush(stderr)
#define RLOG_INFO_STR_NAMED(name, args)   std::cout << "[" << name << "] " << args << std::endl
#endif

#ifdef NWARN
#define RLOG_WARN(...)       (void)0
#define RLOG_WARN_STR(...)   (void)0
#define RLOG_WARN_NAMED(...)   (void)0
#define RLOG_WARN_STR_NAMED(name, args)   (void)0
#else
#define RLOG_WARN(...)      fprintf(stderr, "WARNING: "), fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define RLOG_WARN_STR(args) std::cerr << "WARNING: " << args << std::endl
#define RLOG_WARN_NAMED(...)   fprintf(stderr, "WARNING: "), fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define RLOG_WARN_STR_NAMED(name, args)   std::cerr << "[" << name << "] " << "WARNING: " << args << std::endl
#endif

#define RLOG_ERROR(...)        fprintf(stderr, "ERROR: "), fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define RLOG_ERROR_STR(args)   std::cerr << "ERROR: " << args << std::endl
#define RLOG_ERROR_NAMED(...)   fprintf(stderr, "ERROR: "), fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define RLOG_ERROR_STR_NAMED(name, args)   std::cerr << "[" << name << "] " << "ERROR: " << args << std::endl
#endif

#endif // LOGGING_MACROS_H
