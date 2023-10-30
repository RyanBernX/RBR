/*
 * ==========================================================================
 *
 *       Filename:  api_macro.h
 *
 *    Description:  header defines API macros
 *
 *        Version:  1.0
 *        Created:  08/31/2021 05:15:50 PM
 *       Revision:  none
 *       Compiler:  gcc/g++
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ==========================================================================
 */

#ifndef RBR_API_MACRO_H
#define RBR_API_MACRO_H

// 导出共享库符号的相关宏
#if defined(__WIN32__) || defined(__CYGWIN__)
    #define RBR_HELPER_DLL_IMPORT __declspec(dllimport)
    #define RBR_HELPER_DLL_EXPORT __declspec(dllexport)
    #define RBR_HELPER_DLL_LOCAL
#elif (defined(__GNUC__) && __GNUC__ >= 4) || defined(__clang__)
    #define RBR_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define RBR_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define RBR_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
    #define RBR_HELPER_DLL_IMPORT
    #define RBR_HELPER_DLL_EXPORT
    #define RBR_HELPER_DLL_LOCAL
#endif

#ifdef arrabit_EXPORTS
    #define RBR_API RBR_HELPER_DLL_EXPORT
#else
    #define RBR_API RBR_HELPER_DLL_IMPORT
#endif
#define RBR_LOCAL RBR_HELPER_DLL_LOCAL

#endif

