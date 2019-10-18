/*
 * =====================================================================================
 *
 *       Filename:  debug.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年04月02日 16时35分25秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#define fprintf __fprintf_no
static inline int __fprintf_no(FILE *stream, const char *format, ...)
{
    return 0;
}
