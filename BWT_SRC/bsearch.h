/*
 * =====================================================================================
 *
 *       Filename:  bsearch.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017年10月19日 11时35分19秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef __BSEARCH_H
#define __BSEARCH_H

#include <stdint.h>

#define bs_iter_t int64_t
#define bs_val_t uint32_t

bs_iter_t lower_bound(bs_iter_t first, bs_iter_t last, bs_val_t *a, bs_val_t val);
bs_iter_t uper_bound(bs_iter_t first, bs_iter_t last, bs_val_t *a, bs_val_t val);

#endif
