#ifndef __ROT_BWT_H
#define __ROT_BWT_H

/*
 * =====================================================================================
 *
 *       Filename:  rbwt.h
 *
 *    Description:  rotation bwt
 *
 *        Version:  1.0
 *        Created:  05/19/2020 05:14:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#define __set_bwt(bwt, l, c) ((bwt)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_bwt(bwt, l) ((bwt)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __set_bwt2(bwt, l, c) ((bwt)[(l)>>1] |= (c)<<((~(l)&1)<<2))
#define __get_bwt2(bwt, l) ((bwt)[(l)>>1] >> ((~(l)&1)<<2)&15)

#define __get_col(seq, i) (((seq)>>((i)*4))&3)
#define  MIN_BWT_SIZE 16
#define  NO_BWT_SUM_SIZE 64


#endif
