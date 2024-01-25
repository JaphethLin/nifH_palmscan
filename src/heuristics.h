#pragma once

// 70 < 120 < 130 < 140 > 150 > 160 > 210

const uint MIN_SEG = 70;	// Motifs 的跨度最小值	// Query 序列总长度
const uint MAX_SEG = 210;	// Motifs 的跨度最大值

const uint MIN_SEG1 = 120;	//  70 < Motifs span < 120 或
const uint MAX_SEG1 = 160;	//  160 < Motifs span < 210 则
const double SEG_LENGTH_PENALTY1 = 10;	// 扣 10 分

const uint MIN_SEG2 = 130;	//  120 < Motifs span < 130 或
const uint MAX_SEG2 = 150;	//  150 < Motifs span < 160 则
const double SEG_LENGTH_PENALTY2 = 5;	// 扣 5 分
							//	130 < Motifs span < 150 则标记 good-segment-length.
/* Lim deleted
const uint MIN_AB1 = 65;
const uint MAX_AB1 = 95;
const double AB_PENALTY1 = 5;

const uint MIN_BC1 = 8;
const uint MAX_BC1 = 45;
const double BC_PENALTY1 = 5;
*/

// 66-81

const uint MIN_AB1 = 66;	//  66 < dist(A->B) < 81 不扣分
const uint MAX_AB1 = 81;	
const double AB_PENALTY1 = 5;	// 否则，扣 5 分

// 14-28

const uint MIN_BC1 = 14;		//  8 < dist(B->C) < 45 不扣分
const uint MAX_BC1 = 28;
const double BC_PENALTY1 = 5;	// 否则，扣 5 分

const uint SHORT_QUERY = 100;	// Query length < 100 则标记过短

const double PERMUTED_PENALTY = 10;	// 已弃用

const double myHICONF = 99;	// 鉴定得分阈值
const double MIN_HICONF = 50;	// 不错的置信比对的得分阈值
const double MIN_LOCONF = 30;	// 最低 Motif 得分 < MIN_PSSM_SCORE , 但总得分 > MIN_HICONF , 均赋值为 MIN_LOCONF
const double MIN_POSSIBLE = 5;	// 弃用

const double MIN_PSSM_SCORE = 3;	// 最低 motif 得分 , 低于该值给予惩罚
const double MIN_PSSM_PENALTY = 5;	// 最低 motif 得分低于 MIN_PSSM_SCORE , 总分 - MIN_PSSM_PENALTY

const double MIN_C_SCORE = 2;		// motif C 得最低得分 , 已弃用
const double MIN_PSSM_SCORE_PERMUTED = 4;	// 已弃用

const double REWARD_DDGGDD = 10;	// 已弃用
const double REWARD_DDGSDD = 10;	// 已弃用
const double REWARD_DNMSDD = 10;	// 已弃用
const double REWARD_DNxGDN = 8;		// 已弃用
const double PENALTY_NOT_GDD_SDD_GDN = 10;	// 已弃用

const uint MAX_X = 10;
