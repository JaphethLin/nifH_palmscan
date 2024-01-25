#pragma once

// 70 < 120 < 130 < 140 > 150 > 160 > 210

const uint MIN_SEG = 70;	// Motifs �Ŀ����Сֵ	// Query �����ܳ���
const uint MAX_SEG = 210;	// Motifs �Ŀ�����ֵ

const uint MIN_SEG1 = 120;	//  70 < Motifs span < 120 ��
const uint MAX_SEG1 = 160;	//  160 < Motifs span < 210 ��
const double SEG_LENGTH_PENALTY1 = 10;	// �� 10 ��

const uint MIN_SEG2 = 130;	//  120 < Motifs span < 130 ��
const uint MAX_SEG2 = 150;	//  150 < Motifs span < 160 ��
const double SEG_LENGTH_PENALTY2 = 5;	// �� 5 ��
							//	130 < Motifs span < 150 ���� good-segment-length.
/* Lim deleted
const uint MIN_AB1 = 65;
const uint MAX_AB1 = 95;
const double AB_PENALTY1 = 5;

const uint MIN_BC1 = 8;
const uint MAX_BC1 = 45;
const double BC_PENALTY1 = 5;
*/

// 66-81

const uint MIN_AB1 = 66;	//  66 < dist(A->B) < 81 ���۷�
const uint MAX_AB1 = 81;	
const double AB_PENALTY1 = 5;	// ���򣬿� 5 ��

// 14-28

const uint MIN_BC1 = 14;		//  8 < dist(B->C) < 45 ���۷�
const uint MAX_BC1 = 28;
const double BC_PENALTY1 = 5;	// ���򣬿� 5 ��

const uint SHORT_QUERY = 100;	// Query length < 100 ���ǹ���

const double PERMUTED_PENALTY = 10;	// ������

const double myHICONF = 99;	// �����÷���ֵ
const double MIN_HICONF = 50;	// ��������űȶԵĵ÷���ֵ
const double MIN_LOCONF = 30;	// ��� Motif �÷� < MIN_PSSM_SCORE , ���ܵ÷� > MIN_HICONF , ����ֵΪ MIN_LOCONF
const double MIN_POSSIBLE = 5;	// ����

const double MIN_PSSM_SCORE = 3;	// ��� motif �÷� , ���ڸ�ֵ����ͷ�
const double MIN_PSSM_PENALTY = 5;	// ��� motif �÷ֵ��� MIN_PSSM_SCORE , �ܷ� - MIN_PSSM_PENALTY

const double MIN_C_SCORE = 2;		// motif C ����͵÷� , ������
const double MIN_PSSM_SCORE_PERMUTED = 4;	// ������

const double REWARD_DDGGDD = 10;	// ������
const double REWARD_DDGSDD = 10;	// ������
const double REWARD_DNMSDD = 10;	// ������
const double REWARD_DNxGDN = 8;		// ������
const double PENALTY_NOT_GDD_SDD_GDN = 10;	// ������

const uint MAX_X = 10;
