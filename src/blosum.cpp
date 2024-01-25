#include "myutils.h"

// P(xy|x), so probs sum to 1 one each row.
double g_BLOSUM62_JP_Rows[20][20] = 
	{
//              A       C       D       E       F       G       H       I       K       L       M       N       P       Q       R       S       T       V       W       Y
/* A */	{  0.2901, 0.0216, 0.0297, 0.0405, 0.0216, 0.0783, 0.0148, 0.0432, 0.0445, 0.0594, 0.0175, 0.0256, 0.0297, 0.0256, 0.0310, 0.0850, 0.0499, 0.0688, 0.0054, 0.0175, },  // A sum = 1.00
/* C */	{  0.0650, 0.4837, 0.0163, 0.0163, 0.0203, 0.0325, 0.0081, 0.0447, 0.0203, 0.0650, 0.0163, 0.0163, 0.0163, 0.0122, 0.0163, 0.0407, 0.0366, 0.0569, 0.0041, 0.0122, },  // C sum = 1.00
/* D */	{  0.0410, 0.0075, 0.3974, 0.0914, 0.0149, 0.0466, 0.0187, 0.0224, 0.0448, 0.0280, 0.0093, 0.0690, 0.0224, 0.0299, 0.0299, 0.0522, 0.0354, 0.0243, 0.0037, 0.0112, },  // D sum = 1.00
/* E */	{  0.0552, 0.0074, 0.0902, 0.2965, 0.0166, 0.0350, 0.0258, 0.0221, 0.0755, 0.0368, 0.0129, 0.0405, 0.0258, 0.0645, 0.0497, 0.0552, 0.0368, 0.0313, 0.0055, 0.0166, },  // E sum = 1.00
/* F */	{  0.0338, 0.0106, 0.0169, 0.0190, 0.3869, 0.0254, 0.0169, 0.0634, 0.0190, 0.1142, 0.0254, 0.0169, 0.0106, 0.0106, 0.0190, 0.0254, 0.0254, 0.0550, 0.0169, 0.0888, },  // F sum = 1.00
/* G */	{  0.0783, 0.0108, 0.0337, 0.0256, 0.0162, 0.5101, 0.0135, 0.0189, 0.0337, 0.0283, 0.0094, 0.0391, 0.0189, 0.0189, 0.0229, 0.0513, 0.0297, 0.0243, 0.0054, 0.0108, },  // G sum = 1.00
/* H */	{  0.0420, 0.0076, 0.0382, 0.0534, 0.0305, 0.0382, 0.3550, 0.0229, 0.0458, 0.0382, 0.0153, 0.0534, 0.0191, 0.0382, 0.0458, 0.0420, 0.0267, 0.0229, 0.0076, 0.0573, },  // H sum = 1.00
/* I */	{  0.0471, 0.0162, 0.0177, 0.0177, 0.0442, 0.0206, 0.0088, 0.2710, 0.0236, 0.1679, 0.0368, 0.0147, 0.0147, 0.0133, 0.0177, 0.0250, 0.0398, 0.1767, 0.0059, 0.0206, },  // I sum = 1.00
/* K */	{  0.0570, 0.0086, 0.0415, 0.0708, 0.0155, 0.0432, 0.0207, 0.0276, 0.2781, 0.0432, 0.0155, 0.0415, 0.0276, 0.0535, 0.1071, 0.0535, 0.0397, 0.0328, 0.0052, 0.0173, },  // K sum = 1.00
/* L */	{  0.0445, 0.0162, 0.0152, 0.0202, 0.0547, 0.0213, 0.0101, 0.1154, 0.0253, 0.3755, 0.0496, 0.0142, 0.0142, 0.0162, 0.0243, 0.0243, 0.0334, 0.0962, 0.0071, 0.0223, },  // L sum = 1.00
/* M */	{  0.0522, 0.0161, 0.0201, 0.0281, 0.0482, 0.0281, 0.0161, 0.1004, 0.0361, 0.1968, 0.1606, 0.0201, 0.0161, 0.0281, 0.0321, 0.0361, 0.0402, 0.0924, 0.0080, 0.0241, },  // M sum = 1.00
/* N */	{  0.0427, 0.0090, 0.0831, 0.0494, 0.0180, 0.0652, 0.0315, 0.0225, 0.0539, 0.0315, 0.0112, 0.3169, 0.0202, 0.0337, 0.0449, 0.0697, 0.0494, 0.0270, 0.0045, 0.0157, },  // N sum = 1.00
/* P */	{  0.0568, 0.0103, 0.0310, 0.0362, 0.0129, 0.0362, 0.0129, 0.0258, 0.0413, 0.0362, 0.0103, 0.0233, 0.4935, 0.0207, 0.0258, 0.0439, 0.0362, 0.0310, 0.0026, 0.0129, },  // P sum = 1.00
/* Q */	{  0.0559, 0.0088, 0.0471, 0.1029, 0.0147, 0.0412, 0.0294, 0.0265, 0.0912, 0.0471, 0.0206, 0.0441, 0.0235, 0.2147, 0.0735, 0.0559, 0.0412, 0.0353, 0.0059, 0.0206, },  // Q sum = 1.00
/* R */	{  0.0446, 0.0078, 0.0310, 0.0523, 0.0174, 0.0329, 0.0233, 0.0233, 0.1202, 0.0465, 0.0155, 0.0388, 0.0194, 0.0484, 0.3450, 0.0446, 0.0349, 0.0310, 0.0058, 0.0174, },  // R sum = 1.00
/* S */	{  0.1099, 0.0175, 0.0489, 0.0524, 0.0209, 0.0663, 0.0192, 0.0297, 0.0541, 0.0419, 0.0157, 0.0541, 0.0297, 0.0332, 0.0401, 0.2199, 0.0820, 0.0419, 0.0052, 0.0175, },  // S sum = 1.00
/* T */	{  0.0730, 0.0178, 0.0375, 0.0394, 0.0237, 0.0434, 0.0138, 0.0533, 0.0454, 0.0651, 0.0197, 0.0434, 0.0276, 0.0276, 0.0355, 0.0927, 0.2465, 0.0710, 0.0059, 0.0178, },  // T sum = 1.00
/* V */	{  0.0700, 0.0192, 0.0178, 0.0233, 0.0357, 0.0247, 0.0082, 0.1646, 0.0261, 0.1303, 0.0316, 0.0165, 0.0165, 0.0165, 0.0219, 0.0329, 0.0494, 0.2689, 0.0055, 0.0206, },  // V sum = 1.00
/* W */	{  0.0303, 0.0076, 0.0152, 0.0227, 0.0606, 0.0303, 0.0152, 0.0303, 0.0227, 0.0530, 0.0152, 0.0152, 0.0076, 0.0152, 0.0227, 0.0227, 0.0227, 0.0303, 0.4924, 0.0682, },  // W sum = 1.00
/* Y */	{  0.0405, 0.0093, 0.0187, 0.0280, 0.1308, 0.0249, 0.0467, 0.0436, 0.0312, 0.0685, 0.0187, 0.0218, 0.0156, 0.0218, 0.0280, 0.0312, 0.0280, 0.0467, 0.0280, 0.3178, },  // Y sum = 1.00
	};

double g_AAFreqs[20] =
        {
        0.07410, // A
        0.02460, // C
        0.05360, // D
        0.05430, // E
        0.04730, // F
        0.07410, // G
        0.02620, // H
        0.06790, // I
        0.05790, // K
        0.09880, // L
        0.02490, // M
        0.04450, // N
        0.03870, // P
        0.03400, // Q
        0.05160, // R
        0.05730, // S
        0.05070, // T
        0.07290, // V
        0.01320, // W
        0.03210, // Y
        }; // Sum = 0.99870


// python2 $py/fasta2aafreqs.py gitb_shortlabels.fa

double g_AAFreqs_RdRp[20] =
        {
        0.06412, // A
        0.02028, // C
        0.05610, // D
        0.05840, // E
        0.04536, // F
        0.05778, // G
        0.02448, // H
        0.05934, // I
        0.06129, // K
        0.09363, // L
        0.02618, // M
        0.04622, // N
        0.04563, // P
        0.03399, // Q
        0.05254, // R
        0.07285, // S
        0.06005, // T
        0.06894, // V
        0.01471, // W
        0.03810, // Y
        }; // sum = 1.00000
