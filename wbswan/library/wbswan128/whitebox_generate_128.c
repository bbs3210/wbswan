#include <stdio.h>
#include <math/affine_transform.h>
#include <math/matrix_utils.h>
#include <wbswan/wbswan.h>
#include <time.h>

 #define DELTA 0xb9f4567924888695
//#define DELTA 0x00

static unsigned char S[16] = {0x01, 0x02, 0x0C, 0x05, 0x07, 0x08, 0x0A, 0x0F, 0x04, 0x0D, 0x0B, 0x0E, 0x09, 0x06, 0x00, 0x03};

void InvRotateKeyByte(uint8_t *key, uint16_t keylength)
{
    uint8_t i;
    uint8_t temp[7];
    uint8_t N = keylength / 8 - 1;
    for (i = 0; i < 7; i++)
    {
        temp[6 - i] = key[N - i];
    }

    //Right rotate every byte of the key;
    for (i = N; i >= 7; i--)
    {
        key[i] = key[i - 7];
    }

    //Cyclic the first byte of the key to the MSB;
    for (i = 0; i < 7; i++)
    {
        key[i] = temp[i];
    }
}

void RotateKeyByte(uint8_t *key, uint16_t keylength)
{
    uint8_t i;
    uint8_t temp[7];
    uint8_t N = keylength / 8 - 1;
    for (i = 0; i < 7; i++)
    {
        temp[i] = key[i];
    }

    //Right rotate every byte of the key;
    for (i = 0; i <= N - 7; i++)
    {
        key[i] = key[i + 7];
    }

    //Cyclic the first byte of the key to the MSB;
    for (i = 0; i < 7; i++)
    {
        key[(N + 1) - 7 + i] = temp[i];
    }
}

void AddRoundConstant(uint16_t *subkey, uint64_t sum)
{

    uint64_t *a = (uint64_t *)subkey;
    uint16_t *b = (uint16_t *)a;
    (*a) = (*a) + sum;

    subkey[0] = b[0];
    subkey[1] = b[1];
    subkey[2] = b[2];
    subkey[3] = b[3];
}
int Key_Schedule(unsigned char *Seedkey, int KeyLen, unsigned char Direction, unsigned char *Subkey, int rounds)
{
    int i;
    uint16_t *key;
    uint16_t subkey[4];
    uint64_t round_constant = 0;
    uint16_t(*ekey)[4] = (uint16_t(*)[4])malloc(sizeof(uint16_t) * (rounds)*4);
    key = (uint16_t *)malloc(sizeof(uint16_t) * (KeyLen / 8));
    memcpy(key, Seedkey, (KeyLen / 8) * sizeof(uint8_t));
    MatGf2 mat = make_transposition_64(64);
    MatGf2 mat_back = make_transposition_back_64(64);
    for (i = 0; i < rounds; i++)
    {
        RotateKeyByte(key, KeyLen);
        subkey[0] = key[0];
        subkey[1] = key[1];
        subkey[2] = key[2];
        subkey[3] = key[3];
        round_constant = round_constant + DELTA;
        AddRoundConstant(subkey, round_constant);

        *((uint64_t *)subkey) = ApplyMatToU64(mat, *((uint64_t *)subkey));

        ekey[i][0] = subkey[0];
        ekey[i][1] = subkey[1];
        ekey[i][2] = subkey[2];
        ekey[i][3] = subkey[3];

        *((uint64_t *)subkey) = ApplyMatToU64(mat_back, *((uint64_t *)subkey));

        key[0] = subkey[0];
        key[1] = subkey[1];
        key[2] = subkey[2];
        key[3] = subkey[3];
    }
    memcpy(Subkey, ekey, (rounds)*4 * sizeof(uint16_t));
    free(key);
    free(ekey);
    MatGf2Free(mat);
    MatGf2Free(mat_back);
    mat = NULL;
    mat_back = NULL;
    return 0;
}

int _swan_whitebox_helper_init(swan_whitebox_helper *swh, uint8_t *key, int weak_or_strong)
{

    int block_size = swh->block_size;
    int semi_block = block_size / 2;
    swh->piece_count = semi_block / SWAN_PIECE_BIT;
    swh->key = key;
    swh->weak_or_strong = weak_or_strong;

    return 0;
}

int swan_whitebox_128_strong_helper_init(const uint8_t *key, swan_whitebox_helper *swh, int enc)
// designed for swan_64_128
{

    // uint8_t rk[MAX_RK_SIZE];
    int ret;
    swh->cfg = CFG;
    swh->encrypt = enc;
    swh->rounds = swan_cfg_rounds[swh->cfg];
    swh->block_size = swan_cfg_blocksizes[swh->cfg];
    return _swan_whitebox_helper_init(swh, key, 1);
}

int swan_whitebox_128_weak_helper_init(const uint8_t *key, swan_whitebox_helper *swh, int enc)
// designed for swan_128_128
{

    // uint8_t rk[MAX_RK_SIZE];
    int ret;
    swh->cfg = CFG;
    swh->encrypt = enc;
    swh->rounds = swan_cfg_rounds[swh->cfg];
    swh->block_size = swan_cfg_blocksizes[swh->cfg];
    return _swan_whitebox_helper_init(swh, key, 1);
}

int _swan_whitebox_content_init(swan_whitebox_helper *swh, swan_whitebox_content *swc)
{
    // TODO:
    swc->cfg = swh->cfg;
    swc->rounds = swh->rounds;
    swc->block_size = swh->block_size;
    swc->piece_count = swh->piece_count;

    swc->lut = (swan_wb_semi(*)[4][2][256])malloc(sizeof(swan_wb_semi) * swc->rounds * 256 * 4 * 2);
    swc->P = (CombinedAffine *)malloc((swh->rounds + 2) * sizeof(CombinedAffine));
    swc->B = (CombinedAffine *)malloc(swh->rounds * sizeof(CombinedAffine));
    swc->C = (CombinedAffine *)malloc(swh->rounds * sizeof(CombinedAffine));
    swc->weak_or_strong = swh->weak_or_strong;
    int i;
    CombinedAffine *B_ptr = swc->B;
    CombinedAffine *C_ptr = swc->C;
    CombinedAffine *P_ptr = swc->P;

    if (swc->weak_or_strong)
    {
        swc->PQ = (CombinedAffine *)malloc(swh->rounds * sizeof(CombinedAffine));
        CombinedAffine *PQ_ptr = swc->PQ;
        for (i = 0; i < (swh->rounds); i++)
        {
            combined_affine_init(PQ_ptr++, SWAN_PIECE_BIT, swh->piece_count);
        }
    }
    else
    {
        swc->PQ = NULL;
    }

    for (i = 0; i < (swh->rounds); i++)
    {

        combined_affine_init(B_ptr++, SWAN_PIECE_BIT / 2, swh->piece_count * 2);
        combined_affine_init(C_ptr++, SWAN_PIECE_BIT, swh->piece_count);
    }
    for (i = 0; i < (swh->rounds + 2); i++)
    {
        // combined_affine_init(P_ptr++, SWAN_PIECE_BIT / 2, swh->piece_count * 2);
        combined_affine_init(P_ptr++, SWAN_PIECE_BIT, swh->piece_count);
    }

    return 0;
}

int _swan_whitebox_content_assemble(swan_whitebox_helper *swh, swan_whitebox_content *swc)
{
    int piece_count = swh->piece_count;
    int j, k, r;
    uint16_t i;
    int piece_range = 1 << SWAN_PIECE_BIT;

    /*input encoding and output encoding*/
    CombinedAffine *P_ptr = swc->P;
    CombinedAffine *end_Ptr = swc->P + swh->rounds;
    // int p = 0;
    // for (k = 0; k < 16; k++)
    // {
    //     if (k == 8)
    //     {
    //         P_ptr++;
    //         end_Ptr++;
    //     }
    //     for (p = 0; p < 256; p++)
    //     {

    //         swc->SE[k][p] = ApplyAffineToU8(P_ptr->sub_affine[k % 8], p);
    //         swc->EE[k][p] = ApplyAffineToU8(end_Ptr->sub_affine_inv[k % 8], p);
    //     }
    // }

    /* theta */
    MatGf2 rotate = make_transposition_64(piece_count * SWAN_PIECE_BIT);
    // MatGf2 rotate = make_right_rotate_shift_64(piece_count * SWAN_PIECE_BIT, ROL_A[swh->cfg], ROL_B[swh->cfg], ROL_C[swh->cfg]);
    // MatGf2Mul(special, rotate, &rotate);

    CombinedAffine *B_ptr = swc->B;
    CombinedAffine *C_ptr = swc->C;
    MatGf2 temp;

    CombinedAffine *PQ_ptr = swc->PQ;

    for (i = 0; i < swh->rounds; i++)
    {
        //B * rotate * P'* X + B * rotate * p + b

        temp = GenMatGf2Mul(B_ptr->combined_affine->linear_map, rotate);
        B_ptr->combined_affine->linear_map = GenMatGf2Mul(temp, P_ptr->combined_affine_inv->linear_map);
        B_ptr->combined_affine->vector_translation = GenMatGf2Add(GenMatGf2Mul(temp, P_ptr->combined_affine_inv->vector_translation), B_ptr->combined_affine->vector_translation);
        P_ptr++;
        B_ptr++;
    }

    B_ptr = swc->B;
    P_ptr = swc->P + 1;
    MatGf2 rotate_back = make_transposition_back_64(piece_count * SWAN_PIECE_BIT);
    rotate = make_right_rotate_shift_64(piece_count * SWAN_PIECE_BIT, ROL_A[swh->cfg], ROL_B[swh->cfg], ROL_C[swh->cfg]);
    swan_wb_unit key_schedule[64][4];
    Key_Schedule(swh->key, swan_cfg_keysizes[swh->cfg], 1, key_schedule, swh->rounds);

    srand(time(NULL));
    /* w */
    PQ_ptr = swc->PQ;
    if (swh->encrypt == 1)
    {
        for (r = 0; r < swh->rounds; r++)
        {

            uint16_t randomNum[4][4][2] = {0};

            int j = 0;
            int t = 0;
            // random
            int num = 0;

            int y;
            for (j = 0; j < 4; j++)
            {
                for (t = 0; t < 4; t++)
                {
                    for (y = 0; y < 2; y++)
                    {
                        randomNum[j][t][y] = (uint16_t)(rand() % 50);
                    }
                }
            }

            for (k = 0; k < piece_count; k++)
            {
                uint16_t switchxor = 0x0000;
                for (t = 0; t < 4; t++)
                {

                    if (t != k)
                    {

                        switchxor = switchxor ^ randomNum[t][t][0] ^ randomNum[t][t][1];
                    }
                }

                P_ptr->sub_affine[k].vector_translation = GenMatAddU16(P_ptr->sub_affine[k].vector_translation, ApplyMatToU16(P_ptr->sub_affine[k].linear_map, switchxor));

                for (i = 0; i < 256; i++)
                {

                    int n;
                    for (n = 0; n < 2; n++)
                    {
                        uint8_t t8;
                        if (n == 0)
                        {
                            t8 = ApplyAffineToU8((B_ptr)->sub_affine_inv[2 * k], i) ^ ((key_schedule[r][k]) && 0x00ff);
                        }
                        else
                        {
                            t8 = ApplyAffineToU8((B_ptr)->sub_affine_inv[2 * k + 1], i) ^ ((key_schedule[r][k] >> 8) && 0x00ff);
                        }
                        uint16_t yc[4] = {0};
                        t8 = (S[(t8 >> 4) & 0x0f]) << 4 | (S[t8 & 0x0f]);
                        int t = 0;

                        if (n == 1)
                        {
                            uint16_t at;
                            at = t8;
                            yc[k] = (at << 8) & 0xff00;
                        }
                        else
                        {
                            yc[k] = ((uint16_t)t8) & 0x00ff;
                        }

                        *((uint64_t *)yc) = ApplyMatToU64(rotate_back, *((uint64_t *)yc));
                        *((uint64_t *)yc) = ApplyMatToU64(rotate, *((uint64_t *)yc));

                        if (k != 3)
                        {
                            for (j = 0; j < 4; j++)
                            {
                                yc[j] = yc[j] ^ randomNum[k][j][n];
                            }
                        }
                        else
                        {

                            yc[0] = yc[0] ^ randomNum[1][0][n] ^ randomNum[2][0][n];
                            yc[1] = yc[1] ^ randomNum[0][1][n] ^ randomNum[2][1][n];
                            yc[2] = yc[2] ^ randomNum[0][2][n] ^ randomNum[1][2][n];
                            yc[3] = yc[3] ^ randomNum[0][3][n] ^ randomNum[1][3][n] ^ randomNum[2][3][n] ^ randomNum[3][3][n];
                        }

                        swc->lut[r][k][n][i] = *((uint64_t *)yc);
                    }
                }
            }
            B_ptr++;
            P_ptr++;
        }
    }
    else
    {
        for (r = swh->rounds - 1; r >= 0; r--)
        {

            uint16_t randomNum[4][4][2];

            int j = 0;
            int t = 0;
            // random
            int num = 0;
            int y;
            for (j = 0; j < 4; j++)
            {
                for (t = 0; t < 4; t++)
                {
                    for (y = 0; y < 2; y++)
                    {
                        randomNum[j][t][y] = 0;
                    }
                }
            }

            for (k = 0; k < piece_count; k++)
            {
                uint16_t switchxor = 0;
                for (t = 0; t < 4; t++)
                {

                    if (t != k)
                    {

                        switchxor = switchxor ^ randomNum[t][t][0] ^ randomNum[t][t][1];
                    }
                }

                P_ptr->sub_affine[k].vector_translation = GenMatAddU16(P_ptr->sub_affine[k].vector_translation, ApplyMatToU16(P_ptr->sub_affine[k].linear_map, switchxor));

                for (i = 0; i < 256; i++)
                {
                    int n;

                    for (n = 0; n < 2; n++)
                    {
                        uint8_t t8;
                        if (n == 0)
                        {
                            t8 = ApplyAffineToU8((B_ptr)->sub_affine_inv[2 * k], i) ^ ((key_schedule[r][k]) && 0x00ff);
                        }
                        else
                        {
                            t8 = ApplyAffineToU8((B_ptr)->sub_affine_inv[2 * k + 1], i) ^ ((key_schedule[r][k] >> 8) && 0x00ff);
                        }

                        uint16_t yc[4] = {0};
                        t8 = (S[(t8 >> 4) & 0x0f]) << 4 | (S[t8 & 0x0f]);
                        int t = 0;
                        if (n == 1)
                        {
                            uint16_t at;
                            at = t8;
                            yc[k] = (at << 8) & 0xff00;
                        }
                        else
                        {
                            yc[k] = ((uint16_t)t8) & 0x00ff;
                        }

                        *((uint64_t *)yc) = ApplyMatToU64(rotate_back, *((uint64_t *)yc));
                        *((uint64_t *)yc) = ApplyMatToU64(rotate, *((uint64_t *)yc));

                        if (k != 3)
                        {
                            for (j = 0; j < 4; j++)
                            {
                                yc[j] = yc[j] ^ randomNum[k][j][n];
                            }
                        }
                        else
                        {
                            yc[0] = yc[0] ^ randomNum[1][0][n] ^ randomNum[2][0][n];
                            yc[1] = yc[1] ^ randomNum[0][1][n] ^ randomNum[2][1][n];
                            yc[2] = yc[2] ^ randomNum[0][2][n] ^ randomNum[1][2][n];
                            yc[3] = yc[3] ^ randomNum[0][3][n] ^ randomNum[1][3][n] ^ randomNum[2][3][n] ^ randomNum[3][3][n];
                        }

                        // rotate = make_right_rotate_shift(piece_count * SWAN_PIECE_BIT, ROL_A[swh->cfg], ROL_B[swh->cfg], ROL_C[swh->cfg]);

                        swc->lut[r][k][n][i] = *((uint64_t *)yc);
                    }
                }
            }
            B_ptr++;
            P_ptr++;
        }
    }

    // MatGf2Free(special);
    // special = NULL;
    MatGf2Free(rotate);
    rotate = NULL;
    return 0;
}

int swan_whitebox_128_init(const uint8_t *key, int enc, swan_whitebox_content *swc)
{
    swan_whitebox_helper *swh = (swan_whitebox_helper *)malloc(sizeof(swan_whitebox_helper));
//    swan_whitebox_128_weak_helper_init(key, swh, enc);
    swan_whitebox_128_strong_helper_init(key, swh, enc);
    swan_whitebox_128_content_init(swh, swc);
    swan_whitebox_128_content_assemble(swh, swc);

    return 0;
}
int swan_whitebox_128_content_init(swan_whitebox_helper *swh, swan_whitebox_content *swc)
{
    return _swan_whitebox_content_init(swh, swc);
}
int swan_whitebox_128_content_assemble(swan_whitebox_helper *swh, swan_whitebox_content *swc)
{
    return _swan_whitebox_content_assemble(swh, swc);
}
