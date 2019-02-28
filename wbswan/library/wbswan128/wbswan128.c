#include <stdio.h>
#include <math/affine_transform.h>
#include <math/matrix_utils.h>
#include <wbswan/wbswan.h>
MatGf2 make_swithlane_64(int dim)
{
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    for (i = 0; i < 16; i++)
    {
        MatGf2Set(ind, i, i + 16, 1);
        MatGf2Set(ind, i, i + 32, 1);
        MatGf2Set(ind, i, i + 48, 1);
    }
    for (i = 16; i < 32; i++)
    {
        MatGf2Set(ind, i, i - 16, 1);
        MatGf2Set(ind, i, i + 32, 1);
        MatGf2Set(ind, i, i + 16, 1);
    }
    for (i = 32; i < 48; i++)
    {
        MatGf2Set(ind, i, i + 16, 1);
        MatGf2Set(ind, i, i - 32, 1);
        MatGf2Set(ind, i, i - 16, 1);
    }
    for (i = 48; i < 64; i++)
    {
        MatGf2Set(ind, i, i - 16, 1);
        MatGf2Set(ind, i, i - 32, 1);
        MatGf2Set(ind, i, i - 48, 1);
    }

    return ind;
}

MatGf2 make_right_rotate_shift_64(int dim, int r1, int r2, int r3)
{
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    for (i = 0; i < 16; i++)
    {
        MatGf2Set(ind, i, ((i + 16 - r3) % 16 + 0), 1);
    }
    for (i = 16; i < 32; i++)
    {
        MatGf2Set(ind, i, (i + 16 - r2) % 16 + 16, 1);
    }
    for (i = 32; i < 48; i++)
    {
        MatGf2Set(ind, i, (i + 16 - r1) % 16 + 32, 1);
    }
    for (i = 48; i < 64; i++)
    {
        MatGf2Set(ind, i, (i + 16 - 0) % 16 + 48, 1);
    }

    return ind;
}
MatGf2 make_transposition_64(int dim)
{
    uint8_t rot[16] = {0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51};
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    int j;
    int row;
    for (i = 1; i <= 4; i++)
    {
        row = 64 - 16 * i;

        for (j = 0; j < 16; j++)
        {
            MatGf2Set(ind, row++, rot[j] + 4 * (i - 1), 1);
        }
    }
    return ind;
}

MatGf2 make_transposition_back_64(int dim)
{
    uint8_t rot[16] = {48, 52, 56, 60, 32, 36, 40, 44, 16, 20, 24, 28, 0, 4, 8, 12};
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    int j;
    int row = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 16; j++)
        {
            MatGf2Set(ind, row++, rot[j] + i, 1);
        }
    }
    return ind;
}

void tobin(uint32_t x, int SIZE) //要和函数声明一致，所以后面也要填int x,int a[]
{
    int a[32] = {0};
    int n = 0, t, k;
    do
    {
        a[n] = x % 2;
        x = (unsigned int)x >> 1; //要考虑到参数int x有可能为负数所以填x=x/2是不行的。
        //如果x>=0的话填x=x/2没有问题，实际上我估计这里出题者的本意希望填x/2，但是如果当x为负数的时候
        //会出错的，要么填 x=(unsigned int)x/2也是可以的，不过 x=(unsigned int)x/2的运行效率远远低于x=(unsigned int)x>>1。这里牵涉的东西比较多，三言两语说不清
        //如果想了解原因，建议回去看谭浩强的强制类型转换、正数和负数的2进制表示、移位3个知识点
        n++;
    } while (x != 0);
    //do...while()这个功能就是把这个数的二进制的位存入这个数组中
    for (k = 0; k < n / 2; k++)
    {
        t = a[k];
        a[k] = a[n - k - 1]; //实现数组中2个数交换
        a[n - k - 1] = t;
        //for循环是为了交换顺序，比如x=11是的二进制码是1011这4个码一次存在a[3] a[2] a[1] a[0]中，而输出的时候是按a[0] a[1] a[2] a[3]这样输出的如果没有这个交换屏幕上看到的会是1101
    }
    for (k = SIZE; k > n; k--)
    {
        printf("%d-", a[k]);
    }
    for (k = 0; k < n; k++)
        printf("%d-", a[k]);
    int num = 0;
    for (int i = 0; i < 32; i++)
    {
        num += a[i];
    }
    // printf("%d\n", num);
    printf("\n");
}

int swan_whitebox_encrypt(const uint16_t *in, uint16_t *out, swan_whitebox_content *swc)
{
    int i;
    swan_wb_unit L[4];
    swan_wb_unit R[4];
    swan_wb_unit LT[4];
    swan_wb_unit RT[4];
    swan_wb_semi *tempL;
    swan_wb_semi *tempR;

    CombinedAffine *B_ptr = swc->B;
    CombinedAffine *C_ptr = swc->C;
    CombinedAffine *P_ptr = swc->P;
    CombinedAffine *PQ_ptr = swc->PQ;
    swan_wb_semi(*lut_ptr)[4][2][256] = swc->lut;
    MatGf2 switchmat = make_swithlane_64(64);
    LT[0] = in[0];
    LT[1] = in[1];
    LT[2] = in[2];
    LT[3] = in[3];

    RT[0] = in[4];
    RT[1] = in[5];
    RT[2] = in[6];
    RT[3] = in[7];

    uint16_t tt[4];
    uint16_t dd[4];
    tempL = (swan_wb_semi *)LT;
    tempR = (swan_wb_semi *)RT;

    *tempL = ApplyAffineToU64(*(P_ptr->combined_affine), *tempL);
    *tempR = ApplyAffineToU64(*((P_ptr + 1)->combined_affine), *tempR);

    L[0] = LT[0];
    L[1] = LT[1];
    L[2] = LT[2];
    L[3] = LT[3];

    R[0] = RT[0];
    R[1] = RT[1];
    R[2] = RT[2];
    R[3] = RT[3];
    for (i = 0; i < swc->rounds; i++)
    {

        LT[0] = L[0];
        LT[1] = L[1];
        LT[2] = L[2];
        LT[3] = L[3];

        // start theta
        *tempL = ApplyAffineToU64(*(B_ptr->combined_affine), *tempL);

        // start beta

        uint64_t LC[4][2];
        uint64_t temp = 0x0000000000000000;

        int j;
        int x;
        for (j = 0; j < 4; j++)
        {
            for (x = 0; x < 2; x++)
            {
                LC[j][x] = lut_ptr[i][j][x][x ? ((LT[j] >> 8) & 0x00ff) : LT[j] & 0x00ff];
                temp ^= LC[j][x];
            }
        }

        *tempL = temp;

        //************************************

        //SWITCHLANE
        *tempL = ApplyMatToU64(switchmat, *tempL);


        //******************************************************************

        // P for next

        CombinedAffine *tca = NULL;
        if (i == 0)
        {
            tca = P_ptr + 1;
        }
        else
        {
            tca = P_ptr + 1;
        }
        LT[0] = ApplyAffineToU16((tca)->sub_affine[0], LT[0]);
        LT[1] = ApplyAffineToU16((tca)->sub_affine[1], LT[1]);
        LT[2] = ApplyAffineToU16((tca)->sub_affine[2], LT[2]);
        LT[3] = ApplyAffineToU16((tca)->sub_affine[3], LT[3]);
        #if (DEBUG)
            printf("\t%08X\n", *tempL);
        #endif

        swan_wb_unit tempd[4];

        *((swan_wb_semi *)tempd) = ApplyAffineToU64(*((P_ptr + 2)->combined_affine), ApplyAffineToU64(*((P_ptr)->combined_affine_inv), *((swan_wb_semi *)L)));

        L[0] = R[0] ^ LT[0];
        L[1] = R[1] ^ LT[1];
        L[2] = R[2] ^ LT[2];
        L[3] = R[3] ^ LT[3];

        R[0] = tempd[0];
        R[1] = tempd[1];
        R[2] = tempd[2];
        R[3] = tempd[3];

        swan_wb_semi LP;
        if (swc->weak_or_strong)
        {
            LT[0] = L[0];
            LT[1] = L[1];
            LT[2] = L[2];
            LT[3] = L[3];
            LP = *tempL;
            swan_wb_unit ULP[4];
            swan_wb_semi *ulp_ptr = (swan_wb_semi *)ULP;
            *ulp_ptr = ApplyAffineToU64(*(PQ_ptr->combined_affine), ApplyMatToU64(((P_ptr + 1)->combined_affine_inv->linear_map), *tempL));
            *((swan_wb_semi *)L) = *ulp_ptr;

            MatGf2 mat_x = NULL;
            ReAllocatedMatGf2(64, 1, &mat_x);
            InitVecFromBit(LP, mat_x);
            MatGf2Add(GenMatGf2Mul((tca)->combined_affine->linear_map, (tca)->combined_affine_inv->vector_translation), mat_x, &mat_x);
            LP = get64FromVec(mat_x);

        } else {

            MatGf2 mat_x = NULL;
            ReAllocatedMatGf2(64, 1, &mat_x);
            InitVecFromBit(*((swan_wb_semi *)L), mat_x);
            MatGf2Add(GenMatGf2Mul((tca)->combined_affine->linear_map, (tca)->combined_affine_inv->vector_translation), mat_x, &mat_x);
            *((swan_wb_semi *)L) = get64FromVec(mat_x);
        }



        B_ptr++;
        P_ptr++;
        i++;

        // next round

        if (swc->weak_or_strong)
        {
            *tempL = LP;
        }
        else
        {
            LT[0] = L[0];
            LT[1] = L[1];
            LT[2] = L[2];
            LT[3] = L[3];
        }

        // start theta
        *tempL = ApplyAffineToU64(*(B_ptr->combined_affine), *tempL);

        //start beta
        temp = 0x0000000000000000;
        for (j = 0; j < 4; j++)
        {
            for (x = 0; x < 2; x++)
            {
                LC[j][x] = lut_ptr[i][j][x][x ? ((LT[j] >> 8) & 0x00ff) : LT[j] & 0x00ff];
                temp ^= LC[j][x];
            }
        }

        *tempL = temp;
        // start switch
        *tempL = ApplyMatToU64(switchmat, *tempL);

        LT[0] = ApplyAffineToU16((P_ptr + 1)->sub_affine[0], LT[0]);
        LT[1] = ApplyAffineToU16((P_ptr + 1)->sub_affine[1], LT[1]);
        LT[2] = ApplyAffineToU16((P_ptr + 1)->sub_affine[2], LT[2]);
        LT[3] = ApplyAffineToU16((P_ptr + 1)->sub_affine[3], LT[3]);

        if (swc->weak_or_strong)
        {
            *((swan_wb_semi *)tempd) = ApplyAffineToU64(*((P_ptr + 2)->combined_affine), ApplyAffineToU64(*((PQ_ptr)->combined_affine_inv), *((swan_wb_semi *)L)));
        }
        else
        {
            *((swan_wb_semi *)tempd) = ApplyAffineToU64(*((P_ptr + 2)->combined_affine), ApplyAffineToU64(*((P_ptr)->combined_affine_inv), *((swan_wb_semi *)L)));
        }

        L[0] = R[0] ^ LT[0];
        L[1] = R[1] ^ LT[1];
        L[2] = R[2] ^ LT[2];
        L[3] = R[3] ^ LT[3];

        R[0] = tempd[0];
        R[1] = tempd[1];
        R[2] = tempd[2];
        R[3] = tempd[3];

        MatGf2 mat_x = NULL;
        ReAllocatedMatGf2(64, 1, &mat_x);
        InitVecFromBit(*((uint64_t *)L), mat_x);
        MatGf2Add(GenMatGf2Mul((P_ptr + 1)->combined_affine->linear_map, (P_ptr + 1)->combined_affine_inv->vector_translation), mat_x, &mat_x);
        *((swan_wb_semi *)L) = get64FromVec(mat_x);

        B_ptr++;
        P_ptr++;
        PQ_ptr++;
    }

    *((swan_wb_semi *)L) = ApplyAffineToU64(*((P_ptr + 0)->combined_affine_inv), *((swan_wb_semi *)L));
    *((swan_wb_semi *)R) = ApplyAffineToU64(*((P_ptr + 1)->combined_affine_inv), *((swan_wb_semi *)R));

    out[0] = L[0];
    out[1] = L[1];
    out[2] = L[2];
    out[3] = L[3];

    out[4] = R[0];
    out[5] = R[1];
    out[6] = R[2];
    out[7] = R[3];

    return 0;
}

int swan_whitebox_decrypt(const uint16_t *in, uint16_t *out, swan_whitebox_content *swc)
{
    int i;
    swan_wb_unit L[4];
    swan_wb_unit R[4];
    swan_wb_unit LT[4];
    swan_wb_unit RT[4];
    swan_wb_semi *tempL;
    swan_wb_semi *tempR;

    CombinedAffine *B_ptr = swc->B;
    CombinedAffine *C_ptr = swc->C;
    CombinedAffine *P_ptr = swc->P;
    CombinedAffine *PQ_ptr = swc->PQ;
    MatGf2 switchmat = make_swithlane_64(64);
    swan_wb_semi(*lut_ptr)[4][2][256] = swc->lut;

    LT[0] = in[0];
    LT[1] = in[1];
    LT[2] = in[2];
    LT[3] = in[3];

    RT[0] = in[4];
    RT[1] = in[5];
    RT[2] = in[6];
    RT[3] = in[7];

    tempL = (swan_wb_semi *)LT;
    tempR = (swan_wb_semi *)RT;

    *tempR = ApplyAffineToU64(*(P_ptr->combined_affine), *tempR);
    *tempL = ApplyAffineToU64(*((P_ptr + 1)->combined_affine), *tempL);

    L[0] = LT[0];
    L[1] = LT[1];
    L[2] = LT[2];
    L[3] = LT[3];

    R[0] = RT[0];
    R[1] = RT[1];
    R[2] = RT[2];
    R[3] = RT[3];

    for (i = 0; i < swc->rounds; i++)
    {

        RT[0] = R[0];
        RT[1] = R[1];
        RT[2] = R[2];
        RT[3] = R[3];

        #if (DEBUG)
            printf("from:\t%016X\n", *tempL);
        #endif

        // start theta
        *tempR = ApplyAffineToU64(*(B_ptr->combined_affine), *tempR);

        #if (DEBUG)
            printf("theta:\t%016X\n", *tempL);
        #endif
        //start beta
        uint64_t RC[4][2];

        uint64_t temp = 0x0000000000000000;
        int j;
        int x;
        for (j = 0; j < 4; j++)
        {
            for (x = 0; x < 2; x++)
            {
                RC[j][x] = lut_ptr[swc->rounds - 1 - (i)][j][x][x ? ((RT[j] >> 8) & 0x00ff) : RT[j] & 0x00ff];
                temp ^= RC[j][x];
            }
        }

        *tempR = temp;

        #if (DEBUG)
            printf("beta:\t%016X\n", *tempL);
        #endif

        //SWITCHLANE

        *tempR = ApplyMatToU64(switchmat, *tempR);

        // P for next

        #if (DEBUG)
            printf("\t%08X\n", *tempR);
        #endif
        CombinedAffine * tca = NULL;
        if (i==0) {
            tca = P_ptr+1;
        } else {
            tca = P_ptr+1;
        }
        RT[0] = ApplyAffineToU16((tca)->sub_affine[0], RT[0]);
        RT[1] = ApplyAffineToU16((tca)->sub_affine[1], RT[1]);
        RT[2] = ApplyAffineToU16((tca)->sub_affine[2], RT[2]);
        RT[3] = ApplyAffineToU16((tca)->sub_affine[3], RT[3]);
        #if (DEBUG)
            printf("\t%016X\n", *tempL);
        #endif
        swan_wb_unit tempd[4];

        *((swan_wb_semi *)tempd) = ApplyAffineToU64(*((P_ptr + 2)->combined_affine), ApplyAffineToU64(*((P_ptr)->combined_affine_inv), *((uint64_t *)R)));

        R[0] = L[0] ^ RT[0];
        R[1] = L[1] ^ RT[1];
        R[2] = L[2] ^ RT[2];
        R[3] = L[3] ^ RT[3];

        L[0] = tempd[0];
        L[1] = tempd[1];
        L[2] = tempd[2];
        L[3] = tempd[3];

        swan_wb_semi LP;
        MatGf2 mat_lp = NULL;
        if (swc->weak_or_strong) {
            RT[0] = R[0];
            RT[1] = R[1];
            RT[2] = R[2];
            RT[3] = R[3];
            LP = *tempR;
            swan_wb_unit ULP[4];
            swan_wb_semi *ulp_ptr = (swan_wb_semi *)ULP;
            *ulp_ptr = ApplyAffineToU64(*(PQ_ptr->combined_affine), ApplyMatToU64(((P_ptr + 1)->combined_affine_inv->linear_map), *tempR));
            *((swan_wb_semi *)R) = *ulp_ptr;

            MatGf2 mat_x = NULL;
            ReAllocatedMatGf2(64, 1, &mat_x);
            InitVecFromBit(LP, mat_x);
            MatGf2Add(GenMatGf2Mul((tca)->combined_affine->linear_map, (tca)->combined_affine_inv->vector_translation), mat_x, &mat_x);
            LP = get64FromVec(mat_x);

        } else {

            MatGf2 mat_x = NULL;
            ReAllocatedMatGf2(64, 1, &mat_x);
            InitVecFromBit(*((swan_wb_semi *)R), mat_x);
            MatGf2Add(GenMatGf2Mul((tca)->combined_affine->linear_map, (tca)->combined_affine_inv->vector_translation), mat_x, &mat_x);
            *((swan_wb_semi *)R) = get64FromVec(mat_x);

        }
        #if (DEBUG)
            printf("Round %d:\n", i);
            dump(L, 8);
            printf("\n");
            dump(R, 8);
            printf("\n");
            printf("%016X\n", LP);
        #endif

        P_ptr++;
        B_ptr++;
        i++;

        // next round

        if (swc->weak_or_strong)
        {
            *tempR = LP;
        }
        else
        {
            RT[0] = R[0];
            RT[1] = R[1];
            RT[2] = R[2];
            RT[3] = R[3];
        }

        #if (DEBUG)
            printf("before theta: \n\t%016X\n", *tempL);
        #endif

        *tempR = ApplyAffineToU64(*(B_ptr->combined_affine), *tempR);

        #if (DEBUG)
            printf("after theta: \n\t%016X\n", *tempL);
        #endif
        //start beta
        temp = 0x0000000000000000;
        for (j = 0; j < 4; j++)
        {
            for (x = 0; x < 2; x++)
            {
                RC[j][x] = lut_ptr[swc->rounds - 1 - (i)][j][x][x ? ((RT[j] >> 8) & 0x00ff) : RT[j] & 0x00ff];
                temp ^= RC[j][x];
            }
        }

        *tempR = temp;
        // start switch
        *tempR = ApplyMatToU64(switchmat, *tempR);

        RT[0] = ApplyAffineToU16((P_ptr + 1)->sub_affine[0], RT[0]);
        RT[1] = ApplyAffineToU16((P_ptr + 1)->sub_affine[1], RT[1]);
        RT[2] = ApplyAffineToU16((P_ptr + 1)->sub_affine[2], RT[2]);
        RT[3] = ApplyAffineToU16((P_ptr + 1)->sub_affine[3], RT[3]);

        if (swc->weak_or_strong)
        {
            *((swan_wb_semi *)tempd) = ApplyAffineToU64(*((P_ptr + 2)->combined_affine), ApplyAffineToU64(*((PQ_ptr)->combined_affine_inv), *((swan_wb_semi *)R)));
        }
        else
        {
            *((swan_wb_semi *)tempd) = ApplyAffineToU64(*((P_ptr + 2)->combined_affine), ApplyAffineToU64(*((P_ptr)->combined_affine_inv), *((swan_wb_semi *)R)));
        }


        R[0] = L[0] ^ RT[0];
        R[1] = L[1] ^ RT[1];
        R[2] = L[2] ^ RT[2];
        R[3] = L[3] ^ RT[3];

        L[0] = tempd[0];
        L[1] = tempd[1];
        L[2] = tempd[2];
        L[3] = tempd[3];

        MatGf2 mat_x = NULL;
        ReAllocatedMatGf2(64, 1, &mat_x);
        InitVecFromBit(*((swan_wb_semi *)R), mat_x);
        MatGf2Add(GenMatGf2Mul((P_ptr + 1)->combined_affine->linear_map, (P_ptr + 1)->combined_affine_inv->vector_translation), mat_x, &mat_x);
        *((swan_wb_semi *)R) = get64FromVec(mat_x);

        B_ptr++;
        P_ptr++;
        PQ_ptr++;

        #if (DEBUG)
            printf("Round %d:\n", i);
            dump(L, 8);
            printf("\n");
            dump(R, 8);
            printf("\n");
        #endif
    }

    *((swan_wb_semi *)R) = ApplyAffineToU64(*((P_ptr + 0)->combined_affine_inv), *((swan_wb_semi *)R));
    *((swan_wb_semi *)L) = ApplyAffineToU64(*((P_ptr + 1)->combined_affine_inv), *((swan_wb_semi *)L));

    out[0] = L[0];
    out[1] = L[1];
    out[2] = L[2];
    out[3] = L[3];

    out[4] = R[0];
    out[5] = R[1];
    out[6] = R[2];
    out[7] = R[3];

    return 0;
}

int swan_whitebox_release(swan_whitebox_content *swc)
{
    // TODO:
    // AffineTransformFree(swc->)

    int i, j;
    CombinedAffine *B_ptr = swc->B;
    for (i = 0; i < swc->rounds; i++)
    {
        combined_affine_free(B_ptr++);
    }
    free(swc->B);
    swc->B = NULL;
    CombinedAffine *C_ptr = swc->C;
    for (i = 0; i < swc->rounds; i++)
    {
        combined_affine_free(C_ptr++);
    }
    free(swc->C);
    swc->C = NULL;
    CombinedAffine *P_ptr = swc->P;
    for (i = 0; i < swc->rounds + 2; i++)
    {
        combined_affine_free(P_ptr++);
    }
    free(swc->P);
    swc->P = NULL;

    //free memory
    free(swc->lut);
    swc->lut = NULL;

    return 0;
}
