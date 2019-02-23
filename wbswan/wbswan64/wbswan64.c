#include <stdio.h>
#include <math/affine_transform.h>
#include <math/matrix_utils.h>
#include <wbswan/wbswan.h>
void dump(const uint8_t *li, int len)
{
    int line_ctrl = 16;
    int i;
    for (i = 0; i < len; i++)
    {
        printf("%02X", (*li++));
        if ((i + 1) % line_ctrl == 0)
        {
            printf("\n");
        }
        else
        {
            printf(" ");
        }
    }
}

static MatGf2 make_special_rotate_back(int dim)
{
    uint8_t rot[8] = {24, 28, 16, 20, 8, 12, 0, 4};
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    int j;
    int row = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 8; j++)
        {
            MatGf2Set(ind, row++, rot[j] + i, 1);
        }
    }
    return ind;
}

static MatGf2 make_special_rotate(int dim)
{
    uint8_t rot[8] = {0, 8, 16, 24, 1, 9, 17, 25};
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    int j;
    int row;
    for (i = 1; i <= 4; i++)
    {
        row = 32 - 8 * i;

        for (j = 0; j < 8; j++)
        {
            MatGf2Set(ind, row++, rot[j] + 2 * (i - 1), 1);
        }
    }
    return ind;
}
static MatGf2 make_swithlane_mat(int dim)
{
    MatGf2 ind = GenMatGf2(dim, dim);
    int i;
    for (i = 0; i < 8; i++)
    {
        MatGf2Set(ind, i, i + 8, 1);
        MatGf2Set(ind, i, i + 16, 1);
        MatGf2Set(ind, i, i + 24, 1);
    }
    for (i = 8; i < 16; i++)
    {
        MatGf2Set(ind, i, i - 8, 1);
        MatGf2Set(ind, i, i + 16, 1);
        MatGf2Set(ind, i, i + 8, 1);
    }
    for (i = 16; i < 24; i++)
    {
        MatGf2Set(ind, i, i + 8, 1);
        MatGf2Set(ind, i, i - 16, 1);
        MatGf2Set(ind, i, i - 8, 1);
    }
    for (i = 24; i < 32; i++)
    {
        MatGf2Set(ind, i, i - 8, 1);
        MatGf2Set(ind, i, i - 16, 1);
        MatGf2Set(ind, i, i - 24, 1);
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

int swan_whitebox_encrypt(const uint8_t *in, uint8_t *out, swan_whitebox_content *swc)
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
    swan_wb_semi(*lut_ptr)[4][256] = swc->lut;
    MatGf2 switchmat = make_swithlane_mat(32);
    LT[0] = in[0];
    LT[1] = in[1];
    LT[2] = in[2];
    LT[3] = in[3];

    RT[0] = in[4];
    RT[1] = in[5];
    RT[2] = in[6];
    RT[3] = in[7];

    uint8_t tt[4];
    uint8_t dd[4];
    tempL = (swan_wb_semi *)LT;
    tempR = (swan_wb_semi *)RT;

    *tempL = ApplyAffineToU32(*(P_ptr->combined_affine), *tempL);
    *tempR = ApplyAffineToU32(*((P_ptr + 1)->combined_affine), *tempR);

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
        *tempL = ApplyAffineToU32(*(B_ptr->combined_affine), *tempL);

        //start beta
        //back to 8 in 8 out
        //******************************************************************


        //******************************************************************

        //******************************************************************
        uint32_t LC[4];
    

        LC[0] = lut_ptr[i][0][LT[0]];
        LC[1] = lut_ptr[i][1][LT[1]];
        LC[2] = lut_ptr[i][2][LT[2]];
        LC[3] = lut_ptr[i][3][LT[3]];
    
        *tempL = LC[0] ^ LC[1] ^ LC[2] ^ LC[3];
   
        //SWITCHLANE
        *tempL = MATtoU32(switchmat, *tempL);

        LT[0] = ApplyAffineToU8((C_ptr)->sub_affine[0], LT[0]);
        LT[1] = ApplyAffineToU8((C_ptr)->sub_affine[1], LT[1]);
        LT[2] = ApplyAffineToU8((C_ptr)->sub_affine[2], LT[2]);
        LT[3] = ApplyAffineToU8((C_ptr)->sub_affine[3], LT[3]);

        *tempL = ApplyAffineToU32(*((C_ptr)->combined_affine_inv), *tempL);

        //******************************************************************

        // P for next

        *tempL = ApplyAffineToU32(*((P_ptr + 1)->combined_affine), *tempL);

        swan_wb_unit tempd[4];

        *((swan_wb_semi *)tempd) = ApplyAffineToU32(*((P_ptr + 2)->combined_affine), ApplyAffineToU32(*((P_ptr)->combined_affine_inv), *((uint32_t *)L)));

        L[0] = R[0] ^ LT[0];
        L[1] = R[1] ^ LT[1];
        L[2] = R[2] ^ LT[2];
        L[3] = R[3] ^ LT[3];

        R[0] = tempd[0];
        R[1] = tempd[1];
        R[2] = tempd[2];
        R[3] = tempd[3];

        MatGf2 mat_x = NULL;
        ReAllocatedMatGf2(32, 1, &mat_x);
        InitVecFromBit(*((uint32_t *)L), mat_x);
        MatGf2Add(GenMatGf2Mul((P_ptr + 1)->combined_affine->linear_map, (P_ptr + 1)->combined_affine_inv->vector_translation), mat_x, &mat_x);
        *((swan_wb_semi *)L) = get32FromVec(mat_x);

        B_ptr++;
        P_ptr++;
        C_ptr++;
    }

    *((swan_wb_semi *)L) = ApplyAffineToU32(*((P_ptr + 0)->combined_affine_inv), *((swan_wb_semi *)L));
    *((swan_wb_semi *)R) = ApplyAffineToU32(*((P_ptr + 1)->combined_affine_inv), *((swan_wb_semi *)R));

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

int swan_whitebox_decrypt(const uint8_t *in, uint8_t *out, swan_whitebox_content *swc)
{
    int i;
    swan_wb_unit L[4];
    swan_wb_unit R[4];
    swan_wb_unit LT[4];
    swan_wb_unit RT[4];
    swan_wb_semi *tempL;
    swan_wb_semi *tempR;
 
    uint8_t dd[4];
    uint8_t tt[4];
    CombinedAffine *B_ptr = swc->B;
    CombinedAffine *C_ptr = swc->C;
    CombinedAffine *P_ptr = swc->P;
    MatGf2 switchmat = make_swithlane_mat(32);
    swan_wb_semi(*lut_ptr)[4][256] = swc->lut;

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

    *tempR = ApplyAffineToU32(*(P_ptr->combined_affine), *tempR);
    *tempL = ApplyAffineToU32(*((P_ptr + 1)->combined_affine), *tempL);

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

        // start theta
        *tempR = ApplyAffineToU32(*(B_ptr->combined_affine), *tempR);

        //***************************************************************************

        //**********************************************************************************

        //***********************************************************

        //start beta
        uint32_t RC[4];
        RC[0] = lut_ptr[swc->rounds - 1 - (i)][0][RT[0]];
        RC[1] = lut_ptr[swc->rounds - 1 - (i)][1][RT[1]];
        RC[2] = lut_ptr[swc->rounds - 1 - (i)][2][RT[2]];
        RC[3] = lut_ptr[swc->rounds - 1 - (i)][3][RT[3]];

        *tempR = RC[0] ^ RC[1] ^ RC[2] ^ RC[3];

        //SWITCHLANE

        *tempR = MATtoU32(switchmat, *tempR);

        RT[0] = ApplyAffineToU8((C_ptr)->sub_affine[0], RT[0]);
        RT[1] = ApplyAffineToU8((C_ptr)->sub_affine[1], RT[1]);
        RT[2] = ApplyAffineToU8((C_ptr)->sub_affine[2], RT[2]);
        RT[3] = ApplyAffineToU8((C_ptr)->sub_affine[3], RT[3]);

        *tempR = ApplyAffineToU32(*((C_ptr)->combined_affine_inv), *tempR);

        // P for next
        *tempR = ApplyAffineToU32(*((P_ptr + 1)->combined_affine), *tempR);
        swan_wb_unit tempd[4];

        *((swan_wb_semi *)tempd) = ApplyAffineToU32(*((P_ptr + 2)->combined_affine), ApplyAffineToU32(*((P_ptr)->combined_affine_inv), *((uint32_t *)R)));

        R[0] = L[0] ^ RT[0];
        R[1] = L[1] ^ RT[1];
        R[2] = L[2] ^ RT[2];
        R[3] = L[3] ^ RT[3];

        L[0] = tempd[0];
        L[1] = tempd[1];
        L[2] = tempd[2];
        L[3] = tempd[3];

        MatGf2 mat_x = NULL;
        ReAllocatedMatGf2(32, 1, &mat_x);
        InitVecFromBit(*((swan_wb_semi *)R), mat_x);
        MatGf2Add(GenMatGf2Mul((P_ptr + 1)->combined_affine->linear_map, (P_ptr + 1)->combined_affine_inv->vector_translation), mat_x, &mat_x);
        *((swan_wb_semi *)R) = get32FromVec(mat_x);

        P_ptr++;
        B_ptr++;
        C_ptr++;
    }

    *((swan_wb_semi *)R) = ApplyAffineToU32(*((P_ptr + 0)->combined_affine_inv), *((swan_wb_semi *)R));
    *((swan_wb_semi *)L) = ApplyAffineToU32(*((P_ptr + 1)->combined_affine_inv), *((swan_wb_semi *)L));

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

int main()
{
    uint8_t in[8] = {0x01, 0x02, 0x03, 0x01, 0x05, 0x06, 0x07, 0x08};
    uint8_t key[16] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
    uint8_t out[8] = {0xC8, 0xEB, 0x79, 0x08, 0x85, 0xD3, 0x14, 0xFB};
    swan_whitebox_content swc_enc;
    swan_whitebox_content swc_dec;
    swan_whitebox_64_init(key, 1, &swc_enc);
    swan_whitebox_encrypt(in, out, &swc_enc);
    swan_whitebox_release(&swc_enc);
    dump(out, sizeof(out));
    printf("\n");
    swan_whitebox_64_init(key, 0, &swc_dec);
    swan_whitebox_decrypt(out, in, &swc_dec);
    swan_whitebox_release(&swc_dec);
    dump(in, sizeof(in));
    printf("\n");

}