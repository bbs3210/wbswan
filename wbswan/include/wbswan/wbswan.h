#ifndef _WBSWAN_H_
#define _WBSWAN_H_
#include <math/affine_transform.h>
#include <math/matrix_utils.h>
#define SWAN_PIECE_BIT 8
#define CFG swan_cfg_B128_K128
enum swan_cipher_config_t
{
    swan_cfg_B64_K128,
    swan_cfg_B64_K256,
    swan_cfg_B128_K128,
    swan_cfg_B128_K256,
    swan_cfg_B256_K256,
};

#define MAX_RK_SIZE 64

#if SWAN_PIECE_BIT == 8
typedef uint32_t swan_wb_semi;
typedef uint8_t swan_wb_unit;


#elif SWAN_PIECE_BIT == 16
typedef uint64_t swan_wb_semi;
typedef uint16_t swan_wb_unit;

#elif SWAN_PIECE_BIT == 32
typedef uint128_t swan_wb_semi;
typedef uint32_t swan_wb_unit;

#endif

static int swan_cfg_rounds[] = {64, 1, 1, 64, 64};

static int swan_cfg_keysizes[] = {128, 256, 128, 256, 256};

static int swan_cfg_blocksizes[] = {64, 64, 128, 128, 256};

static MatGf2 make_special_rotate(int dim);



static int ROL_A[] = {1, 1, 1};
static int ROL_B[] = {3, 3, 3};
static int ROL_C[] = {5, 5, 5};

typedef struct swan_whitebox_helper
{
    enum swan_cipher_config_t cfg;
    int aff_in_round;
    int block_size;
    int rounds;
    int piece_count;
    int encrypt;
    CombinedAffine *P;
    CombinedAffine *B;
    CombinedAffine *C;
    CombinedAffine *D;
    uint8_t *key;
} swan_whitebox_helper;

typedef struct swan_wb_t
{
    enum swan_cipher_config_t cfg;
    uint32_t rounds;
    uint32_t block_size;
    uint32_t piece_count; // piece_count = block_size / 8, every 8 bit combined as a piece
    uint64_t (*lut)[4][2][256];
    CombinedAffine *P;
    CombinedAffine *B;
    CombinedAffine *C;



} swan_whitebox_content;

#define ROL16(x, n) ((x >> n) | (x << (16 - n)))

int swan_whitebox_64_init(const uint8_t *key, int enc, swan_whitebox_content *swc);



#endif
