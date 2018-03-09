/* copied from Jack Mott                                             *
 * https://gist.github.com/jackmott/b7268510227303bb28756596a09cd1e2 */

// A header file to get you set going with Intel SIMD instrinsic programming. 
// All necessary header files are inlucded for SSE2, SSE41, and AVX2
// Macros make the intrinsics easier to read and  generic so you can compile to 
// SSE2 or AVX2 with the flip of a #define

 

//#define SSE2  //indicates we want SSE2
//#define SSE41 //indicates we want SSE4.1 instructions (floor and blend is available)
//#define AVX2  //indicates we want AVX2 instructions (double speed!) 


#ifndef __AVX2__
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE 2
#endif

#ifdef __SSE4_1__
#include <smmintrin.h> // SSE4.1
#endif

#ifdef __AVX2__
#include <immintrin.h> // AVX2
#endif

// #include <zmmintrin.h> //avx512 the world is not yet ready...SOON


// create types we can use in either the 128 or 256 case
#ifndef __AVX2__
// m128 will be our base type
typedef __m128 SIMD;   //for floats
typedef __m128i SIMDi; //for integers

// we process 4 at a time
//#define VECTOR_SIZE 4
//#define MEMORY_ALIGNMENT 16

// intrinsic functions
#define Store(x,y) _mm_store_ps(x,y)
#define Load(x) _mm_load_ps(x)
#define Storei(x,y) _mm_store_si128(x,y)
#define Loadi(x) _mm_load_si128(x)
#define Seti64(x,y) __mm_set_epi64(x,y)
#define Seti32(x1,x2,x3,x4) __mm_set_epi32(x1,x2,x3,x4)
#define Seti16(x1,x2,x3,x4,x5,x6,x7,x8) __mm_set_epi16(x1,x2,x3,x4,x5,x6,x7,x8)
#define Seti8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16) __mm_set_epi8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)
#define SetOne(x) _mm_set1_ps(x)
#define SetZero() _mm_setzero_ps()
#define SetOnei64(x) _mm_set1_epi64(x)
#define SetOnei32(x) _mm_set1_epi32(x)
#define SetOnei16(x) _mm_set1_epi16(x)
#define SetOnei8(x) _mm_set1_epi8(x)
#define SetZeroi(x) _mm_setzero_epi32(x)
#define Add(x,y) _mm_add_ps(x,y)
#define Sub(x,y) _mm_sub_ps(x,y)
#define Addi64(x,y) _mm_add_epi64(x,y)
#define Subi64(x,y) _mm_sub_epi64(x,y)
#define Addi32(x,y) _mm_add_epi32(x,y)
#define Subi32(x,y) _mm_sub_epi32(x,y)
#define Addi16(x,y) _mm_add_epi16(x,y)
#define Subi16(x,y) _mm_sub_epi16(x,y)
#define Addi8(x,y) _mm_add_epi8(x,y)
#define Subi8(x,y) _mm_sub_epi8(x,y)
#define Mul(x,y) _mm_mul_ps(x,y)
#define Muli(x,y) _mm_mul_epi32(x,y)
#define And(x,y) _mm_and_ps(x,y)
#define Andi(x,y) _mm_and_si128(x,y)
#define AndNot(x,y) _mm_andnot_ps(x,y)
#define AndNoti(x,y) _mm_andnot_si128(x,y)
#define Or(x,y) _mm_or_ps(x,y)
#define Ori(x,y) _mm_or_si128(x,y)
#define ShiftLeft(x,y) _mm_slli_si128(x,y) // shift whole x by y bits
#define ShiftRight(x,y) _mm_srli_si128(x,y)
#define CastToFloat(x) _mm_castsi128_ps(x)
#define CastToInt(x) _mm_castps_si128(x)
#define ConvertToInt(x) _mm_cvtps_epi32(x)
#define ConvertToFloat(x) _mm_cvtepi32_ps(x)
#define ConvertFrom32To128 _mm_cvtsi32_si128(x)
#define ConvertFrom64To128 _mm_cvtsi64_si128(x)
#define ConvertFrom128To64 _mm_cvtsi128_si64(x)
#define ConvertFrom128To32 _mm_cvtsi128_si32(x)
#define Equal(x,y)  _mm_cmpeq_ps(x,y) 
#define Equali64(x,y) _mm_cmpeq_epi64(x,y)
#define Equali32(x,y) _mm_cmpeq_epi32(x,y)
#define Equali16(x,y) _mm_cmpeq_epi16(x,y)
#define Equali8(x,y) _mm_cmpeq_epi8(x,y)
#define GreaterThan(x,y) _mm_cmpgt_ps(x,y)
#define GreaterThani64(x,y) _mm_cmpgt_epi64(x,y)
#define GreaterThani32(x,y) _mm_cmpgt_epi32(x,y)
#define GreaterThani16(x,y) _mm_cmpgt_epi16(x,y)
#define GreaterThani8(x,y) _mm_cmpgt_epi8(x,y)
#define GreaterThanOrEq(x,y) _mm_cmpge_ps(x,y)
#define LessThan(x,y) _mm_cmplt_ps(x,y)
#define LessThani64(x,y) _mm_cmplt_epi64(y,x) 
#define LessThani32(x,y) _mm_cmplt_epi32(y,x) 
#define LessThani16(x,y) _mm_cmplt_epi16(y,x) 
#define LessThani8(x,y) _mm_cmplt_epi8(y,x) 
#define LessThanOrEq(x,y) _mm_cmple_ps(x,y)
#define LessThanOrEqi64(x,y) _mm_cmple_epi64(x,y)
#define LessThanOrEqi32(x,y) _mm_cmple_epi32(x,y)
#define LessThanOrEqi16(x,y) _mm_cmple_epi16(x,y)
#define LessThanOrEqi8(x,y) _mm_cmple_epi8(x,y)
#define NotEqual(x,y) _mm_cmpneq_ps(x,y)
#define NotEquali64(x,y) _mm_cmpneq_epi64(x,y)
#define NotEquali32(x,y) _mm_cmpneq_epi32(x,y)
#define NotEquali16(x,y) _mm_cmpneq_epi16(x,y)
#define NotEquali8(x,y) _mm_cmpneq_epi8(x,y)
#define Max(x,y) _mm_max_ps(x,y)
#define Maxi64(x,y) _mm_max_epi64(x,y)
#define Maxi32(x,y) _mm_max_epi32(x,y)
#define Maxi16(x,y) _mm_max_epi16(x,y)
#define Maxi8(x,y) _mm_max_epi8(x,y)
#define Min(x,y) _mm_min_ps(x,y)
#define Mini64(x,y) _mm_min_epi64(x,y)
#define Mini32(x,y) _mm_min_epi32(x,y)
#define Mini16(x,y) _mm_min_epi16(x,y)
#define Mini8(x,y) _mm_min_epi8(x,y)
#ifndef __SSE4_1__
#define BlendV(x,y,z) Or(AndNot(z,x), And(z,y))   //if we don't have sse4
#else  // __SSE4_1__
#define BlendV(x,y,z) _mm_blendv_ps(x,y,z)	
#define BlendVi32(x,y,z) _mm_blend_epi32(x,y,z)	
#define BlendVi16(x,y,z) _mm_blend_epi16(x,y,z)	
#define BlendVi8(x,y,z) _mm_blendv_epi8(x,y,z)	
#define Floor(x) _mm_floor_ps(x)
#endif // __SSE4_1__
	
#endif // no __AVX2__
#ifdef __AVX2__

// m256 will be our base type
typedef __m256 SIMD;  //for floats
typedef __m256i SIMDi; //for integers

//process 8 at t time
//#define VECTOR_SIZE 8
//#define MEMORY_ALIGNMENT 32

//intrinsic functions
//TODO storeu, loadu
#define Store(x,y) _mm256_store_ps(x,y)
#define Load(x) _mm256_load_ps(x)
#define Storei(x,y) _mm256_store_si256(x,y)
#define Loadi(x) _mm256_load_si256(x)
#define Set(x,y,z,w,a,b,c,d) _mm256_set_ps(x,y,z,w,a,b,c,d);
#define Seti128(x,y) __mm256_set_m128(x,y)
#define Seti64(x1,x2,x3,x4) __mm256_set_epi64x(x1,x2,x3,x4)
#define Seti32(x1,x2,x3,x4,x5,x6,x7,x8) __mm256_set_epi32(x1,x2,x3,x4,x5,x6,x7,x8)
#define Seti16(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16) __mm256_set_epi16(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)
#define Seti8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32) __mm256_set_epi8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32)
#define SetOne(x) _mm256_set1_ps(x)
#define SetZero() _mm256_setzero_ps()
#define SetOnei64(x) _mm256_set1_epi64x(x)
#define SetOnei32(x) _mm256_set1_epi32(x)
#define SetOnei16(x) _mm256_set1_epi16(x)
#define SetOnei8(x) _mm256_set1_epi8(x)
#define SetZeroi(x) _mm256_setzero_epi32(x)
#define Add(x,y) _mm256_add_ps(x,y)
#define Sub(x,y) _mm256_sub_ps(x,y)
#define Addi64(x,y) _mm256_add_epi64(x,y)
#define Subi64(x,y) _mm256_sub_epi64(x,y)
#define Addi32(x,y) _mm256_add_epi32(x,y)
#define Subi32(x,y) _mm256_sub_epi32(x,y)
#define Addi16(x,y) _mm256_add_epi16(x,y)
#define Subi16(x,y) _mm256_sub_epi16(x,y)
#define Addi8(x,y) _mm256_add_epi8(x,y)
#define Subi8(x,y) _mm256_sub_epi8(x,y)
#define Mul(x,y) _mm256_mul_ps(x,y)
#define Muli(x,y) _mm256_mul_epi32(x,y)
#define And(x,y) _mm256_and_ps(x,y)
#define Andi(x,y) _mm256_and_si256(x,y)
#define AndNot(x,y) _mm256_andnot_ps(x,y)
#define AndNoti(x,y) _mm256_andnot_si256(x,y)
#define Or(x,y) _mm256_or_ps(x,y)
#define Ori(x,y) _mm256_or_si256(x,y)
#define ShiftLeft(x,y) _mm256_slli_si256(x,y)
#define ShiftRight(x,y) _mm256_srli_si256(x,y)
#define CastToFloat(x) _mm256_castsi256_ps(x)
#define CastToInt(x) _mm256_castps_si256(x)
#define ConvertToInt(x) _mm256_cvtps_epi32(x)
#define ConvertToFloat(x) _mm256_cvtepi32_ps(x)
#define Equal(x,y)  _mm256_cmp_ps(x,y,_CMP_EQ_OQ) 
#define Equali64(x,y) _mm256_cmpeq_epi64(x,y)
#define Equali32(x,y) _mm256_cmpeq_epi32(x,y)
#define Equali16(x,y) _mm256_cmpeq_epi16(x,y)
#define Equali8(x,y) _mm256_cmpeq_epi8(x,y)
#define GreaterThan(x,y) _mm256_cmp_ps(x,y,_CMP_GT_OQ)
#define GreaterThani64(x,y) _mm256_cmpgt_epi64(x,y) 
#define GreaterThani32(x,y) _mm256_cmpgt_epi32(x,y)
#define GreaterThani16(x,y) _mm256_cmpgt_epi16(x,y)
#define GreaterThani8(x,y) _mm256_cmpgt_epi8(x,y)
#define LessThan(x,y) _mm256_cmp_ps(x,y,_CMP_LT_OQ)
#define LessThani(x,y) _mm256_cmpgt_epi32(y,x) 
#define LessThanOrEq(x,y) _mm256_cmp_ps(x,y,_CMP_LE_OQ)
#define LessThanOrEqi64(x,y) _mm256_cmple_epi64_mask(x,y)
#define LessThanOrEqi32(x,y) _mm256_cmple_epi32_mask(x,y)
#define LessThanOrEqi16(x,y) _mm256_cmple_epi16_mask(x,y)
#define LessThanOrEqi8(x,y) _mm256_cmple_epi8_mask(x,y)
#define GreaterThanOrEq(x,y) _mm256_cmp_ps(x,y,_CMP_GE_OQ)
#define GreaterThanOrEqi64(x,y) _mm256_cmpge_epi64_mask(x,y)
#define GreaterThanOrEqi32(x,y) _mm256_cmpge_epi32_mask(x,y)
#define GreaterThanOrEqi16(x,y) _mm256_cmpge_epi16_mask(x,y)
#define GreaterThanOrEqi8(x,y) _mm256_cmpge_epi8_mask(x,y)
#define NotEqual(x,y) _mm256_cmp_ps(x,y,_CMP_NEQ_OQ)
#define Floor(x) _mm256_floor_ps(x)
#define Max(x,y) _mm256_max_ps(x,y)
#define Maxi64(x,y) _mm256_max_epi64(x,y)
#define Maxi32(x,y) _mm256_max_epi32(x,y)
#define Maxi16(x,y) _mm256_max_epi16(x,y)
#define Maxi8(x,y) _mm256_max_epi8(x,y)
#define Min(x,y) _mm256_min_ps(x,y)
#define Mini64(x,y) _mm256_min_epi64(x,y)
#define Mini32(x,y) _mm256_min_epi32(x,y)
#define Mini16(x,y) _mm256_min_epi16(x,y)
#define Mini8(x,y) _mm256_min_epi8(x,y)
#define Gather(x,y,z) _mm256_i32gather_epi32(x,y,z)
#define Gatherf(x,y,z) _mm256_i32gather_ps(x,y,z)
#define BlendV(x,y,z) _mm256_blendv_ps(x,y,z)
#define BlendVi32(x,y,z) _mm256_blend_epi32(x,y,z)
#define BlendVi16(x,y,z) _mm256_blend_epi16(x,y,z)
#define BlendVi8(x,y,z) _mm256_blendv_epi8(x,y,z)
#endif // __AVX2__
