#include "mechanism.H"
const int rmap[235] = {
  10,  22,  35,  82,  121, 133, 139, 168, 224, 225, 155, 0,   4,   6,   7,
  8,   12,  13,  95,  96,  105, 138, 144, 149, 1,   2,   3,   5,   9,   11,
  14,  15,  16,  17,  18,  19,  20,  21,  23,  24,  25,  26,  27,  28,  29,
  30,  31,  32,  33,  34,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,
  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
  76,  77,  78,  79,  80,  81,  83,  84,  85,  86,  87,  88,  89,  90,  91,
  92,  93,  94,  97,  98,  99,  100, 101, 102, 103, 104, 106, 107, 108, 109,
  110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125,
  126, 127, 128, 129, 130, 131, 132, 134, 135, 136, 137, 140, 141, 142, 143,
  145, 146, 147, 148, 150, 151, 152, 153, 154, 156, 157, 158, 159, 160, 161,
  162, 163, 164, 165, 166, 167, 169, 170, 171, 172, 173, 174, 175, 176, 177,
  178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,
  193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
  208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222,
  223, 226, 227, 228, 229, 230, 231, 232, 233, 234};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 235; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[235] = {
    2, 4, 4, 4, 2, 3, 3, 3, 2, 4, 3, 3, 3, 3, 4, 4, 4, 3, 4, 2, 4, 4, 2, 4,
    4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
    2, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
    4, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 5, 3, 3, 4, 4, 4, 4,
    3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 3, 4, 4, 4,
    2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 3, 3, 2, 4, 4, 4, 4, 4, 5, 5, 4, 4, 4, 2, 4, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 4, 4, 3};
  const int kiv[1175] = {
    1,  2,  0,  0,  0, 1,  4,  2,  6,  0, 1,  4,  2,  6,  0,  1,  6,  2,  5,  0,
    2,  1,  0,  0,  0, 2,  3,  6,  0,  0, 2,  4,  6,  0,  0,  2,  6,  5,  0,  0,
    4,  3,  0,  0,  0, 2,  3,  4,  6,  0, 2,  6,  5,  0,  0,  6,  5,  4,  0,  0,
    2,  4,  6,  0,  0, 2,  4,  7,  0,  0, 5,  7,  5,  6,  0,  1,  7,  1,  6,  0,
    0,  7,  0,  6,  0, 6,  7,  6,  0,  0, 2,  7,  2,  6,  0,  7,  6,  0,  0,  0,
    3,  7,  3,  6,  0, 4,  7,  4,  6,  0, 8,  6,  0,  0,  0,  2,  8,  5,  6,  0,
    2,  8,  1,  9,  0, 8,  4,  9,  6,  0, 8,  6,  5,  9,  0,  8,  6,  5,  9,  0,
    2,  9,  1,  3,  0, 2,  9,  6,  0,  0, 9,  4,  3,  6,  0,  9,  6,  5,  3,  0,
    9,  6,  5,  3,  0, 9,  8,  3,  0,  0, 9,  8,  3,  0,  0,  2,  3,  9,  0,  0,
    2,  11, 1,  12, 0, 11, 4,  12, 6,  0, 9,  11, 8,  12, 0,  1,  13, 2,  12, 0,
    12, 4,  2,  23, 0, 12, 4,  13, 6,  0, 12, 4,  13, 6,  0,  12, 6,  5,  13, 0,
    9,  12, 5,  23, 0, 9,  12, 5,  23, 0, 9,  12, 5,  24, 0,  12, 3,  21, 4,  0,
    12, 3,  23, 6,  0, 12, 13, 11, 0,  0, 13, 12, 14, 11, 0,  14, 12, 2,  0,  0,
    23, 12, 11, 10, 0, 26, 12, 11, 25, 0, 12, 25, 21, 10, 0,  12, 25, 5,  30, 0,
    2,  13, 1,  14, 0, 13, 4,  2,  10, 0, 13, 6,  2,  23, 0,  13, 6,  5,  14, 0,
    13, 3,  23, 4,  0, 13, 3,  10, 6,  0, 13, 14, 12, 0,  0,  14, 13, 2,  0,  0,
    13, 10, 2,  30, 0, 13, 10, 0,  6,  0, 26, 13, 12, 25, 0,  13, 25, 30, 6,  0,
    13, 25, 23, 10, 0, 14, 6,  2,  10, 0, 14, 3,  10, 4,  0,  14, 10, 0,  4,  0,
    15, 2,  0,  0,  0, 2,  15, 1,  0,  0, 15, 4,  2,  30, 0,  15, 4,  0,  6,  0,
    15, 4,  13, 10, 0, 15, 6,  5,  0,  0, 15, 3,  9,  0,  0,  13, 15, 0,  12, 0,
    12, 15, 0,  11, 0, 15, 10, 23, 0,  0, 20, 12, 6,  0,  0,  2,  20, 1,  22, 0,
    2,  20, 1,  21, 0, 20, 4,  22, 6,  0, 20, 4,  21, 6,  0,  20, 6,  5,  22, 0,
    20, 6,  21, 5,  0, 12, 20, 22, 11, 0, 12, 20, 21, 11, 0,  13, 20, 22, 12, 0,
    13, 20, 21, 12, 0, 9,  20, 8,  22, 0, 9,  20, 21, 8,  0,  21, 2,  23, 0,  0,
    21, 22, 0,  0,  0, 2,  21, 1,  23, 0, 2,  21, 12, 6,  0,  21, 4,  23, 6,  0,
    21, 6,  5,  23, 0, 21, 9,  8,  23, 0, 21, 12, 23, 11, 0,  21, 10, 23, 0,  0,
    21, 25, 23, 26, 0, 22, 2,  23, 0,  0, 2,  22, 12, 6,  0,  2,  22, 1,  23, 0,
    22, 4,  23, 6,  0, 22, 4,  23, 6,  0, 22, 6,  5,  23, 0,  22, 9,  8,  23, 0,
    22, 9,  20, 3,  0, 22, 3,  23, 9,  0, 22, 12, 23, 11, 0,  22, 25, 23, 26, 0,
    2,  23, 1,  10, 0, 23, 4,  10, 6,  0, 23, 9,  27, 6,  0,  23, 5,  30, 0,  0,
    9,  10, 25, 6,  0, 10, 4,  25, 0,  0, 25, 4,  10, 3,  0,  9,  25, 26, 3,  0,
    9,  25, 27, 3,  0, 25, 10, 3,  0,  0, 25, 10, 28, 0,  0,  26, 4,  25, 6,  0,
    26, 6,  5,  25, 0, 26, 25, 29, 10, 0, 26, 5,  10, 25, 0,  27, 4,  25, 6,  0,
    27, 6,  5,  25, 0, 25, 4,  28, 0,  0, 2,  28, 25, 6,  0,  28, 4,  25, 3,  0,
    28, 6,  9,  25, 0, 9,  28, 25, 3,  6, 28, 10, 3,  0,  0,  25, 6,  29, 0,  0,
    2,  29, 1,  28, 0, 2,  29, 5,  25, 0, 2,  29, 26, 6,  0,  29, 6,  5,  28, 0,
    9,  10, 31, 0,  0, 2,  31, 1,  28, 0, 2,  31, 5,  25, 0,  2,  31, 26, 6,  0,
    31, 6,  5,  28, 0, 24, 2,  10, 0,  0, 2,  24, 2,  23, 0,  2,  24, 13, 6,  0,
    24, 4,  10, 6,  0, 24, 6,  2,  26, 0, 24, 3,  25, 6,  0,  30, 0,  4,  0,  0,
    2,  30, 0,  6,  0, 2,  30, 0,  7,  0, 30, 4,  0,  3,  0,  30, 4,  10, 0,  0,
    30, 6,  9,  0,  0, 30, 6,  23, 10, 0, 30, 10, 0,  25, 0,  12, 1,  18, 0,  0,
    12, 1,  19, 0,  0, 13, 12, 2,  18, 0, 22, 12, 17, 6,  0,  22, 12, 19, 5,  0,
    12, 16, 0,  0,  0, 16, 1,  19, 0,  0, 2,  16, 1,  17, 0,  2,  16, 12, 11, 0,
    16, 4,  17, 6,  0, 16, 4,  5,  18, 0, 16, 6,  5,  17, 0,  16, 12, 17, 11, 0,
    16, 10, 23, 17, 0, 16, 25, 26, 17, 0, 16, 25, 27, 17, 0,  2,  17, 1,  18, 0,
    17, 4,  18, 6,  0, 17, 4,  23, 12, 0, 17, 4,  2,  12, 10, 17, 6,  5,  18, 0,
    17, 6,  19, 5,  0, 17, 6,  23, 11, 0, 9,  17, 8,  18, 0,  9,  17, 16, 3,  0,
    17, 12, 18, 11, 0, 17, 12, 19, 11, 0, 17, 13, 18, 12, 0,  2,  18, 1,  15, 0,
    18, 4,  15, 6,  0, 18, 4,  12, 10, 0, 18, 6,  5,  15, 0,  18, 12, 11, 15, 0,
    18, 13, 12, 15, 0, 18, 10, 30, 12, 0, 19, 2,  15, 0,  0,  19, 1,  0,  0,  0,
    19, 18, 0,  0,  0, 2,  19, 1,  15, 0, 2,  19, 2,  18, 0,  19, 4,  15, 6,  0,
    19, 4,  12, 10, 0, 19, 6,  5,  15, 0, 19, 6,  2,  12, 10, 19, 9,  12, 10, 6,
    19, 9,  8,  15, 0, 19, 3,  12, 25, 0, 19, 12, 11, 15, 0,  27, 26, 0,  0,  0,
    2,  25, 10, 6,  0, 10, 6,  26, 0,  0, 10, 6,  27, 0,  0,  2,  26, 1,  25, 0,
    2,  27, 1,  25, 0, 2,  26, 5,  10, 0, 23, 6,  2,  26, 0,  23, 6,  5,  10, 0,
    1,  30, 5,  0,  0, 11, 6,  5,  12, 0, 2,  27, 5,  10, 0,  2,  27, 23, 6,  0,
    2,  12, 11, 0,  0, 2,  10, 23, 0,  0, 9,  12, 21, 6,  0,  12, 10, 5,  0,  0,
    12, 10, 15, 6,  0, 21, 3,  23, 9,  0, 9,  12, 11, 3,  0,  18, 2,  15, 0,  0,
    23, 25, 26, 10, 0, 23, 3,  9,  10, 0, 17, 2,  18, 0,  0};
  const int nuv[1175] = {
    -1, 2,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  0, 0, 0, -2, -1, 2, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -2, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 2,  0, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 2,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 235) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[32]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 32; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 235; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[32];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 14.007000; // N
  awt[1] = 1.008000;  // H
  awt[2] = 15.999000; // O
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 3;
  // Zero ncf
  for (int id = 0; id < kd * 32; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 0] = 2; // N

  // H2
  ncf[1 * kd + 1] = 2; // H

  // H
  ncf[2 * kd + 1] = 1; // H

  // O2
  ncf[3 * kd + 2] = 2; // O

  // O
  ncf[4 * kd + 2] = 1; // O

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 2] = 1; // O

  // OH
  ncf[6 * kd + 1] = 1; // H
  ncf[6 * kd + 2] = 1; // O

  // OHV
  ncf[7 * kd + 1] = 1; // H
  ncf[7 * kd + 2] = 1; // O

  // H2O2
  ncf[8 * kd + 1] = 2; // H
  ncf[8 * kd + 2] = 2; // O

  // HO2
  ncf[9 * kd + 1] = 1; // H
  ncf[9 * kd + 2] = 2; // O

  // NO
  ncf[10 * kd + 0] = 1; // N
  ncf[10 * kd + 2] = 1; // O

  // NH3
  ncf[11 * kd + 1] = 3; // H
  ncf[11 * kd + 0] = 1; // N

  // NH2
  ncf[12 * kd + 1] = 2; // H
  ncf[12 * kd + 0] = 1; // N

  // NH
  ncf[13 * kd + 1] = 1; // H
  ncf[13 * kd + 0] = 1; // N

  // N
  ncf[14 * kd + 0] = 1; // N

  // NNH
  ncf[15 * kd + 1] = 1; // H
  ncf[15 * kd + 0] = 2; // N

  // N2H4
  ncf[16 * kd + 1] = 4; // H
  ncf[16 * kd + 0] = 2; // N

  // N2H3
  ncf[17 * kd + 1] = 3; // H
  ncf[17 * kd + 0] = 2; // N

  // N2H2
  ncf[18 * kd + 1] = 2; // H
  ncf[18 * kd + 0] = 2; // N

  // H2NN
  ncf[19 * kd + 1] = 2; // H
  ncf[19 * kd + 0] = 2; // N

  // NH2OH
  ncf[20 * kd + 1] = 3; // H
  ncf[20 * kd + 0] = 1; // N
  ncf[20 * kd + 2] = 1; // O

  // H2NO
  ncf[21 * kd + 1] = 2; // H
  ncf[21 * kd + 0] = 1; // N
  ncf[21 * kd + 2] = 1; // O

  // HNOH
  ncf[22 * kd + 1] = 2; // H
  ncf[22 * kd + 0] = 1; // N
  ncf[22 * kd + 2] = 1; // O

  // HNO
  ncf[23 * kd + 1] = 1; // H
  ncf[23 * kd + 0] = 1; // N
  ncf[23 * kd + 2] = 1; // O

  // HON
  ncf[24 * kd + 1] = 1; // H
  ncf[24 * kd + 0] = 1; // N
  ncf[24 * kd + 2] = 1; // O

  // NO2
  ncf[25 * kd + 0] = 1; // N
  ncf[25 * kd + 2] = 2; // O

  // HONO
  ncf[26 * kd + 1] = 1; // H
  ncf[26 * kd + 0] = 1; // N
  ncf[26 * kd + 2] = 2; // O

  // HNO2
  ncf[27 * kd + 1] = 1; // H
  ncf[27 * kd + 0] = 1; // N
  ncf[27 * kd + 2] = 2; // O

  // NO3
  ncf[28 * kd + 0] = 1; // N
  ncf[28 * kd + 2] = 3; // O

  // HONO2
  ncf[29 * kd + 1] = 1; // H
  ncf[29 * kd + 0] = 1; // N
  ncf[29 * kd + 2] = 3; // O

  // N2O
  ncf[30 * kd + 0] = 2; // N
  ncf[30 * kd + 2] = 1; // O

  // HNO3
  ncf[31 * kd + 1] = 1; // H
  ncf[31 * kd + 0] = 1; // N
  ncf[31 * kd + 2] = 3; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(3);
  ename[0] = "N";
  ename[1] = "H";
  ename[2] = "O";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(32);
  kname[0] = "N2";
  kname[1] = "H2";
  kname[2] = "H";
  kname[3] = "O2";
  kname[4] = "O";
  kname[5] = "H2O";
  kname[6] = "OH";
  kname[7] = "OHV";
  kname[8] = "H2O2";
  kname[9] = "HO2";
  kname[10] = "NO";
  kname[11] = "NH3";
  kname[12] = "NH2";
  kname[13] = "NH";
  kname[14] = "N";
  kname[15] = "NNH";
  kname[16] = "N2H4";
  kname[17] = "N2H3";
  kname[18] = "N2H2";
  kname[19] = "H2NN";
  kname[20] = "NH2OH";
  kname[21] = "H2NO";
  kname[22] = "HNOH";
  kname[23] = "HNO";
  kname[24] = "HON";
  kname[25] = "NO2";
  kname[26] = "HONO";
  kname[27] = "HNO2";
  kname[28] = "NO3";
  kname[29] = "HONO2";
  kname[30] = "N2O";
  kname[31] = "HNO3";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (Jac[33 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[33 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[33 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 33;
    int offset_col = nc * 33;
    for (int k = 0; k < 33; k++) {
      for (int l = 0; l < 33; l++) {
        if (Jac[33 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (Jac[33 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (Jac[33 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[33 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[33 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 33 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[33 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 33 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 33; l++) {
      for (int k = 0; k < 33; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[33 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 33; l++) {
      for (int k = 0; k < 33; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[33 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
