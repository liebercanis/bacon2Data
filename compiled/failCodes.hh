 enum FAILURECODES
{
  PASS = 0,
  BASEFAIL = 0x1,  // 2^0
  EARLYCUT = 0x2,  // 2^1
  FIRSTTIME = 0x4, // 2^2
  COSMIC = 0x8,    //// 2^3 cosmic event cut
  GAMMA = 0x10,    // // 2^4 overlap event cut
  TRIGFAIL = 0x20, // 2^5
  TRIGTIME = 0x40, // 2^6
  TRIANGLE = 0x80, // 2^7
  MULTI = 0x100,   // 2^8
  TOTALCODES = 2 * MULTI
};

enum
{
  FAILBITS = 11
};
