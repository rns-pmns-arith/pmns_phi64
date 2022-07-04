//~ representations of polynomials Pi, used for conversion into the AMNS
//~ Note: each Pi is a representation of ((1 << RHO_LOG2)^i * phi^2)
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{0x62da9bfdcbb58L, -0x37bbf0cb33ad0L, -0x236f0fbe11de5L, -0x24e494c3e9fbdL, 0x2631baa420ddL, -0x42d9846bf878L, 0x5211b68f2c1cdL, 0x4277214c6f978L},
	{0x55c3f162a8d01L, -0x2836329fd11d4L, -0x24cd727cd070dL, -0x1aaf08b9e9a4bL, 0x2affa46010e64L, 0x1a27c5ec08a01L, 0x41c10d752ad13L, 0x4847d3e2d6e07L},
	{0x63dd9ebd3fbf9L, -0x398d50ddc59cfL, -0x6761c01f8a37fL, -0x582c3678ae0cfL, -0x2d05225035597L, -0x78fdb8ffb617L, 0x29d417ec5973fL, 0x5b6221d647f1dL},
	{0x59414adba2783L, -0x28fef2e7e418eL, -0x47481619e551fL, -0x57ed97a2dfd79L, -0x4a8ba39fcd5dL, -0x173f4fd6dae95L, 0x200ea6f049175L, 0x69052641398f0L},
	{0x733bb60eace14L, -0x48271ea536387L, -0x16793976037dbL, -0x46caedab08daaL, 0xa974be86d135L, 0x2e9e022654d78L, 0x1d102f8ba9a7bL, 0x7db1a847d3be3L},
	{0x84310e0722eabL, -0x64296f61bbf73L, -0x33638a5c35e0aL, -0x2e4fd9901ffe3L, 0x2318706a9a202L, 0x44069a715a25L, 0x614a1618cce2bL, 0x5f062343f6f8cL},
	{0x816b025437c92L, -0x4d6b129459b7fL, -0x4fd23fa870016L, -0x21105bd60329L, -0x14be39642a3f5L, 0x9bb2bc5e436eL, 0x53240761aa7aaL, 0x58c26f68c3f31L},
	{0x6600bd1de8c28L, -0x39f0abd80ae79L, -0x5b2760be842a1L, -0x2750d459c2a3eL, -0x1e6a8eb9bda53L, 0xc0de36e86562L, 0x4229cdd935ae6L, 0x5f0f2167455b7L}};
mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];
