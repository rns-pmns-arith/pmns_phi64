//~ representations of polynomials Pi, used for conversion into the AMNS
//~ Note: each Pi is a representation of ((1 << RHO_LOG2)^i * phi^2)
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{0x40af4cbba1e4L, 0xb90406cdb762L, 0x3619e21bfb54L, 0x9a6015065f26L, 0x7a2a888a7e87L, 0x35f67d733837L, 0xb34cf3ee6121L, 0x373ff085007fL, 0xf77653e36511L, 0x3c39eb37a053L, 0x77dce6289cf1L, 0x73a298d736dcL, 0x932c8df7701fL, 0x2afbdd1287f6L, 0x2f9c323b699cL, 0x3ee061c2b108L, 0x7fb3017543dcL},
	{0xf179c6b4a11dL, 0x3307e5bd588fL, 0x4fb4774f63f7L, 0x50ef27d45ac5L, 0x10c27e0463243L, -0x4d0b2a34b763L, 0x79445b471b2fL, 0x82e265737e39L, 0x8800be5bfa8eL, 0xd2861b397bf8L, -0x5e3bba9fa55fL, 0x16326548ef20aL, -0x277ebb06d274L, 0x10b2f4c0feb17L, -0x3c6cf218dd07L, 0x9e63e3c0af20L, 0x3a92d894405eL},
	{0x96a568a31e89L, 0x6ed7ba5c8016L, -0x15caf4af26fcL, 0xbbd8c3bd75e4L, 0x7913e29bd2c9L, -0x1e2bdb976941L, 0x74a9411e1b76L, 0xad016d788c76L, 0x22b8a2c311daL, 0x1b79fcd841bL, 0x7f9636b37115L, 0xe25f524624c8L, -0x384c8b283077L, 0xa38c4504b939L, 0x6e4e85b10a82L, 0x5d81faa875f1L, 0x2337a56d7e0aL},
	{0xbf17d7e23fb4L, 0xefcc69755bc1L, 0x7a61daee2962L, -0xbef370bd0d2L, 0xe6cf46352432L, 0xa01a2ca6fea6L, -0x7cd520dd5d4L, 0xb95d6d5b4e4eL, 0x67628cea163eL, 0x148d87c074d97L, -0x4f2d35fbc966L, 0xdd1e6d37fcf6L, 0x33f29fb731feL, 0x6b50c447cc61L, 0x3c1bed89a88cL, 0x3d519c7cc6ddL, 0x47ef61f4b434L},
	{0x673cd4c03ea1L, 0xaca0cbd8e843L, 0xfe98b11eebL, 0x18e84bca4843L, 0x922afbac9725L, 0x17655a19ae72L, 0x3977da34614dL, 0x766a5dd35b3L, 0x157ac9494f7c2L, -0x67ccbba9805cL, 0x9a2106dbe831L, 0x216f6c65e5b4L, 0xd6e045878929L, 0x44a7a62122aL, 0x518a2c235598L, 0x861f6b544878L, 0x643d1a45cbebL},
	{0x12ba685d02b5cL, 0x1fdd349ded7L, 0x12a076ca1657aL, 0x2aa94339e262L, 0xe9fd2822211bL, -0xb9e41579c6aL, 0xae054e971080L, 0x5f3560e35443L, 0xabc05b36058aL, -0x27eec6d7fd13L, 0xde4f071865dbL, 0x30d9056bcb72L, 0x81474ccade42L, -0x5ab26287df1L, 0x7268cb446514L, 0x6d6ba87669e9L, 0x288a15d16d1cL},
	{0x11c9357d24cf4L, -0xcc14f6dd28fL, 0x83af9d5fa055L, 0x10dbd4ec6489L, 0x3e8d8403670aL, 0x6edde64925d8L, 0x8003602473f2L, 0xbc9d63521566L, -0x8e283a9c51d7L, 0xa9e2cd767a55L, 0x13b6508d5fa3L, 0x55818891c22cL, -0x4770b0ba1e52L, 0x9bf1de12018dL, 0x722e65d8d6adL, 0x4ce384b61059L, 0x31f366c7686L},
	{0x52317ffc96deL, -0x6c693cec7a88L, 0x140ea7eb7e0d4L, -0x1f4826cd3f62L, 0x7d5feb337397L, 0xd533e7514977L, 0x944a4f21efc5L, 0x1aa5f295c69eL, 0x4f243144583cL, 0x11889e48674ecL, 0x1183bc940a0aL, 0x24c9780cb1f7L, 0x1024b69be12ffL, 0x92220e220299L, -0x142b4bbc78d5L, 0x71f9542f1ffaL, 0x645112a2a063L},
	{0x17e40308f025dL, 0x7a2024f3312dL, 0x9a68824dc375L, 0x5475947ae17aL, -0x259023f44c05L, 0x56fb539722fL, -0x23fa348fd56aL, 0xa4186407a138L, 0x94d68b69044eL, 0x95781e4ff1f1L, 0x96ec3a4a915aL, 0x29780baa4553L, 0x3852736ac211L, 0x37886a28d7acL, 0x6d12d543260L, 0x42df49b00cbcL, 0x9bf77a88da4cL},
	{0x718122a1c372L, 0xea3240d081e5L, 0x98436fde2babL, 0x31f43ad41d73L, 0x3408bb6308d9L, 0x9b612314b3d5L, 0x9dabc6f7e7d5L, 0xba94404435b8L, -0x1ea38672da38L, 0xad696dd380cfL, 0x76bfb7271b9aL, 0xbc7de29d494cL, -0x50e1cbed8a19L, 0x888f8de7fb00L, 0x97f3de45ea0bL, 0x6ffbc6276e45L, 0x241d1959249cL},
	{0x1092aea33f816L, 0x58d9e7bf0532L, 0x8f7b08e57526L, 0x69cce9e6f639L, 0x550e073d43beL, 0xa933c97a3698L, 0x1c0b7af51a29L, 0xc571e387fe60L, -0x4c3ac243fbbL, 0xf07ae1c75285L, 0x2c9bca5b30e1L, 0xb77b0276864bL, 0x815b7c8d6cf4L, 0x9886df59d5ceL, 0x55f6fcb8b1c7L, 0x877f4be9a4ceL, 0x4d78fd5da147L},
	{0x10d52108c37ffL, 0xaf3237a16032L, 0x86d74c5b1e15L, 0x669c413156a1L, 0x63af273d1be7L, 0x366a120e808L, 0x12aea41531334L, -0x4d321d2b1a2L, 0x108eb2e025b74L, 0x9870b75cf4c0L, 0x567388dd0a37L, 0x592dffcc6b17L, -0x181de5b59414L, 0xb4c8cd5fdb93L, -0x18bc988b8a1L, 0x62865bb495cdL, 0x76682817738cL},
	{0x8a4c6569dd4eL, 0x9aa40c02f565L, -0x149082d8801bL, 0x875dc548839eL, 0x1aa547f2a795L, 0xac7e76c8f826L, 0x6cf7022c41a0L, 0x36ee9d587bb5L, 0x69be17da51b9L, 0x36e4455ae23eL, 0xd10888e8754eL, -0x5f5b84d76c5dL, 0x5b6578f3757bL, 0x61eb2a25d9a3L, 0x69523c5c761aL, 0x80cd82028c1L, 0x2d7650f421b0L},
	{0x27ba01c5df9cL, 0xdfafb445a66cL, 0x7567861fc8acL, 0x13a727bb74720L, 0x395e39a299caL, 0x157e9b37628fL, 0x89706ec56539L, 0xf12cf17e348dL, 0x625cc5997f78L, -0x3301864a0b5L, 0xc50d29584c8aL, 0xf20f26106786L, 0x43294ab26f80L, 0x484a93ad2f44L, 0x45271f910d43L, 0x86574eb086f0L, 0x501ef51178e3L},
	{0xb3980c0f0429L, 0x7ac9ddffee98L, 0xc9aae1db3640L, 0x901f9cc6af82L, 0x1780f34068f7L, 0x885006ed2fe3L, 0x695cb04a5454L, 0x7e336672bbafL, 0x6c2608b92d9dL, 0x9e4859122366L, 0xa1340ec59bb3L, 0x3eba54666740L, 0x9ffd2ba7eaa0L, 0x72b9656af016L, 0x40631e743901L, 0x68bf265d7e64L, 0x8dcc9a171da1L},
	{0x720761a625ddL, 0x223c37572cd2L, 0xf35e58535673L, 0x886c714c8963L, 0x9d4802209d25L, -0x4ed556dda4c4L, 0x99ff814b4fb5L, 0x3d93d6ae5509L, 0x6b249f1ee6f3L, 0x579c6aa27b98L, 0x40caf64134baL, 0xc2e0be91b774L, 0x4073eb83885fL, 0x9538c49c78c9L, -0x201e574f4251L, 0x6888578d1483L, 0x7872e9305591L},
	{0x2ddb9bc9681eL, 0x10aaeb4cf65b5L, -0x180e820a925L, 0x74c76ccd8f4fL, -0x1fa51e364696L, 0x824a644d22bdL, 0x68059d66ed84L, 0xb52c8af50149L, 0xb4eafed8ee8fL, 0x204ba531c5ebL, 0xaff754fa5e7fL, 0x4534f08276a3L, 0x96c179c3ce97L, -0x1a65dcc50eccL, 0xc1b9af3fb125L, 0x9a317614da6cL, 0x997b70343a61L}};
mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];