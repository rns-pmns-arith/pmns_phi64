//~ representations of polynomials Pi, used for conversion into the AMNS
//~ Note: each Pi is a representation of ((1 << RHO_LOG2)^i * phi^2)
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{0x2ac9ff7ed436L, 0x1479346bb18eL, -0x29d25e8e3559L, 0x34f644c09506L, 0x4ecabd7c9561L, -0x47dacbded1L, 0x12bd123d26a4L, 0x9f0cd13cf83L, 0x47aadabe047aL, 0x41a430bd043fL, 0xbcb4d9976adL, 0x12a0b59f62bcL, 0x487ffccb2ce0L, 0x88e3e799f0d8L, 0x196146eb209dL, 0x270f90acfa25L, 0x31610fd4aec9L, 0x3fd5c2ed5caL, 0x4b6c01590b8L, -0x19daad272efaL, 0x1c712d904bf5L, 0x231670c809f5L, -0x44e4b6618a23L, 0x2a8182a2353aL, 0x40aa105f97d7L, 0x3b4ab0f68daL, -0x17c34b506da4L, 0x2a22c2c7c795L, 0xb0a9a4b606dL, -0x3dff9d9e8becL, 0x34e97836b1f2L, -0x7050df7ee86L, -0x15cdbb42bb7L, 0x1bca75ec707eL, -0x2443049a1e78L},
	{0x1171615563L, -0x73d0551dd71L, -0x11eab20d35d0L, 0x2d570b93814dL, 0x1a3337a49df9L, -0x4117f1e2b94L, 0x3b822338c96dL, 0x1c4d512c60d9L, 0x438f43e15e54L, 0x4ae9ec66c922L, 0x3e99cab19a8dL, 0x5e4dd8df03b0L, 0x6508421fac6aL, 0x2d04a3e13e3bL, -0x19ddce78a0bcL, 0x168af20c6ce0L, 0x8adce9d4b035L, 0xc886f937f30L, -0x5c1888e382a4L, 0x6057108cad1dL, 0x1600481c111bL, -0x33edb5a89474L, 0x3e538444755dL, 0x11f9077dba2dL, -0x4e89efcfa3c0L, 0x3d97531d60a3L, 0x2196d963192L, -0x37d604df503bL, 0x5aca40a9df3aL, -0x9152aa290e9L, -0x69807d9ebcbcL, 0x1d6cac0f106bL, -0x284ffdd53d2bL, -0x92b6edb220dL, -0x800fd4c2133L},
	{0x34488b299264L, 0x121150ee7a63L, 0x3286ff96ac2L, 0x1f31d5a15dd9L, 0x3db02227c3e5L, 0x12c20037a42eL, 0x48c7daf7e3a3L, 0x506bde9fa0e2L, 0xcedfc4b9951L, 0x38dde327bdbbL, 0x23f482391e44L, 0x5fc3c879e2bL, 0x501b5869630eL, 0x3f97dc9b9283L, -0x2cd9dcecee90L, -0x1d8ab16847c2L, 0x162488709b15L, 0xf41d3efe89bL, -0x9c68934891fL, 0x6e922e85978L, -0xc95f215fdb9L, 0x293bf1fba464L, 0x36e4069b8c94L, 0x1837ffa5128fL, -0x2a1669855f64L, 0x2599cd82ec7aL, -0x45f5f18fafdL, 0x28f994c0ed70L, -0x34fd1805d5f5L, 0x1291eed07018L, 0x285b3dadac77L, 0x141421545656L, -0x9393aeb9e84L, -0x376000acf1d1L, 0x110f98a52f9eL},
	{0xc851ef6b457L, 0x194d31c9d4c5L, -0xe5964c85e09L, -0x1b691536547aL, 0x26c2c2e66d83L, -0xaa17b4c95b0L, 0x2b3016f73feL, -0x29bff1e15df2L, 0x285ec731b4eL, 0x23c3ebce454bL, 0x24dd036ed66dL, 0x360524fd2e89L, 0x3fe7e20ad0fdL, 0x2490a7c3452dL, -0x1067f5ee8d6aL, 0x773a415921d3L, 0x2ee26078f9abL, 0x7845fc4f92cL, -0xb5047f8c366L, 0x20afd1a62c31L, 0x6be31cbc6ebL, -0x15c9e43aa6b6L, 0x16ca0681152cL, -0x5becfc6223f9L, 0x14da66974c61L, 0x56ed1838ef1L, -0x1b90e39a4d53L, 0x6865b15ab150L, -0x186cb59a7d9bL, 0x7b0f6e13c568L, 0x16564bcd6047L, 0x13ec46fdf038L, 0x1d8c2fe1fbabL, 0x3f72522c98efL, -0x3cb87e8a7fecL},
	{-0x8c27dddd61L, -0x1dba999ee9a2L, -0x8c2c236b79dL, 0x4d0aa7b163eeL, 0x4cb76d992330L, 0xd221db894feL, 0x33dc8156a91eL, 0x68c746b31780L, 0x63bba5c9db28L, 0x39689babff8dL, 0x3ecabe7f2f06L, 0x4f9844b00fc5L, 0x27a1bb53a83aL, 0x49525646e2a9L, 0x48423971c648L, 0xb850f6b7950L, 0x3288a9c1e8c4L, 0x238c6819da38L, -0x994cea45547L, -0x219b2979db18L, -0x24a4a50595beL, -0x15d0fffa154bL, 0x1778493d9f14L, -0x1f57af679107L, -0x14f123ec298bL, -0xac11bb5dd6aL, 0x2cf5f45a7ff6L, 0x3969f0f37ceeL, -0x99320f2262dL, 0x132e7b0a93aaL, -0x466e37fbaa51L, 0x3be816cdbac8L, -0xa4e726ca2d7L, -0x1a3b65e572dfL, 0x109a76eeb663L},
	{-0x10b95976a1eaL, 0x5de1570f227L, 0x2b0eeca6f63L, 0x13b77644a38fL, 0x1d776a4116ffL, 0x1be099767b60L, 0x42435cef5ed3L, 0x2718b14490c7L, 0x115cebe604afL, 0x6feceb9dffffL, 0x597397be9614L, -0x3cf78e06545fL, 0x5401ffe870eL, 0x2313f6c37980L, -0xe31ebbceb3aL, -0x9714c63969dL, 0x167846e2f443L, -0x481399ad974L, 0x15ebde62aae4L, 0xfeb25569b6cL, 0x1ab61e20d9d7L, 0x314d37bd86ecL, -0xd428565503aL, 0x7c796d947f5L, -0x2f190832c908L, 0x2582a0c13599L, -0x225bf682a778L, 0x1ac24ba33654L, 0xd64d62e33a0L, -0x2d6e147be9d6L, 0x71d4cc56d342L, -0x33f5886de263L, 0x1122c2a3911eL, 0x334b7cd97e04L, -0x97fe3b9cfecL},
	{0x21576e830ea8L, -0x4394c45385ecL, -0x43f9341a4b1bL, 0x3ab436ced96aL, 0x3fc8f50a0f43L, 0x26f1e66bcb49L, 0x5a0edb255ac1L, 0x45e7e12d5500L, 0x33388d3323fbL, 0x6a2169380614L, 0x339bc420ab70L, 0x35c6a00e76dcL, 0x4a6111ca1ca9L, 0x4e9450e4b69bL, 0x19387d724f94L, 0x2bd61e063da7L, 0x5658b274c0eL, 0x46172a548202L, 0x2729bd3e3a7L, -0x3ca78502b896L, 0x248e5de95301L, 0x2cfc1caa638eL, -0x1957d8bd3f34L, -0x448dc8a92caL, 0x245bdd5837e1L, 0x314f07714e37L, -0x3286e1dd9584L, 0xd2d5e78e0e6L, 0x47356eadb66cL, -0x2e840b46fa2dL, 0xcf0c0d6208aL, -0xcbeff387d65L, 0x30a550666358L, -0x6b3f549d84edL, 0x2f3743ef073dL},
	{0x4e0a1dce48bL, 0x65a61efb9aeL, -0x3db5906e9e68L, -0x3410ab7e9148L, 0x408c8d436b48L, 0x6aeda261f04L, 0x2c91ad81a01L, 0x4c04eafd70afL, 0x366a548f4009L, -0x31d9aaedae90L, 0x63b58d9a9ec2L, 0x9f6fbaafadb2L, 0x11762cceb604L, 0x51b6f8f3078dL, 0x41d4c15b1551L, 0x12c992024cd8L, 0x3a6f77e99c46L, 0x36227b55db84L, -0x3b4828b3ce54L, 0x9c25bf1215cL, -0x97d4afdefd1L, 0x35c0bf38abdaL, 0xb4916739642L, 0x38a51069a5L, 0xf7639b1da45L, 0x32a180226446L, -0xfd19ba3c5e6L, -0x2b46b675c1a4L, 0x15c2aaaaaea7L, 0x5b986a6b20fbL, -0x219c7aa1ef99L, -0x19b9418e4b5bL, -0x388502a5b70L, 0x3b2a728c89eaL, -0x4417801f1362L},
	{0x56cc5048a1cbL, 0x59859624cc0dL, 0x607a7687e69L, 0x4a172ac61adfL, 0x4b14e50b0169L, -0xb4d9f3c4349L, 0x29f716b1d3f8L, 0x5ff10977425aL, -0x1ae90fd60205L, -0x1540acf71bf9L, 0x13f96a56c1bcL, 0xebd86a75656L, 0x16ac2b1df5fcL, -0x171db0b5903bL, 0x1c71b738921aL, 0x93364d6dc693L, -0x7b87a3ebb15L, -0x3a97e2cbb598L, 0x44cbff01eca2L, 0x28da37a7d6bbL, -0x15bbbd4f2591L, 0xa33c3ca966cL, -0x6d778b7ab07L, -0x3da992b9f858L, 0x1501d4dbe8edL, -0x18ab18bc5642L, -0x45afa0ef0d8eL, 0xffc66ae575dL, -0x25ed52ddba63L, -0x790e7b95cb1L, 0x5bcbe5285021L, -0x2703ece3bf0fL, 0x7afd24cb4049L, 0xe863802bdd8L, -0x1da5302aed13L},
	{-0x4718f651ca99L, -0xe7e2b39fa76L, 0x2a5c2cbef95eL, -0x530740f41b8eL, -0xbc06fd820aeL, 0x76b0a1582cfeL, 0x7715c363d095L, 0xa144b6f10faL, 0x19bb602b3073L, 0x5eff2898066dL, 0x26274f39ec3bL, 0x51104b298e01L, 0x19650f7db150L, 0x191a069df169L, 0x674d15968e33L, 0x1a7c151ff753L, -0x28e9880d23c9L, -0x25ac2a9ed75L, 0xc9da7698936L, -0x20b7e67dd50aL, 0x35ea69ec92eL, -0x35642f8caf4L, -0xfb3dd762dcaL, 0x2bc4b24f4bafL, 0xc8656a7f8efL, -0x3223cad49908L, -0x318358d875efL, 0x53287e9deb2aL, -0x2755c2a27a9aL, 0x21c9f4e465e9L, 0x45b1e3bc0a47L, -0xb22b0a2b808L, 0x23243d97b006L, 0x24f30c84c4e9L, -0x1a6cda552d96L},
	{-0x5703be1e03a4L, -0x9101d45f1ef8L, -0x726c7538b33L, -0x54fa6e73746L, -0x12d5b1f156edL, 0x31dc0f2dd00aL, 0x39615bdb2e08L, 0x2fba21f22b7eL, 0x442dc4ea4cd7L, 0x203b22fbcc7cL, 0x189ac4df997dL, 0x264cf719a6c9L, 0x183d9dea42a0L, 0x3ce76bcf2427L, -0x26e03096b952L, 0x4f49fa8a45e4L, 0x4f00b411a359L, 0x1746a3e10628L, 0x2b17baba57abL, 0x13d473ea5f61L, -0x2823b0d1349fL, -0x146dec95d21bL, 0x5904a42bfbdaL, -0x44492274ff52L, -0x2038636e9d42L, 0x20b0e98d7834L, 0x13b6d8d1aec7L, 0x5044e274d8a0L, -0x96f0f67ca94L, 0xb071e60986dL, 0x165300a8db60L, 0x3d45fc202c32L, -0x33f393ffb985L, -0x88211fa939dL, 0x29d30357795fL},
	{0x17cfb92d0fc4L, 0x3e48d591e7eeL, 0xdc329187440L, 0x39617d98c90bL, 0x2e6738eecfe7L, 0x5339ce8f481L, 0x75eddbe889d8L, 0x3aff1cee50d9L, 0x4d930c2b33caL, 0x42cd3a8913b8L, -0x3207c92a497L, 0x525cbfced128L, 0x4c09af8aebfdL, 0x406c0ca9be7bL, 0x12c4a15fd0faL, 0x138f2d3a5eL, 0x36769bdc15d3L, 0x1b9ca35bb48fL, 0xec839a989f7L, -0x39b8e2237defL, 0x8a6be97251eL, 0x26d4e33ec1ccL, -0x16c9491f4900L, 0x30582b4bee4fL, 0x1753242b5691L, -0x13ed733e472fL, -0x19a6b82b85dL, 0x3ac6764f1ea5L, -0x1f8bd4b3e7b8L, -0x43d803622e41L, 0x2333cad42fe8L, -0x12e3be321a27L, -0x12999dd85748L, 0x1289a451e583L, -0x32877daf172fL},
	{-0x52779ba4f318L, -0x14a3b2d333e9L, 0x19936ab7968aL, 0x40dab64726aL, 0x426c4f8d16b6L, 0x287dc5cb9b0dL, 0x720516898045L, 0x2eaa7e2718b7L, 0x5a3b8b04e55aL, 0x43a5b4af51ebL, 0x31b459a21ac3L, 0x47670d0d9831L, 0x3237a1d2e8d5L, 0x88e9d675704eL, -0x19ff6180d91bL, 0x338321cd4ad7L, 0x3a1c79c6de15L, 0x35a24bff0584L, 0xf8a1c960b6eL, -0xed05fd0f22fL, -0x3c17e2a05ab4L, 0x17056c380d87L, 0x13d80653ee93L, -0x153bd5b89172L, 0x10a369bb02f3L, 0x4fcef938911L, -0x14d91e31e728L, -0xc782c8e6045L, 0x1cd76f5b83dfL, 0xf748973c375L, -0x4a0974e2f397L, 0xb5419cbb959L, 0x756fae99d75L, 0x25e50ef443eeL, 0xd04782f6ff9L},
	{0x1d3a9634fb2fL, 0x239e336aea90L, 0x2d39589a2aceL, 0x5b59b50fe918L, 0x175a73abc5a6L, 0x59f2c31f24f8L, 0x35ec5a92cbe9L, -0x343f5f2b2e28L, 0x1f70bb51c197L, 0x274b313d8b88L, 0x29cb7e9dbee7L, -0x1173e43d60b1L, -0x28000ebc7392L, 0x211644dc5503L, 0x43110464f80bL, 0x159d8dee922dL, 0xa673a39423L, 0x6c27bb604955L, 0x120f83847bfdL, -0x209384898d3dL, 0x379711d6d779L, -0x3570693a676eL, -0x95664ef051aL, -0x8c7a9d8b0f4L, -0x47f5a6d111c5L, 0xbe5a0ac9114L, 0x2e471a27b7c8L, 0xef670bf8856L, 0x43acab54d3d1L, -0x4d7236e680L, 0x26e6efca5e06L, -0x996d8c866eeL, 0xe902ba4e7dfL, 0x128d7efafc01L, -0xac1cba0ff84L},
	{-0x18e801efd340L, -0x494c2d1e5befL, 0x118e1d309da5L, 0x34c0025e3313L, 0x176e8146854bL, 0x1886a039a00dL, 0x14fbcb271b69L, 0x355a1b450f50L, 0x788c5ae57e74L, 0x4e3f31aaac7eL, 0x2c072cd8a319L, 0x452545a2e1c8L, 0x361269c28554L, 0x7266e0c7f0cdL, 0x18b986d3c189L, 0x368b0ea7d1cfL, -0xcf17cf88a7cL, 0x50c8280e9be4L, 0x3a84dc051f9bL, -0x2858b184986eL, -0x10f65233ff9fL, -0x2032cf9cec0dL, 0x20aaa02ea4ccL, -0x2eed92a81f0dL, 0x1d72a62a17baL, 0x29d72015d00L, -0x57395ce0103cL, 0x16ad46ce35ccL, 0x2b2f5fb899fbL, -0xdb429a10d5eL, -0x3612a2c8d300L, -0x165bc64d7677L, 0x5afad2100fdL, 0x7d7c6638c12L, 0x257163d05347L},
	{0x228f0d27f4d7L, 0x9edcaa343edL, 0x20eff83042b3L, 0x49104cbaa3e4L, 0x2488963b039eL, 0x1ca91fa9b041L, 0x5c1937c3ea38L, 0x323c97367549L, 0x34155b8696ecL, 0x2ff7c7ad965cL, 0x2e6d1fbb8b63L, 0x21505ead86b9L, 0x5f7f004fba1L, 0x725e9506b17L, 0xa59ea867faaL, 0x5d669a6ec532L, -0x16e14f1e79daL, -0x693977eb151L, 0x5d03b5821fa6L, 0x6f9e7087e8f5L, -0x45c810710dc0L, 0x51b33a6ec85L, -0x4a8d3e5edcdL, -0x3caba52bd3aeL, -0xcea66ac3579L, -0xc052f52e4bdL, -0x25be54d3c42L, -0x12b0d900b1e7L, 0x26eaa757ca2fL, 0x370c769314a4L, 0x14d4909b0e36L, -0x1c35e2aa43a5L, -0x14bdc40e3656L, 0x14eeed144b89L, -0xd2f3076db3eL},
	{-0xa18b3204f43L, -0xec344d062d0L, -0x2977cc1f8d67L, 0x49a6e832098L, 0x47d26a0f9f26L, 0x24142a1c5530L, 0x2191f64a98a1L, 0x421583d842ceL, 0x4d26dbf13a3bL, 0xe5701fc6e22L, 0x3200194e9482L, 0x475f44970408L, 0x572df6eeb1e7L, 0x476b1a21fdb8L, 0x16b870f3b4f2L, 0x475f8967f1c4L, 0x1a84025fff6L, 0x2613aa8eae7cL, 0x37a3411946eL, -0x34778b969111L, -0x530275a8a43L, 0x4724ac2887ffL, -0x2e556474d9aL, 0x415e62c764dL, 0x399a683dc94cL, -0x1d6c09636826L, -0x325e7875b913L, -0x3740fe77c7feL, -0xb99a4a41617L, 0x53e5ba1e4effL, -0x43f6b58a68e4L, 0x3223879d13ecL, 0x5ce2963f2baL, 0x17efd360a2d2L, 0xb0eb2f21f2bL},
	{0xe7062fdf757L, 0x2f2d00b81d4eL, 0x195ce8c6d152L, -0x19d7d91beaa3L, 0x1d2c2a1bc640L, 0x3e1864458abaL, 0x17ff6a41e9c6L, 0x5ebc71ee86c7L, -0x14276645633fL, -0x4ba4795485b1L, 0x4a9a6d0e143bL, 0x2a9a6f54c0c9L, 0x11186f2a4747L, 0x35960f286cf9L, 0x221f834509daL, -0x1bb749e6443eL, 0x3bd8b6689dbbL, -0x23d3e2b93f3L, -0x2d7608b42682L, -0x18ab304aafacL, 0x2c64a18e9457L, 0x31dbe1cbfecaL, 0x16ad890ef665L, 0x2ceef3bd0a93L, 0xf9a00630689L, 0x2b37b827b41eL, -0x6482549eedaL, -0x32459a409a83L, -0x191b468abe05L, 0x2a2988eed5faL, -0x22ef802b2630L, 0x18f3b9184aeL, 0x4290398e4802L, -0x162cfd23c4e9L, 0xd020665d4c3L},
	{0x380d04bcfed5L, 0x2187d434f20cL, -0xd1931e416e4L, 0x6167b23d3a0cL, 0x46c949210a6cL, 0xe4b5c63e7d2L, 0x4ab65750502eL, 0x55768e166beeL, 0x4b6faa8b72d8L, 0x4b6872beb0f1L, 0x67a466bb00c9L, 0x4de44127cf71L, 0x7ae8b5d87c95L, 0x33e22a3da005L, -0x3ccece1f376dL, 0x6f6a26e8b4L, -0x9927835578eL, -0x2918587ba6c0L, -0xbba4d437f52L, -0x1bea221e27d8L, -0x74f0d3efc29L, -0x93e6ec580fL, 0xb3accbb2bbeL, 0x2270de3c6986L, 0x3ec0167ad18L, 0x239a7fc1edfeL, -0x21706c17fbc1L, 0x7224518ea723L, -0x15ec185b4f8cL, -0xbd07acba025L, 0x1fe9e877faaeL, -0x389b6e885f0aL, -0x149c13b8cc54L, -0x29a7906339f7L, -0x202360b8fe5fL},
	{-0x127fad08fee8L, -0x189af8e88e9dL, 0x25801a34446bL, 0x42d03d76ba3eL, 0x75b646000c5cL, 0xc72ac1213b3L, -0x2f8dd27fb814L, 0x535dbf62b9f5L, 0x141dd9fb77cbL, -0x62640aa1019L, 0x38a5f04b66b9L, 0x2af3b8228c9bL, 0x14ca9bf5025eL, 0x31941438991fL, 0x14d2eaf30dcbL, 0xf50d30d7b60L, 0x23541c119b8fL, -0x2a0f9894944eL, 0x19a68752acadL, 0x40090e542efcL, -0x4a41bfc5dedL, -0x2ddfca175bbL, 0x83f092cf132L, 0x2943241b74feL, 0x140602464fe4L, -0x3068ee3913f6L, -0x83db104bd02L, 0x6e908ebfbdc9L, 0x19f33c5a0e13L, 0x256d4828063eL, -0xb2e11eb360fL, -0x3378b5fc3dcaL, 0x2a20c2dcfd67L, -0xdcb555672e4L, -0xab1277bf088L},
	{0x29558bdf320aL, -0x186f6fb202a0L, -0x2cf72786b2bbL, 0x45a03cc694faL, 0x3d452ce18554L, 0x1f70022b3bfcL, 0x1d7bcb8f1e86L, 0x972e132047cbL, 0x3ff995b38dbdL, 0x2349058b251cL, 0x5f60fe8f8a0fL, 0x50e5ff7fa5c3L, 0x4ca4abd0704fL, -0x17113e7ff261L, 0x2d29f6c762e1L, -0x34d1cea735adL, 0x1996cc9f5609L, 0x1210b4f3a9a8L, -0x25d8a6303a9dL, 0x1d6e359004f3L, -0x6a3bfe6dbcL, 0x1ea71b8dc030L, 0x37eb94f9cc2fL, -0x2d6d2f180235L, 0xc193766db87L, 0x3445df2e80dL, -0x1c442c2d7c86L, -0xdc1e276d624L, -0x1693fd0f77f2L, 0x2f3eab03bcc8L, 0x28540c69b9dcL, -0x8e3adb99ea8L, 0xab7b8b43865L, 0xdac5a00cc01L, 0xa323e3827d6L},
	{-0x2b81150439d7L, -0x4e5c1a0626aL, 0x743676034f17L, 0xc50493484fbL, 0x386b4dd02fdL, 0x2036ff55d0bdL, 0x1c24814b57a3L, 0xb45f3c2934fL, -0x7759ba10bbaL, -0x28d1b9c56a2L, 0x45c156dfa697L, 0x428aff2a5842L, 0x1f05e3f198aL, 0x72b0c966fc0cL, 0x1a403a00c442L, -0x933691cfc9bL, -0xb738a96996aL, 0x106e95e39f30L, 0x6cfd8fc4fe8L, 0x50e627ce2bcL, -0x12a6f31b6180L, -0xf775e29ed0L, 0x5958344faa7fL, -0x2772569ba5f2L, -0x169df2b227ceL, 0x6fa0ef95991dL, -0x14528ef6ae94L, 0x3f4c87064e3L, -0x30f71e45cadL, 0x37ce4663f431L, 0x157ec5dde3b5L, 0x116757901d86L, -0x41f98682acc0L, 0x89424b69c74L, 0x1fa8c2a4f268L},
	{0x315e27c0e392L, 0x58898481de69L, 0x3e206cb39ffdL, 0x1df83533d432L, 0x3ffa4beed4e9L, 0x453545d7bdf0L, -0x1331f774e564L, 0x309e3a641640L, 0x295c3925033aL, -0x2684bf8b127dL, 0x35366288e797L, 0x8a0b904db7b2L, 0x2b71c23d2316L, 0x1f47ef816aa3L, 0x48f7c85774b0L, 0xd86b12147d2L, 0x86ffaad7cadeL, -0xd7a19ade7f8L, -0xba68ec2bcaL, 0x41b4d02f0c4L, -0x264f87e18f20L, -0x1e21113a5d3bL, -0x5149a2eed087L, -0x1c071499537cL, -0x1ad6082396f8L, 0x32b8430d1cf4L, -0x1cbf0bc592aL, 0x277e12e0bddaL, 0x68373f2688cL, 0x2db67d1523dfL, 0x16a3077e6e89L, -0x2020a1e1ca3cL, 0x185246d8f0c1L, -0x42ed403b128bL, 0x18c9a4d2b392L},
	{0x1c818f05b653L, -0x7f180ea22a6L, 0x36f852f924c5L, 0x484bc3ae2b12L, 0x1387f6c07941L, 0x2232213233a0L, 0x125f49013deeL, 0x65f281bcd28aL, 0xcbe7a204fe3L, -0x1ca4cbd777a0L, 0x8332046e8e35L, 0x3ed1a81e95f6L, 0x3c5e8756bdb1L, 0x301a7b822760L, 0x1642bfdbe557L, -0x18304bda8ff0L, 0x2133fe633541L, 0x1ddf44d5f29fL, -0x55d047fc06bL, 0x35d2e746a233L, -0xc88ac0d9e69L, 0x38ae206f8ba9L, 0x2e51ec603232L, 0x8a4b68301f4L, -0x180e901ebde6L, -0x72fe9096bc7L, -0x8bc7bf1dd67L, 0xcc990d181aeL, -0x19db7c09ee49L, -0xb98d1cf4087L, 0x295c93ea1fbfL, -0x8c6df9022f7L, -0xda00146922L, -0x25d7ee3d3e43L, -0x598e7a31005L},
	{-0xf6dd6895d7eL, -0x791b032718dL, -0x1a2dd813ea06L, -0x22c86c71c898L, 0x57a26b4faec8L, 0x45662c67c94aL, 0x387ee46f7b95L, 0x66ff4aa35238L, 0x6240dde358aeL, 0x355cf3cf24f4L, 0x256cbbd8d632L, 0x3afb9e831e54L, 0x20285066b736L, 0x2964a6a16213L, 0x14aedc630e41L, 0x3b06d638d8b3L, 0x1d0b5ae98c89L, 0x166773d8f6ebL, 0x2d9fdb7314bdL, -0x2c8315392438L, 0x2d75419cbca1L, 0x1791546bdb3cL, 0x57864b1c6b9L, -0xb370464ef35L, 0x20432b179dL, 0x9323c70f687L, -0x673dbf4d4e1bL, 0x3b30c81b9bf8L, -0x166674e8df94L, 0xad7b5aea25cL, 0x29b945758efbL, -0x1ecc699023e0L, 0x4d00801fbf93L, 0x3b5907f7a179L, -0x293b372dde60L},
	{0x55b70a88420L, 0x447822f35455L, -0xc569355d36L, -0x26bfc112012L, 0x469594372173L, 0x270af37f85ccL, 0x1c18946fde68L, 0x2cf85ea1584fL, 0x58a6d9a6a728L, 0x36f5e55dd3caL, 0x48971463587dL, 0x611f00e6d072L, 0xf344f39a378L, 0x2b178f94c5c8L, 0x45fda8900f27L, -0x250582dd0618L, 0xf3cf7f7e761L, 0x190f2a863147L, 0x23a24d1047a2L, 0x291f39de67a4L, 0xd3ca0bd802eL, -0x89295848bdeL, 0x3fff5f561aa7L, 0x2343054e5486L, -0x45052b954df4L, -0x1ef7873fd49bL, 0x2ebf43cff486L, -0x26b344dfc4c1L, 0x68a0214adaL, 0x12725d9e02b4L, -0x44353c6fc06L, -0x721756072dcL, -0x1963a5b8701eL, 0x5deb00220c76L, -0x2579197acdf4L},
	{0x870cfc27b9fL, -0x1c97f75e0bccL, -0x114e0c9a33bdL, 0x3c47ffacc6b1L, 0x1cdc27dc9c9dL, -0x4d98a25ca1L, 0x42d05f89583bL, -0x22836d8e2b80L, -0x6f78ba8f6d6L, 0x255b009e51ddL, 0x2fb404ba8c2bL, -0x1841c00826b0L, -0x25cbb0062bb1L, 0x76e28b7a3e32L, 0x5071cb7b0de3L, 0x2180962eec10L, 0x27f718bffcd4L, 0x4132e1590ddeL, -0x163862fbe46dL, -0x8fd5e11df83L, -0xd84f87539feL, 0x35a37b695928L, -0x35aba4df2c43L, 0x29433d0c2485L, 0x41d556772989L, -0x2846b1754976L, 0x4c3a3d035826L, -0x6ce8de978e8L, 0x12a52106b2aaL, -0x2cd140c30dc0L, -0xd607a6e2ad1L, 0x62984e6e9562L, -0x6a8721373e13L, 0x384481db3d10L, 0x1c3d9ea5ef6dL},
	{-0xa72cf41af90L, -0x1e540238aa9dL, 0x1f4a7acc0b25L, 0x64f6e2469416L, 0x57e01a5a2acaL, -0xfb2b113a1b5L, -0x2522816a5788L, 0x12266e713a47L, 0x5b20b23912aaL, 0x310b05d1c75L, -0x1dabff13297fL, 0x4a722f72f797L, 0x435f61c6b5ccL, 0x3eb4207c39c7L, 0x3a2809a38d4eL, 0x866c9e51d47cL, 0x11679bf5594eL, 0x1e87011b90c7L, 0x3ab6a27be3a3L, 0x2c741a6022c4L, 0xa10787bcf7L, -0x368de038025eL, 0x1277a99078d8L, -0x297a64f05fceL, -0x1e0e320c2e7fL, 0x177fa88008a4L, -0x15fdb3eef9eeL, 0x4754edd72a9cL, -0x27f6fc4d0204L, 0x383c46943b27L, 0x2361e880e76cL, -0x3337d3fe5604L, 0xb824865a6d9L, 0x1afb76ea712dL, -0x8c1b6c67abcL},
	{0x840bb8ce834L, 0x4f78a4cf4c4L, 0x93d56981372L, 0x325683e03407L, 0x24dafefbd99aL, 0x342ac415b50aL, 0x4d4d88f88f95L, 0x3b1572249a0aL, 0x2e0cfdc10093L, 0x929fc63927d8L, 0x424a0cd7665bL, -0x1dd01fbb9a44L, 0x58fee2fe4a9dL, 0x3905fc26417dL, 0x2b27376a10fbL, 0xefd5f5d6caaL, 0xfca8ac97a8eL, 0x1ea7244743bfL, 0x31f955ea3f66L, -0x6a490eb7524L, -0x24f09631a5cfL, 0x30e9a42876f8L, 0x220f001c8b90L, -0x3d9557de2573L, -0x24ebdc368df3L, -0x38b990fafb8L, 0x1c1c4331c5dcL, -0x22ee4b70da70L, -0x28e539fe73ecL, 0x68b90717b85L, 0x1425ff926a83L, 0xd84b36b7435L, -0x124cbe9a72baL, 0x65d56b80fb3L, 0x1b901dca92ebL},
	{-0xb30d7fc8feeL, 0x12e92b28954bL, 0x1f40826adffeL, 0x224d05c995c6L, 0x41c209173eb1L, 0x1cee6c00418L, 0x4f96d1b464bfL, 0x7fffde263a78L, 0x574dc31f2bcfL, 0x108c74894115L, 0xad54bbfeae2L, 0x3eb545da9873L, -0x54b2f32a2d3L, -0x12b746521ff9L, -0xbf986d74370L, 0x10964044db11L, 0xf7d655a3471L, 0x22f73125b973L, 0x2be90b374c1L, 0xe3940d29cf0L, 0x2260cb226d0dL, 0x22626e5f6316L, 0x121ccf2656fL, 0xcb2e598d8deL, -0x23d89807903dL, -0x5486dfae5baL, -0x2a06d7ff162eL, 0x36ad133fa691L, -0x56c7f8ed129bL, -0x16bbc5bc8381L, 0x3bb8d1ab7dcfL, -0x2cf8ca7fba88L, 0x2dd9d8cd31c7L, 0x2fc5d1101fcdL, -0x22862bf09940L},
	{0x2389a62371deL, 0xeb270de02d6L, 0x20db8f904c3bL, 0x298e127e1c47L, 0x14f820f4ed3aL, 0x1bbab487f1a1L, 0x1baa95b012a3L, 0x28c17fab56efL, 0x4134983f916dL, 0x489301b35581L, 0x1d23a0fc66eeL, 0x2e8caa07635bL, 0x63218fad1562L, 0x167eda4af224L, 0xacad5163fc9L, 0x69a1b911bdeL, -0x494cbdaa26cL, -0x30e1d5613234L, -0x539aaa5cdd33L, 0x20b5bd13e9fdL, -0x1e03248f18ffL, -0x4cd0fb8a920L, -0x2ec467ca9e19L, -0xcc86b79fb24L, -0x101afee6c01bL, 0x599b9fcc6819L, 0x25345f3d167dL, -0x238fc73e75f6L, 0x3320220a9927L, 0x15bfbae69f91L, 0x566aa9f7be40L, 0x7f428e06803L, -0x1144acf25abcL, -0xae6baf7dbfeL, -0xa36ad754dd6L},
	{-0x1c92a59ea866L, 0x298c78ebd5deL, 0x2601e9d9f508L, -0x3ddb4a674c5fL, 0x629fa1888955L, 0x4bc7c5d765dcL, 0x3c227bed1bfL, 0x700b8eb58cddL, 0x3c113298bed7L, 0x46fbbadf3e5L, 0x6ed9053859bdL, 0x28cb54170fcbL, -0x180dcf98a900L, 0x6fcfb4a90a25L, 0x3daa15ae33fL, -0x28240a398688L, 0x1ed782b412a5L, 0x7917af1e7b6L, -0x1d457d44b90aL, 0x4e6116ed2b2cL, -0x634c3b0224L, -0x15078ea74663L, 0x35d49e973ba0L, -0x8c3f8bd4569L, 0x39dd294f8adaL, 0x2d6edbde4805L, -0x26dd8606d534L, -0x799898eff92L, 0x3434ad457656L, -0xf229bf8bdb6L, 0x2a56ee12b11dL, -0x65995eb3795L, -0x1d18178cfb39L, 0xca1bde64da2L, 0x6fa63df4167L},
	{0x7e2ee70ed41L, -0x8475f3cb95aL, 0x5f432c7bca0L, 0x2ab9cad1805cL, -0x1a24366f36d6L, 0x1c49483765f3L, 0x4b5bb0b69cfaL, 0x702177e9313L, -0x17bd84a9314L, 0x64151abf2608L, -0x3355e7a1590L, 0x1cb7ac7eb3fcL, 0x5c7809f968daL, 0xcd7f633b8a1L, 0x1071ecc17538L, 0x27392131eb6fL, -0xd813d52291aL, 0x19992c9555edL, -0xaa3e838736bL, -0x8dfba10c069L, 0x37f9b7116174L, -0x145d829b95fcL, 0xde0fdbccfbbL, 0x2b389860ed13L, -0x1a2d3c813ed3L, -0x435c5dd463d5L, 0x54db2f73fe2bL, -0x194bf2d91f68L, -0x1864397f6a20L, 0x68c9613daddL, 0x244ed283c919L, 0x16b05e88569L, 0x5cdab7052ff9L, -0x2d3a25b06440L, 0x29fa7fc7e1b8L},
	{-0x4ce096e20096L, -0x2bcf3c6b40baL, -0x1c70565a1111L, -0x66e474429a16L, -0x22fd35a279b9L, 0x4cc7d8c49430L, 0x20b48f7c4a13L, -0x155b57618da0L, 0x5d211f127ce8L, 0x7c4476f9df9fL, 0x1308fb7e2ef4L, 0x5f29060e6c82L, 0x5a662f96ff48L, 0x4af99ed98d8dL, 0x3be5be79aeafL, 0x34174fdf8ed2L, 0x4464a8144bacL, 0x225fda37ec9cL, 0x14f432ed2727L, -0x1d0b266f6e4aL, 0x4336751de2afL, -0x415bc080c857L, 0x1e1111e9a603L, 0x951e480455cL, 0x2689baa42690L, 0x455c6bea0ee2L, -0x1e355d6ac4abL, 0x2c789043c9b8L, 0x2f5a18289e91L, 0x27f7d446fd9cL, -0x692f893f24L, -0x226845129818L, -0x21641c50c0dL, 0x7fda858eee1L, -0xd3e8556fe1eL},
	{-0x2e662656c792L, -0x8c178e9781fL, 0x384e886f54edL, 0x4b030681e1d8L, 0x49f48852c48dL, 0x26e06944a11cL, 0x2d74b9513f8eL, -0x1621822c0a5eL, -0x1d73f67bdb15L, 0x252b1767548bL, -0x34d96857c102L, -0x18995270c6f1L, -0x1303df16a9e8L, 0x813eef7d1fdL, 0x2238123f2a28L, -0x16190e3a0879L, 0x2a52e2061a57L, 0x38afb49f9eb5L, -0x1ca2ae749a98L, 0xbc9ef25c650L, 0x4563e127c3c8L, -0x19b3088de436L, 0x118309edfe75L, 0x259e79d46befL, -0x16b0d17f6b12L, -0x2c887a355c79L, 0x2b707a86d9deL, 0xde32ad1daedL, 0xce45aaaa407L, 0x2497eff8cc11L, 0x45f242a694aaL, 0xb13d5c0eb36L, 0x101d78ba2a6dL, 0x2ab9edd2a7ddL, 0x40767c92c653L}};
mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];
