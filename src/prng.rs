//! Project-owned deterministic PRNG primitives for future permutation work.

use crate::error::MageckError;

const SPLITMIX64_GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;
const SPLITMIX64_MUL_1: u64 = 0xbf58_476d_1ce4_e5b9;
const SPLITMIX64_MUL_2: u64 = 0x94d0_49bb_1331_11eb;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(SPLITMIX64_GAMMA);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(SPLITMIX64_MUL_1);
        z = (z ^ (z >> 27)).wrapping_mul(SPLITMIX64_MUL_2);
        z ^ (z >> 31)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Xoshiro256PlusPlus {
    state: [u64; 4],
}

impl Xoshiro256PlusPlus {
    pub(crate) fn from_seed(seed: u64) -> Self {
        let mut expansion = SplitMix64::new(seed);
        Self {
            state: [
                expansion.next_u64(),
                expansion.next_u64(),
                expansion.next_u64(),
                expansion.next_u64(),
            ],
        }
    }

    pub(crate) fn from_state(state: [u64; 4]) -> Result<Self, MageckError> {
        if state == [0; 4] {
            return Err(MageckError::Math(
                "xoshiro256++ initial state must not be all-zero".to_owned(),
            ));
        }
        Ok(Self { state })
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        let result = self.state[0]
            .wrapping_add(self.state[3])
            .rotate_left(23)
            .wrapping_add(self.state[0]);

        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);

        result
    }

    /// One ordinary null percentile from the upper 53 output bits.
    ///
    /// The exact power-of-two scale rounds the largest half-unit numerator to
    /// 1.0. The immediately smaller scale retains one draw and the half-unit
    /// construction while keeping every representable result in (0, 1).
    pub(crate) fn next_f64_open(&mut self) -> f64 {
        let bits53 = self.next_u64() >> 11;
        let exact_scale = 2.0_f64.powi(-53);
        let endpoint_safe_scale = f64::from_bits(exact_scale.to_bits() - 1);
        (bits53 as f64 + 0.5) * endpoint_safe_scale
    }

    /// Draws an unbiased pseudo-random integer in `[0, upper_exclusive)` using
    /// Lemire's (2019) multiply-shift rejection sampling algorithm.
    pub(crate) fn next_bounded(&mut self, upper_exclusive: u64) -> Result<u64, MageckError> {
        if upper_exclusive == 0 {
            return Err(MageckError::Config(
                "upper_exclusive bound must be positive (> 0)".to_owned(),
            ));
        }

        let mut product = (self.next_u64() as u128) * (upper_exclusive as u128);
        let mut low = product as u64;

        if low < upper_exclusive {
            let threshold = upper_exclusive.wrapping_neg() % upper_exclusive;
            while low < threshold {
                product = (self.next_u64() as u128) * (upper_exclusive as u128);
                low = product as u64;
            }
        }

        Ok((product >> 64) as u64)
    }
}

#[cfg(test)]
mod tests {
    use super::{SplitMix64, Xoshiro256PlusPlus};
    use crate::error::MageckError;
    use crate::stats::permutation_seed;

    type Transition = (u64, [u64; 4]);
    type TransitionCase = (&'static str, [u64; 4], &'static [Transition]);

    // Canonical JSON Part B: four SplitMix64 outputs and expanded state per seed.
    #[test]
    fn splitmix64_part_b_vectors() {
        let cases = [
            (
                0x0000000000000000,
                [
                    0xe220a8397b1dcdaf,
                    0x6e789e6aa1b965f4,
                    0x06c45d188009454f,
                    0xf88bb8a8724c81ec,
                ],
            ),
            (
                0x0000000000000001,
                [
                    0x910a2dec89025cc1,
                    0xbeeb8da1658eec67,
                    0xf893a2eefb32555e,
                    0x71c18690ee42c90b,
                ],
            ),
            (
                0x9e3779b97f4a7c15,
                [
                    0x6e789e6aa1b965f4,
                    0x06c45d188009454f,
                    0xf88bb8a8724c81ec,
                    0x1b39896a51a8749b,
                ],
            ),
            (
                0xffffffffffffffff,
                [
                    0xe4d971771b652c20,
                    0xe99ff867dbf682c9,
                    0x382ff84cb27281e9,
                    0x6d1db36ccba982d2,
                ],
            ),
            (
                0xf519f86ee2385b79,
                [
                    0x61cc42d4094b9c40,
                    0x0ef74aa7dd7f4aec,
                    0x0a45f286b9e7d382,
                    0x51ce3318362240f6,
                ],
            ),
            (
                0xf519f86ee2385b72,
                [
                    0x07dc930b302278a8,
                    0x12c015a97019e937,
                    0xcc06c31652ebf438,
                    0xecee65630a691e37,
                ],
            ),
            (
                0x0123456789abcdef,
                [
                    0x157a3807a48faa9d,
                    0xd573529b34a1d093,
                    0x2f90b72e996dccbe,
                    0xa2d419334c4667ec,
                ],
            ),
        ];

        for (seed, expected) in cases {
            let mut generator = SplitMix64::new(seed);
            for value in expected {
                assert_eq!(generator.next_u64(), value, "seed {seed:#018x}");
            }
        }
    }

    // Canonical JSON Part B: the four sequential SplitMix64 outputs populate slots 0..=3.
    #[test]
    fn xoshiro_seed_expansion_part_b_vectors() {
        let cases = [
            (
                0x0000000000000000,
                [
                    0xe220a8397b1dcdaf,
                    0x6e789e6aa1b965f4,
                    0x06c45d188009454f,
                    0xf88bb8a8724c81ec,
                ],
            ),
            (
                0x0000000000000001,
                [
                    0x910a2dec89025cc1,
                    0xbeeb8da1658eec67,
                    0xf893a2eefb32555e,
                    0x71c18690ee42c90b,
                ],
            ),
            (
                0x9e3779b97f4a7c15,
                [
                    0x6e789e6aa1b965f4,
                    0x06c45d188009454f,
                    0xf88bb8a8724c81ec,
                    0x1b39896a51a8749b,
                ],
            ),
            (
                0xffffffffffffffff,
                [
                    0xe4d971771b652c20,
                    0xe99ff867dbf682c9,
                    0x382ff84cb27281e9,
                    0x6d1db36ccba982d2,
                ],
            ),
            (
                0xf519f86ee2385b79,
                [
                    0x61cc42d4094b9c40,
                    0x0ef74aa7dd7f4aec,
                    0x0a45f286b9e7d382,
                    0x51ce3318362240f6,
                ],
            ),
            (
                0xf519f86ee2385b72,
                [
                    0x07dc930b302278a8,
                    0x12c015a97019e937,
                    0xcc06c31652ebf438,
                    0xecee65630a691e37,
                ],
            ),
            (
                0x0123456789abcdef,
                [
                    0x157a3807a48faa9d,
                    0xd573529b34a1d093,
                    0x2f90b72e996dccbe,
                    0xa2d419334c4667ec,
                ],
            ),
        ];

        for (seed, expected_state) in cases {
            assert_eq!(Xoshiro256PlusPlus::from_seed(seed).state, expected_state);
        }
    }

    // Canonical JSON Part C: every output and post-transition state, grouped by source case.
    #[test]
    fn xoshiro_part_c_transition_vectors() {
        let cases: &[TransitionCase] = &[
            (
                "seed_0",
                [
                    0xe220a8397b1dcdaf,
                    0x6e789e6aa1b965f4,
                    0x06c45d188009454f,
                    0xf88bb8a8724c81ec,
                ],
                &[
                    (
                        0x53175d61490b23df,
                        [
                            0x74d38efba8e829b7,
                            0x8a9c6b4b5aaded14,
                            0xd831b65330fc88e0,
                            0xbc8312de64d85a7e,
                        ],
                    ),
                    (
                        0x61da6f3dc380d507,
                        [
                            0x42ccf76e969d9edd,
                            0x267e53e3c2b94c43,
                            0x7a748df3423ca157,
                            0xb6ed46c3ef32a7ce,
                        ],
                    ),
                    (
                        0x5c0fdf91ec9a7bfc,
                        [
                            0xd25fe24ebb167550,
                            0x1ec6297e161873c9,
                            0x9f7fffef4c273f8a,
                            0x7d71b21262a405b1,
                        ],
                    ),
                    (
                        0x02eebf8c3bbe5e1a,
                        [
                            0xb1e87922cfaa0328,
                            0x53e634dfe1293913,
                            0x1fdc319110a34ada,
                            0x8ecf0c76f36d8e97,
                        ],
                    ),
                    (
                        0x7eca04ebaf4a5eea,
                        [
                            0x6cc1418bddeeb4ac,
                            0xfdd27c6c3e2070e1,
                            0xc78b8ae1ad2f49f2,
                            0x96f09ba527152248,
                        ],
                    ),
                    (
                        0x0543c37757f08d9a,
                        [
                            0x07e3a642c4dbe605,
                            0x5698b7064ee18dbf,
                            0x5392b72a9103fd5e,
                            0xaa552d645cf92326,
                        ],
                    ),
                    (
                        0xdb7490c75ab5026e,
                        [
                            0xfb2e3c20d6c3489c,
                            0x02e9a66e1b3996e4,
                            0x3a7d8cab4ea61b5b,
                            0x15d33f99b34c4243,
                        ],
                    ),
                    (
                        0xd87343e6464bc959,
                        [
                            0xec14a5d77eb69c3b,
                            0xc3ba16e5835cc523,
                            0x8d8f86f8b5ad53c7,
                            0xba94e2e7533ef50e,
                        ],
                    ),
                    (
                        0x4b7da0a02389f0ff,
                        [
                            0x953a51d5aed4ac16,
                            0xa22135ca48470adf,
                            0x4c502596415dcffc,
                            0x4605af25de805a0c,
                        ],
                    ),
                    (
                        0x1300fc58c0424c16,
                        [
                            0x711ecb3a3813fcc5,
                            0x7b4b4189a7ce6935,
                            0xb2fee4cdfa3763ea,
                            0xea1a7c84935df2d8,
                        ],
                    ),
                ],
            ),
            (
                "seed_p42_7",
                [
                    0x61cc42d4094b9c40,
                    0x0ef74aa7dd7f4aec,
                    0x0a45f286b9e7d382,
                    0x51ce3318362240f6,
                ],
                &[
                    (
                        0x57ebf9c2a4a5697a,
                        [
                            0x3ef53b6be216965a,
                            0x657efaf56dd3052e,
                            0xfec60aac25744fc2,
                            0xa1434be72f37fd6b,
                        ],
                    ),
                    (
                        0xe87de2b5c506b29d,
                        [
                            0xfac88a79a0f26e1f,
                            0xa54dcb32aab1dcb6,
                            0x35d9ea61cd3ed998,
                            0x9f08b887b622485c,
                        ],
                    ),
                    (
                        0x7b7414d4debf56c0,
                        [
                            0xc08df9ccbc61faf5,
                            0x6a5cab2ac77d6b31,
                            0x5974357bd4a0b787,
                            0x729d4748ae76a392,
                        ],
                    ),
                    (
                        0x4b43661bfffb9095,
                        [
                            0xd84c15aed56a3256,
                            0xf3a5679dafbc2643,
                            0xcfac424dbea34d72,
                            0x791463183d8c4d21,
                        ],
                    ),
                    (
                        0x3bd590ee9112e292,
                        [
                            0x52fd112b475a5934,
                            0xe445307ec4755967,
                            0xd8db089b274f7f24,
                            0x0d6c51562090b246,
                        ],
                    ),
                    (
                        0x93b106b1048a8de5,
                        [
                            0xbbd47003a3bfb215,
                            0x6e6329cea4607f77,
                            0xeadb915ad2db2610,
                            0xbd643d252c251c9c,
                        ],
                    ),
                    (
                        0x503c626afc7c4e6b,
                        [
                            0x68d364e82bfad1fe,
                            0x3f6cc897d504eb72,
                            0x0292a9998f8a9405,
                            0xac7d7a60e29d7108,
                        ],
                    ),
                    (
                        0x0d5ab109af057a6d,
                        [
                            0xfbc2d61f1c634b84,
                            0x552d05e67174ae89,
                            0xfb6e6778729445fb,
                            0x334f5262365ee6f3,
                        ],
                    ),
                    (
                        0x3c6c373857fad498,
                        [
                            0x9da0819b5b4903fe,
                            0x5581b4811f83a0f6,
                            0x0b60538e33e50e7f,
                            0x490f4ccc4af088e5,
                        ],
                    ),
                    (
                        0xd1739e61cd3c5be5,
                        [
                            0x812e79d60e3a2bed,
                            0xc3416694772fad77,
                            0xffc2ed1229400d81,
                            0x65026391df09aaae,
                        ],
                    ),
                ],
            ),
            (
                "manual_state_1_2_3_4",
                [1, 2, 3, 4],
                &[
                    (
                        0x0000000002800001,
                        [
                            0x0000000000000007,
                            0x0000000000000000,
                            0x0000000000040002,
                            0x0000c00000000000,
                        ],
                    ),
                    (
                        0x0000000003800067,
                        [
                            0x0000c00000000007,
                            0x0000000000040005,
                            0x0000000000040005,
                            0x0000000018000000,
                        ],
                    ),
                    (
                        0x000cc00003800067,
                        [
                            0x0000c00018040002,
                            0x0000c00000000007,
                            0x0000c008000e0002,
                            0x8000a00000000300,
                        ],
                    ),
                    (
                        0x000cc201994400b2,
                        [
                            0x8000a00018040305,
                            0x0000c008180a0007,
                            0x8000000818040000,
                            0x0060f0000c000000,
                        ],
                    ),
                    (
                        0x8012a2019ac433cd,
                        [
                            0x806090080c0e0302,
                            0x00006000180a0302,
                            0x8010901c000e0305,
                            0x4000e00c06010281,
                        ],
                    ),
                    (
                        0x8a69978acdee33ba,
                        [
                            0xc060100412050281,
                            0x00706014140a0305,
                            0xc07030000a040007,
                            0x60306800100183c1,
                        ],
                    ),
                    (
                        0xc271134733154abd,
                        [
                            0xa0201810160e8245,
                            0x006040100c0b0183,
                            0xc03808101e0b0286,
                            0x70188c0801028081,
                        ],
                    ),
                    (
                        0xac2ba09179169e97,
                        [
                            0xd058d4081b070347,
                            0x60785010040e8140,
                            0xe03808160b0380c3,
                            0x30204e0f198301a1,
                        ],
                    ),
                    (
                        0xdbf3190a8f073fd8,
                        [
                            0x8000ca17068a83a6,
                            0x50188c0e140a02c4,
                            0x9040d40312848384,
                            0xb01c2a0b03c3e3b1,
                        ],
                    ),
                    (
                        0x9105f14ab2229220,
                        [
                            0x60046c12114362d3,
                            0x4058921a000402e6,
                            0x085c360011860022,
                            0x3c2ebc0094c0a2f9,
                        ],
                    ),
                ],
            ),
            (
                "seed_0x0123456789abcdef",
                [
                    0x157a3807a48faa9d,
                    0xd573529b34a1d093,
                    0x2f90b72e996dccbe,
                    0xa2d419334c4667ec,
                ],
                &[
                    (
                        0xb2f2a310e96bd1c5,
                        [
                            0x62dd73afdc681de2,
                            0xef99ddb20943b6b0,
                            0x9fdce66a9cc46623,
                            0xf6efeef4e9750f1c,
                        ],
                    ),
                    (
                        0xb54062465b950493,
                        [
                            0x7bab40e93c5ea44e,
                            0x1298487749efcd71,
                            0x466587422dcc7bc1,
                            0xd735832ec668dc06,
                        ],
                    ),
                    (
                        0x87aca4a9668814b0,
                        [
                            0xbe068bb0b3d9b539,
                            0x2f568fdc587d12fe,
                            0xad2054748b70df8f,
                            0xe22ef8b5b96b31f0,
                        ],
                    ),
                    (
                        0xf13d2e2448a9cffb,
                        [
                            0x737efcd952cf9637,
                            0x3c70501860d47848,
                            0x0c9e6f3e1d556ab6,
                            0xc461d9af0eed3c22,
                        ],
                    ),
                    (
                        0xb7afdb427f6b86a2,
                        [
                            0x8b6f756e3cf6d25d,
                            0x4390c3ff2f4e84c9,
                            0xdfd0524fbf0afc81,
                            0x288d5f023136edc7,
                        ],
                    ),
                    (
                        0xc3a68c4e4f50d0c7,
                        [
                            0xe072e993228ebb53,
                            0x172fe4deacb2aa15,
                            0xd34179bc8a6e2edc,
                            0x0d21cd63b39fa3cf,
                        ],
                    ),
                    (
                        0x5bde00c2b40585ae,
                        [
                            0xfa7cc02e3da3b289,
                            0x241c74f104523f9a,
                            0xfa8ec94afcca958f,
                            0xa13b4341c537a3e5,
                        ],
                    ),
                    (
                        0xb27e2dd974f18e8a,
                        [
                            0x7f5bf79efcc62ef6,
                            0x24ee7d95c53b189c,
                            0xe91001c0be5d2706,
                            0xb38ff0a4e6f6182c,
                        ],
                    ),
                    (
                        0xa14dd5c28ddfa4ea,
                        [
                            0xe83a7aafdf0b2e46,
                            0xb2a58bcb87a0116c,
                            0x6d607c2873a309f0,
                            0xa01612ec31a62479,
                        ],
                    ),
                    (
                        0xb642d3593ecf568c,
                        [
                            0xfa89e388690d1b53,
                            0x37ff8d4c2b0836da,
                            0x92cd09c78e7027b6,
                            0xc6a2a2567324f6c0,
                        ],
                    ),
                    (
                        0xe9f7fc9172edb195,
                        [
                            0x0bd4cc923121db49,
                            0x5fbb6703cc750a3f,
                            0x72dcbc5f8ac93ce5,
                            0x98035e2ba5e34b05,
                        ],
                    ),
                    (
                        0x6ac04f255873c75e,
                        [
                            0xcc6cf5ba58b79a73,
                            0x26b317ce779ded93,
                            0xb70fe827af96e7ac,
                            0xc82758f707250d32,
                        ],
                    ),
                ],
            ),
        ];

        for (case_name, initial_state, expected_steps) in cases {
            let mut generator = Xoshiro256PlusPlus::from_state(*initial_state).unwrap();
            for (step, (expected_output, expected_state)) in expected_steps.iter().enumerate() {
                assert_eq!(
                    generator.next_u64(),
                    *expected_output,
                    "{case_name}, output step {step}"
                );
                assert_eq!(
                    generator.state, *expected_state,
                    "{case_name}, state step {step}"
                );
            }
        }
    }

    #[test]
    fn all_zero_state_is_rejected_and_seed_expansion_is_nonzero() {
        assert!(matches!(
            Xoshiro256PlusPlus::from_state([0; 4]),
            Err(MageckError::Math(message)) if message == "xoshiro256++ initial state must not be all-zero"
        ));
        assert!(Xoshiro256PlusPlus::from_state([1, 0, 0, 0]).is_ok());
        assert_ne!(Xoshiro256PlusPlus::from_seed(0).state, [0; 4]);
    }

    #[test]
    fn engines_are_reproducible_and_locally_owned() {
        let mut first = Xoshiro256PlusPlus::from_seed(0x0123_4567_89ab_cdef);
        let mut second = Xoshiro256PlusPlus::from_seed(0x0123_4567_89ab_cdef);
        let mut distinct = Xoshiro256PlusPlus::from_seed(0);
        for _ in 0..10 {
            let first_draw = first.next_u64();
            let _ = distinct.next_u64();
            assert_eq!(first_draw, second.next_u64());
        }
    }

    #[test]
    fn bounded_zero_bound_is_typed_error_without_state_change() {
        let mut generator = Xoshiro256PlusPlus::from_seed(0);
        let before = generator.state;

        assert!(matches!(
            generator.next_bounded(0),
            Err(MageckError::Config(message)) if message == "upper_exclusive bound must be positive (> 0)"
        ));
        assert_eq!(generator.state, before);
    }

    #[test]
    fn bounded_standard_vectors_and_bound_one_post_state() {
        let cases = [
            (1, 0),
            (2, 0),
            (3, 0),
            (4, 1),
            (5, 1),
            (7, 2),
            (10, 3),
            (16, 5),
            (17, 5),
            (100, 32),
            (4_294_967_295, 1_394_040_160),
            (9_223_372_036_854_775_809, 2_993_678_451_015_520_751),
            (u64::MAX, 5_987_356_902_031_041_502),
        ];
        let post_state = [
            0x74d3_8efb_a8e8_29b7,
            0x8a9c_6b4b_5aad_ed14,
            0xd831_b653_30fc_88e0,
            0xbc83_12de_64d8_5a7e,
        ];

        for (bound, expected) in cases {
            let mut generator = Xoshiro256PlusPlus::from_seed(0);
            assert_eq!(
                generator.next_bounded(bound).unwrap(),
                expected,
                "bound {bound}"
            );
            assert_eq!(generator.state, post_state, "bound {bound}");
        }
    }

    #[test]
    fn bounded_nonzero_threshold_without_rejection() {
        let mut generator =
            Xoshiro256PlusPlus::from_state([0, 1, 2, 0x2492_4a92_4924_9249]).unwrap();

        assert_eq!(generator.next_bounded(7).unwrap(), 2);
        assert_eq!(
            generator.state,
            [
                0x2492_4a92_4924_9248,
                0x0000_0000_0000_0003,
                0x0000_0000_0002_0002,
                0x9249_0492_4952_4924,
            ]
        );
    }

    #[test]
    fn bounded_single_rejection_vector() {
        let mut generator =
            Xoshiro256PlusPlus::from_state([0, 1, 2, 0xb6db_6edb_6db6_db6d]).unwrap();

        assert_eq!(generator.next_bounded(7).unwrap(), 4);
        assert_eq!(
            generator.state,
            [
                0x6db6_f800_006d_b6d9,
                0xb6db_6edb_6db4_db6d,
                0xb6db_6edb_6db2_db6e,
                0x6db6_bb6d_b2db_6dbb,
            ]
        );
    }

    #[test]
    fn bounded_four_rejection_trace_and_reproducibility() {
        let state = [
            0x42cc_f76e_969d_9edd,
            0x267e_53e3_c2b9_4c43,
            0x7a74_8df3_423c_a157,
            0xb6ed_46c3_ef32_a7ce,
        ];
        let bound = 0x8000_0000_0000_0001;
        let post_state = [
            0xfb2e_3c20_d6c3_489c,
            0x02e9_a66e_1b39_96e4,
            0x3a7d_8cab_4ea6_1b5b,
            0x15d3_3f99_b34c_4243,
        ];
        let mut first = Xoshiro256PlusPlus::from_state(state).unwrap();
        let mut second = Xoshiro256PlusPlus::from_state(state).unwrap();

        assert_eq!(first.next_bounded(bound).unwrap(), 0x6dba_4863_ad5a_8137);
        assert_eq!(first.state, post_state);
        assert_eq!(second.next_bounded(bound).unwrap(), 0x6dba_4863_ad5a_8137);
        assert_eq!(second.state, post_state);
    }

    #[test]
    fn bounded_outputs_are_strictly_below_bound() {
        let bounds = [1, 2, 3, 7, 10, 17, 100, u64::MAX];
        for bound in bounds {
            let mut generator = Xoshiro256PlusPlus::from_seed(0x0123_4567_89ab_cdef);
            for _ in 0..64 {
                assert!(generator.next_bounded(bound).unwrap() < bound);
            }
        }
    }

    #[test]
    fn existing_permutation_seed_anchors_are_preserved() {
        let cases = [
            (0, 0, 0),
            (0, 1, 1),
            (1, 0, 0x9e3779b97f4a7c15),
            (1, 1, 0x9e3779b97f4a7c16),
            (0x3039, 0x1a85, 0xa12ce22b4ed9ab32),
            (u64::MAX, 0, 0x61c8864680b583eb),
            (u64::MAX, 1, 0x61c8864680b583ec),
            (0x8000000000000000, 0x8000000000000000, 0),
            (u64::MAX, u64::MAX, 0x61c8864680b583ea),
            (0x2a, 0, 0xf519f86ee2385b72),
        ];
        for (base, task, expected) in cases {
            assert_eq!(permutation_seed(base, task), expected);
        }
    }

    #[test]
    fn p4f6_open_draw_excludes_both_endpoints_and_advances_once() {
        for (state, expected_bits) in [([0, 1, 0, 0], 0_u64), ([0, 1, 0, u64::MAX], u64::MAX)] {
            let mut open = Xoshiro256PlusPlus::from_state(state).unwrap();
            let mut integer = open.clone();
            assert_eq!(integer.next_u64(), expected_bits);
            let value = open.next_f64_open();
            assert!(value.is_finite() && 0.0 < value && value < 1.0);
            assert_eq!(open, integer);
        }
    }

    #[test]
    fn p4f6_open_draw_is_repeatable_for_seed_vectors() {
        for seed in [0, 1, 0x0123_4567_89ab_cdef, u64::MAX] {
            let mut first = Xoshiro256PlusPlus::from_seed(seed);
            let mut second = Xoshiro256PlusPlus::from_seed(seed);
            for _ in 0..256 {
                let value = first.next_f64_open();
                assert!(value.is_finite() && 0.0 < value && value < 1.0);
                assert_eq!(value.to_bits(), second.next_f64_open().to_bits());
            }
        }
    }

    #[test]
    fn p4f6_open_draw_exact_seed_vectors() {
        let cases = [
            (0, [0x3fd4_c5d7_5852_42c8, 0x3fd8_769b_cf70_e034]),
            (1, [0x3fe9_f8ba_0fed_e077, 0x3fe7_e848_2652_c7fb]),
            (
                0x0123_4567_89ab_cdef,
                [0x3fe6_5e54_621d_2d79, 0x3fe6_a80c_48cb_729f],
            ),
        ];
        for (seed, expected) in cases {
            let mut rng = Xoshiro256PlusPlus::from_seed(seed);
            for bits in expected {
                assert_eq!(rng.next_f64_open().to_bits(), bits, "seed {seed:#x}");
            }
        }
    }
}
