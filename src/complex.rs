
macro_rules! rect_common_impls {
    ($t:ty) => {
        impl Rect {
            #[inline(always)]
            pub fn zero() -> Rect { 0.0.into() }
            #[inline(always)]
            pub fn one() -> Rect { 1.0.into() }
            #[inline(always)]
            pub fn i() -> Rect { Rect { real: 0.0, imag: 1.0 } }

            #[inline]
            pub fn from_phase(radians: $t) -> Rect {
                Rect {
                    real: radians.cos(),
                    imag: radians.sin(),
                }
            }

            #[inline(always)]
            pub fn sqnorm(self) -> $t { self.real * self.real + self.imag * self.imag }
            #[inline]
            pub fn abs(self) -> $t { <$t>::hypot(self.real, self.imag) }
            #[inline(always)] // (otherwise it doesn't get inlined into `Ket::overlap`)
            pub fn conj(self) -> Rect {
                Rect {
                    real:  self.real,
                    imag: -self.imag,
                }
            }

            #[inline]
            pub fn lexical_cmp(self, other: Rect) -> Option<::std::cmp::Ordering> {
                (self.real, self.imag).partial_cmp(&(other.real, other.imag))
            }
        }

        impl From<$t> for Rect {
            #[inline(always)]
            fn from(x: $t) -> Rect { Rect { real: x, imag: 0.0 } }
        }

        impl ::std::ops::Mul<Rect> for Rect {
            type Output = Rect;

            #[inline(always)]
            fn mul(self, other: Rect) -> Rect {
                Rect {
                    real: self.real * other.real - self.imag * other.imag,
                    imag: self.real * other.imag + self.imag * other.real,
                }
            }
        }

        impl ::std::ops::Add<Rect> for Rect {
            type Output = Rect;

            #[inline(always)] // (otherwise it doesn't get inlined into `Ket::overlap`)
            fn add(self, other: Rect) -> Rect {
                Rect {
                    real: self.real + other.real,
                    imag: self.imag + other.imag,
                }
            }
        }

        impl ::std::ops::Sub<Rect> for Rect {
            type Output = Rect;

            #[inline(always)]
            fn sub(self, other: Rect) -> Rect {
                Rect {
                    real: self.real - other.real,
                    imag: self.imag - other.imag,
                }
            }
        }
    };
}

pub(crate) mod lossless {
    // N.B. PartialOrd is very deliberately not implemented because
    //      it could be a huge footgun. Generic code will have to use
    //      function-based APIs and `lexical_cmp`.
    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct Rect {
        pub real: f64,
        pub imag: f64,
    }

    rect_common_impls!(f64);
}

pub(crate) mod compact {
    use ::std::ops::Mul;

    // N.B. PartialOrd is very deliberately not implemented because
    //      it could be a huge footgun. Generic code will have to use
    //      function-based APIs and `lexical_cmp`.
    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct Rect {
        pub real: f32,
        pub imag: f32,
    }

    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct Polar {
        pub abs: f32,
        pub phase: u8,
    }

    impl Polar {
        #[inline(always)]
        pub fn zero() -> Polar { 0_f32.into() }
        #[inline(always)]
        pub fn one() -> Polar { 1_f32.into() }

        #[inline(always)]
        pub fn from_phase_byte(byte: u8) -> Polar {
            Polar { abs: 1.0, phase: byte }
        }

        #[inline(always)]
        pub fn sqnorm(self) -> f32 { self.abs * self.abs }
        #[inline(always)]
        pub fn conj(self) -> Polar {
            Polar {
                abs: self.abs,
                phase: self.phase.wrapping_neg(),
            }
        }
        #[inline]
        pub fn to_rect(self, table: &PhaseTable) -> Rect {
            Rect {
                real: self.abs * table.cos(self.phase),
                imag: self.abs * table.sin(self.phase),
            }
        }
    }

    impl Mul<Polar> for Polar {
        type Output = Polar;

        #[inline(always)]
        fn mul(self, other: Polar) -> Polar {
            Polar {
                abs: self.abs * other.abs,
                phase: self.phase.wrapping_add(other.phase),
            }
        }
    }

    impl From<f32> for Polar {
        #[inline]
        fn from(x: f32) -> Polar { Polar { abs: x, phase: 0 } }
    }

    pub struct PhaseTable {
        // TODO: these should probably be [f64; 256] for bounds check elimination
        tau: Vec<f32>,
        radians: Vec<f32>,
        sin: Vec<f32>,
        cos: Vec<f32>,
    }

    lazy_static! {
        static ref PHASE_TABLE: PhaseTable = PhaseTable::compute();
    }

    impl PhaseTable {
        pub fn compute() -> PhaseTable {
            let tau: Vec<_> = (0u16..256).map(|i| i as f32 / 256.0).collect();
            let radians: Vec<_> = tau.iter().map(|&x| x * 2.0 * ::std::f32::consts::PI).collect();
            let sin: Vec<_> = radians.iter().map(|&x| x.sin()).collect();
            let cos: Vec<_> = radians.iter().map(|&x| x.cos()).collect();
            PhaseTable { tau, radians, sin, cos }
        }

        #[inline]
        pub fn get<'a>() -> &'a PhaseTable { &*PHASE_TABLE }

        #[inline]
        pub fn sin(&self, phase: u8) -> f32 { self.sin[phase as usize] }
        #[inline]
        pub fn cos(&self, phase: u8) -> f32 { self.cos[phase as usize] }
        #[inline]
        pub fn radians(&self, phase: u8) -> f32 { self.radians[phase as usize] }
        #[inline]
        pub fn fraction(&self, phase: u8) -> f32 { self.tau[phase as usize] }

        pub fn nearest_phase(&self, phase: f64) -> u8 {
            let x = (phase as f32 / self.radians(1)).round();
            let x = ((x % 256.) + 256.) % 256.;
            x as u8
        }
    }

    rect_common_impls!(f32);
}
