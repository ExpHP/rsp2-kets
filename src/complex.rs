
macro_rules! rect_common_impls {
    ($t:ty) => {
        impl Rect {
            pub fn zero() -> Rect { 0.0.into() }
            pub fn one() -> Rect { 1.0.into() }
            pub fn i() -> Rect { Rect { real: 0.0, imag: 1.0 } }

            pub fn from_phase(radians: $t) -> Rect {
                Rect {
                    real: radians.cos(),
                    imag: radians.sin(),
                }
            }

            pub fn sqnorm(self) -> $t { self.real * self.real + self.imag * self.imag }
            pub fn conj(self) -> Rect {
                Rect {
                    real:  self.real,
                    imag: -self.imag,
                }
            }
        }

        impl From<$t> for Rect {
            fn from(x: $t) -> Rect { Rect { real: x, imag: 0.0 } }
        }

        impl ::std::ops::Mul<Rect> for Rect {
            type Output = Rect;
            fn mul(self, other: Rect) -> Rect {
                Rect {
                    real: self.real * other.real - self.imag * other.imag,
                    imag: self.real * other.imag + self.imag * other.real,
                }
            }
        }

        impl ::std::ops::Add<Rect> for Rect {
            type Output = Rect;
            fn add(self, other: Rect) -> Rect {
                Rect {
                    real: self.real + other.real,
                    imag: self.imag + other.imag,
                }
            }
        }
    };
}

pub(crate) mod lossless {
    pub struct Rect {
        pub real: f64,
        pub imag: f64,
    }

    rect_common_impls!(f64);
}

pub(crate) mod compact {
    use ::std::ops::Mul;

    pub struct Rect {
        pub real: f32,
        pub imag: f32,
    }

    pub struct Polar {
        pub norm: f32,
        pub phase: u8,
    }

    impl Polar {
        pub fn zero() -> Polar { 0_f32.into() }
        pub fn one() -> Polar { 1_f32.into() }

        pub fn from_phase_byte(byte: u8) -> Polar {
            Polar { norm: 1.0, phase: byte }
        }

        pub fn sqnorm(self) -> f32 { self.norm * self.norm }
        pub fn conj(self) -> Polar {
            Polar {
                norm: self.norm,
                phase: self.phase.wrapping_neg(),
            }
        }
        pub fn to_rect(self, table: &PhaseTable) -> Rect {
            Rect {
                real: self.norm * table.cos(self.phase),
                imag: self.norm * table.sin(self.phase),
            }
        }
    }

    impl Mul<Polar> for Polar {
        type Output = Polar;
        fn mul(self, other: Polar) -> Polar {
            Polar {
                norm: self.norm * other.norm,
                phase: self.phase.wrapping_add(other.phase),
            }
        }
    }

    impl From<f32> for Polar {
        fn from(x: f32) -> Polar { Polar { norm: x, phase: 0 } }
    }

    pub struct PhaseTable {
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

        pub fn get<'a>() -> &'a PhaseTable { &*PHASE_TABLE }

        pub fn sin(&self, phase: u8) -> f32 { self.sin[phase as usize] }
        pub fn cos(&self, phase: u8) -> f32 { self.cos[phase as usize] }
        pub fn radians(&self, phase: u8) -> f32 { self.radians[phase as usize] }
        pub fn fraction(&self, phase: u8) -> f32 { self.tau[phase as usize] }

        pub fn nearest_phase(&self, phase: f64) -> u8 {
            let x = (phase as f32 / self.radians(1)).round();
            let x = ((x % 256.) + 256.) % 256.;
            x as u8
        }
    }

    rect_common_impls!(f32);
}
