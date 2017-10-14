


/// Full double-precision rectangular representation,
/// for applications where precision matters.
pub mod lossless {
    pub use self::basis::Basis;
    pub use self::basis::Raw as RawBasis;
    pub mod basis {
        use super::KetRef;

        pub type Iter<'a> = Box<Iterator<Item=KetRef<'a>> + 'a>;

        // invariants:
        //  - data.len() is divisible by 2 * width
        // contracts that aren't strictly protected as invariants:
        //  - eigenvectors SHOULD be orthogonal
        //  - eigenvectors SHOULD be normalized
        #[derive(Debug, Clone, PartialEq)]
        pub struct Basis {
            width: usize,
            // reals of ket 1, then imags of ket 1, then reals of ket 2...
            data: Vec<f64>,
        }

        impl Basis {
            pub fn new(data: Vec<f64>, width: usize) -> Basis {
                Raw { data, width }.validate()
            }

            // takes a tuple to be forward-compatible with
            // a possible overload in the future for KetRef
            pub fn insert(&mut self, (real, imag): (&[f64], &[f64])) {
                assert_eq!(self.width, real.len());
                assert_eq!(self.width, imag.len());
                self.data.extend_from_slice(real);
                self.data.extend_from_slice(imag);
            }

            /// Number of dimensions in a ket.
            pub fn width(&self) -> usize { self.width }
            /// Number of kets
            pub fn rank(&self) -> usize { self.data.len() / (2 * self.width) }
            pub fn ket(&self, i: usize) -> KetRef {
                let w = self.width;
                KetRef {
                    real: &self.data[w * (2 * i + 0) .. w * (2 * i + 1)],
                    imag: &self.data[w * (2 * i + 1) .. w * (2 * i + 2)],
                }
            }

            pub fn iter(&self) -> Iter {
                Box::new((0..self.rank()).map(move |i| self.ket(i)))
            }

            pub fn lossy_compress(&self) -> ::basis::compact::Basis {
                use ::complex::compact::PhaseTable;

                let table = PhaseTable::get();
                let mut norm = vec![];
                let mut phase = vec![];
                for ket in self.iter() {
                    for c in ket.iter() {
                        phase.push(table.nearest_phase(c.imag.atan2(c.real)));
                        norm.push(c.sqnorm().sqrt() as f32);
                    }
                }
                let width = self.width;
                ::basis::compact::basis::Raw { width, phase, norm }.validate()
            }
        }

        /// Raw data type with no invariants, for serialization
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
        #[derive(Debug, Clone, PartialEq)]
        pub struct Raw {
            pub width: usize,
            pub data: Vec<f64>,
        }

        impl Raw {
            pub fn validate(self) -> Basis {
                let Raw { width, data } = self;
                assert_eq!(data.len() % (2 * width), 0);
                Basis { width, data }
            }
        }

        impl Basis {
            pub fn raw(self) -> Raw {
                let Basis { width, data } = self;
                Raw { width, data }
            }
        }
    }

    pub use self::ket::Ket;
    pub use self::ket::KetRef;
    pub use self::ket::AsKetRef;
    pub mod ket {
        use ::complex::lossless::Rect;
        pub type IntoIter = Box<Iterator<Item=Rect>>;
        pub type Iter<'a> = Box<Iterator<Item=Rect> + 'a>;

        /// An owned ket.
        #[derive(Debug, Clone)]
        pub struct Ket {
            pub(crate) real: Vec<f64>,
            pub(crate) imag: Vec<f64>,
        }

        impl Ket {
            pub fn as_ref(&self) -> KetRef {
                let Ket { ref real, ref imag } = *self;
                KetRef { real, imag }
            }

            pub fn at(&self, i: usize) -> Rect { self.as_ref().at(i) }
            pub fn overlap<K: AsKetRef>(self, other: &K) -> f64 { self.as_ref().overlap(other) }
            pub fn iter(&self) -> Iter { self.as_ref().iter() }
        }

        impl IntoIterator for Ket {
            type Item = Rect;
            type IntoIter = IntoIter;
            fn into_iter(self) -> Self::IntoIter {
                let Ket { real, imag } = self;
                Box::new(real.into_iter().zip(imag).map(|(real,imag)| Rect { real, imag }))
            }
        }

        pub trait AsKetRef { fn as_ket_ref(&self) -> KetRef; }
        impl<'a> AsKetRef for KetRef<'a> { fn as_ket_ref(&self) -> KetRef { *self } }
        impl AsKetRef for Ket { fn as_ket_ref(&self) -> KetRef { self.as_ref() } }

        #[derive(Debug, Copy, Clone)]
        pub struct KetRef<'a> {
            pub(crate) real: &'a [f64],
            pub(crate) imag: &'a [f64],
        }

        impl<'a> KetRef<'a> {
            fn at(&self, i: usize) -> Rect {
                Rect {
                    real: self.real[i],
                    imag: self.imag[i],
                }
            }

            pub fn overlap<K: AsKetRef>(self, other: &K) -> f64 {
                let other = other.as_ket_ref();
                assert_eq!(self.real.len(), other.real.len());
                (0..self.real.len())
                    .map(|i| (self.at(i).conj() * other.at(i)))
                    .fold(Rect::zero(), |a,b| a + b)
                    .sqnorm()
            }

            pub fn iter(&self) -> Iter<'a> {
                let &KetRef { real, imag } = self;
                Box::new(real.iter().zip(imag).map(|(&real, &imag)| Rect { real, imag }))
            }
        }
    }
}

/// A lossily-compressed form that is amenable to further compression.
///
/// Suitable for e.g. band uncrossing.
pub mod compact {
    pub use self::basis::Basis;
    pub use self::basis::Raw as RawBasis;
    pub mod basis {
        use super::KetRef;

        pub type Iter<'a> = Box<Iterator<Item=KetRef<'a>> + 'a>;

        // invariants:
        //  - norm.len() == phase.len()
        //  - norm.len() is divisible by width
        // contracts that aren't strictly protected as invariants:
        //  - eigenvectors SHOULD be orthogonal
        //  - eigenvectors SHOULD be normalized
        #[derive(Debug, Clone)]
        #[derive(PartialEq)]
        pub struct Basis {
            width: usize,
            norm:  Vec<f32>,
            phase: Vec<u8>,
        }

        impl Basis {
            pub fn new(norm: Vec<f32>, phase: Vec<u8>, width: usize) -> Basis {
                Raw { norm, phase, width }.validate()
            }

            // takes a tuple to be forward-compatible with
            // a possible overload in the future for KetRef
            pub fn insert(&mut self, (norm, phase): (&[f32], &[u8])) {
                assert_eq!(self.width, norm.len());
                assert_eq!(self.width, phase.len());
                self.norm.extend_from_slice(norm);
                self.phase.extend_from_slice(phase);
            }

            pub fn rank(&self) -> usize { self.norm.len() / self.width }
            pub fn width(&self) -> usize { self.width }
            pub fn ket(&self, i: usize) -> KetRef {
                let w = self.width;
                KetRef {
                    norm:  &self.norm [w * i .. w * (i + 1)],
                    phase: &self.phase[w * i .. w * (i + 1)],
                }
            }

            pub fn iter(&self) -> Iter {
                Box::new((0..self.rank()).map(move |i| self.ket(i)))
            }
        }

        /// Raw data type with no invariants, for serialization
        #[derive(Debug, Clone, PartialEq)]
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
        pub struct Raw {
            pub width: usize,
            pub norm:  Vec<f32>,
            pub phase: Vec<u8>,
        }

        impl Raw {
            pub fn validate(self) -> Basis {
                let Raw { width, norm, phase } = self;
                assert_eq!(norm.len(), phase.len());
                assert_eq!(norm.len() % width, 0);
                Basis { width, norm, phase }
            }
        }

        impl Basis {
            pub fn raw(self) -> Raw {
                let Basis { width, norm, phase } = self;
                Raw { width, norm, phase }
            }
        }
    }

    pub use self::ket::Ket;
    pub use self::ket::KetRef;
    pub use self::ket::AsKetRef;
    pub mod ket {
        use ::complex::compact::{Rect, Polar, PhaseTable};
        pub type IntoIter = Box<Iterator<Item=Polar>>;
        pub type Iter<'a> = Box<Iterator<Item=Polar> + 'a>;

        /// An owned ket.
        #[derive(Debug, Clone)]
        pub struct Ket {
            pub(crate) norm:  Vec<f32>,
            pub(crate) phase: Vec<u8>,
        }

        impl Ket {
            pub fn as_ref(&self) -> KetRef {
                let Ket { ref norm, ref phase } = *self;
                KetRef { norm, phase }
            }

            pub fn at(&self, i: usize) -> Polar { self.as_ref().at(i) }
            pub fn overlap<K: AsKetRef>(self, other: &K) -> f64 { self.as_ref().overlap(other) }
        }

        impl IntoIterator for Ket {
            type Item = Polar;
            type IntoIter = IntoIter;
            fn into_iter(self) -> Self::IntoIter {
                let Ket { norm, phase } = self;
                Box::new(norm.into_iter().zip(phase)
                    .map(|(norm, phase)| Polar { norm, phase }))
            }
        }

        pub trait AsKetRef { fn as_ket_ref(&self) -> KetRef; }
        impl<'a> AsKetRef for KetRef<'a> { fn as_ket_ref(&self) -> KetRef { *self } }
        impl AsKetRef for Ket { fn as_ket_ref(&self) -> KetRef { self.as_ref() } }

        #[derive(Debug, Copy, Clone)]
        pub struct KetRef<'a> {
            pub(crate) norm:  &'a [f32],
            pub(crate) phase: &'a [u8],
        }

        impl<'a> KetRef<'a> {
            fn at(&self, i: usize) -> Polar {
                Polar {
                    norm: self.norm[i],
                    phase: self.phase[i],
                }
            }

            pub fn overlap<K: AsKetRef>(self, other: &K) -> f64 {
                let other = other.as_ket_ref();
                assert_eq!(self.norm.len(), other.norm.len());
                let table = PhaseTable::get();
                (0..self.norm.len())
                    .map(|i| (self.at(i).conj() * other.at(i)).to_rect(table))
                    .fold(Rect::zero(), |a,b| a + b)
                    .sqnorm() as f64
            }

            pub fn iter(&self) -> Iter<'a> {
                let KetRef { norm, phase } = *self;
                Box::new(norm.into_iter().zip(phase)
                    .map(|(&norm, &phase)| Polar { norm, phase }))
            }
        }
    }
}
