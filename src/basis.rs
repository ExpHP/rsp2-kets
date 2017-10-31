macro_rules! impl_common_trash {
    (
        types: [$Basis:ident, $Ket:ident, $KetRef:ident]
        traits: [$AsKetRef:ident]
        elements: [$Complex:ident { $a:ident : $A:path, $b:ident : $B:path }]
    ) => {
        pub type IntoIter = Box<Iterator<Item=$Complex>>;
        pub type Iter<'a> = Box<Iterator<Item=$Complex> + 'a>;

        /// An owned ket.
        #[derive(Debug, Clone)]
        pub struct $Ket {
            pub(crate) $a: Vec<$A>,
            pub(crate) $b: Vec<$B>,
        }

        impl $Ket {
            pub fn new($a: Vec<$A>, $b: Vec<$B>) -> Self {
                $Ket { $a, $b }
            }

            pub fn as_ref(&self) -> $KetRef {
                let $Ket { ref $a, ref $b } = *self;
                $KetRef { $a, $b }
            }

            pub fn len(&self) -> usize { self.as_ref().len() }

            pub fn $a(&self) -> &[$A] { &self.$a }
            pub fn $b(&self) -> &[$B] { &self.$b }

            // can't do Index because we can't return a borrow
            pub fn at(&self, i: usize) -> $Complex { self.as_ref().at(i) }
            pub fn overlap<K: $AsKetRef>(self, other: &K) -> f64 { self.as_ref().overlap(other) }
            pub fn iter(&self) -> Iter { self.as_ref().iter() }
        }

        impl IntoIterator for $Ket {
            type Item = $Complex;
            type IntoIter = IntoIter;
            fn into_iter(self) -> Self::IntoIter {
                let $Ket { $a, $b } = self;
                Box::new($a.into_iter().zip($b).map(|($a, $b)| $Complex { $a, $b }))
            }
        }

        impl ::std::iter::FromIterator<$Complex> for $Ket {
            #[inline]
            fn from_iter<I>(iter: I) -> Self
            where I: IntoIterator<Item=$Complex>
            { iter.into_iter().map(|$Complex { $a, $b }| ($a, $b)).collect() }
        }

        impl ::std::iter::FromIterator<($A, $B)> for $Ket {
            #[inline]
            fn from_iter<I>(iter: I) -> Self
            where I: IntoIterator<Item=($A, $B)>
            {
                let ($a, $b) = iter.into_iter().unzip();
                $Ket { $a, $b }
            }
        }

        pub trait $AsKetRef { fn as_ket_ref(&self) -> $KetRef; }
        impl<'a> $AsKetRef for $KetRef<'a> { fn as_ket_ref(&self) -> $KetRef { *self } }
        impl $AsKetRef for $Ket { fn as_ket_ref(&self) -> $KetRef { self.as_ref() } }

        /// A not-owned ket.
        #[derive(Debug, Copy, Clone)]
        pub struct $KetRef<'a> {
            pub(crate) $a: &'a [$A],
            pub(crate) $b: &'a [$B],
        }

        impl<'a> $KetRef<'a> {
            pub fn new($a: &'a [$A], $b: &'a [$B]) -> Self {
                $KetRef { $a, $b }
            }

            // can't do Index because we can't return a borrow
            pub fn at(&self, i: usize) -> $Complex {
                $Complex {
                    $a: self.$a[i],
                    $b: self.$b[i],
                }
            }

            pub fn len(&self) -> usize { self.$a.len() }
            pub fn $a(&self) -> &[$A] { self.$a }
            pub fn $b(&self) -> &[$B] { self.$b }

            pub fn iter(&self) -> Iter<'a> {
                let &$KetRef { $a, $b } = self;
                Box::new($a.iter().zip($b).map(|(&$a, &$b)| $Complex { $a, $b }))
            }
        }

        impl<'a> IntoIterator for $KetRef<'a> {
            type Item = $Complex;
            type IntoIter = Iter<'a>;
            fn into_iter(self) -> Self::IntoIter { self.iter() }
        }
    };
}

macro_rules! forward_serde_impls {
    (
        serialize: [$Type:ident :: $cereal:ident]
        deserialize: [$Cereal:ident :: $validate:ident]
    ) => {
        #[cfg(feature = "serde")]
        impl ::serde::Serialize for $Type
        {
            fn serialize<S>(&self, s: S) -> Result<S::Ok, S::Error>
            where S: ::serde::Serializer,
            { self.clone().$cereal().serialize(s) }
        }

        #[cfg(feature = "serde")]
        impl<'de> ::serde::Deserialize<'de> for $Type
        {
            fn deserialize<D>(d: D) -> Result<$Type, D::Error>
            where D: ::serde::Deserializer<'de>
            {
                let cereal = <$Cereal as ::serde::Deserialize>::deserialize(d)?;
                Ok(cereal.$validate())
            }
        }
    };
}

/// Full double-precision rectangular representation,
/// for applications where precision matters.
pub(crate) mod lossless {
    pub(crate) mod basis {
        use super::ket::KetRef;

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
                Cereal { data, width }.validate()
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

            pub fn lossy_compress(&self) -> ::compact::Basis {
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
                ::basis::compact::basis::Cereal { width, phase, norm }.validate()
            }
        }

        /// Raw data type with no invariants, for serialization
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
        #[derive(Debug, Clone, PartialEq)]
        pub(crate) struct Cereal {
            pub width: usize,
            pub data: Vec<f64>,
        }

        impl Cereal {
            pub fn validate(self) -> Basis {
                let Cereal { width, data } = self;
                assert_eq!(data.len() % (2 * width), 0);
                Basis { width, data }
            }
        }

        impl Basis {
            /// Prepare for serialization.
            #[cfg(feature = "serde")]
            pub(crate) fn cereal(self) -> Cereal {
                let Basis { width, data } = self;
                Cereal { width, data }
            }
        }

        forward_serde_impls!{
            serialize: [Basis::cereal]
            deserialize: [Cereal::validate]
        }
    }

    pub(crate) mod ket {
        use ::complex::lossless::Rect;

        impl_common_trash! {
            types: [Basis, Ket, KetRef]
            traits: [AsKetRef]
            elements: [Rect { real: f64, imag: f64 }]
        }

        impl<'a> KetRef<'a> {
            pub fn overlap<K: AsKetRef>(self, other: &K) -> f64 {
                let other = other.as_ket_ref();
                assert_eq!(self.real.len(), other.real.len());
                (0..self.real.len())
                    .map(|i| (self.at(i).conj() * other.at(i)))
                    .fold(Rect::zero(), |a,b| a + b)
                    .sqnorm()
            }
        }
    }
}

/// A lossily-compressed form that is amenable to further compression.
///
/// Suitable for e.g. band uncrossing.
pub(crate) mod compact {
    pub(crate) mod basis {
        use super::ket::KetRef;

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
                Cereal { norm, phase, width }.validate()
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
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
        #[derive(Debug, Clone, PartialEq)]
        pub(crate) struct Cereal {
            pub width: usize,
            pub norm:  Vec<f32>,
            pub phase: Vec<u8>,
        }

        impl Cereal {
            pub fn validate(self) -> Basis {
                let Cereal { width, norm, phase } = self;
                assert_eq!(norm.len(), phase.len());
                assert_eq!(norm.len() % width, 0);
                Basis { width, norm, phase }
            }
        }

        impl Basis {
            #[cfg(feature = "serde")]
            pub(crate) fn cereal(self) -> Cereal {
                let Basis { width, norm, phase } = self;
                Cereal { width, norm, phase }
            }
        }


        forward_serde_impls!{
            serialize: [Basis::cereal]
            deserialize: [Cereal::validate]
        }
    }

    pub(crate) mod ket {
        use ::complex::compact::{Rect, Polar, PhaseTable};

        impl_common_trash! {
            types: [Basis, Ket, KetRef]
            traits: [AsKetRef]
            elements: [Polar { norm: f32, phase: u8 }]
        }

        impl<'a> KetRef<'a> {
            pub fn overlap<K: AsKetRef>(self, other: &K) -> f64 {
                let other = other.as_ket_ref();
                assert_eq!(self.norm.len(), other.norm.len());
                let table = PhaseTable::get();
                (0..self.norm.len())
                    .map(|i| (self.at(i).conj() * other.at(i)).to_rect(table))
                    .fold(Rect::zero(), |a,b| a + b)
                    .sqnorm() as f64
            }
        }
    }
}
