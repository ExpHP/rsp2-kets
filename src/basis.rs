macro_rules! impl_common_trash {
    (
        types: [$Basis:ident, $Ket:ident, $KetRef:ident]
        traits: [$AsKetRef:ident]
        elements: [$Complex:ident { $a:ident : $A:path, $b:ident : $B:path }]
        real: [$Real:ident]
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
            #[inline]
            pub fn new($a: Vec<$A>, $b: Vec<$B>) -> Self {
                $Ket { $a, $b }
            }
        }

        impl $Ket {
            #[inline]
            pub fn as_ref(&self) -> $KetRef {
                let $Ket { ref $a, ref $b } = *self;
                $KetRef { $a, $b }
            }

            // methods forwarded to KetRef

            #[inline]
            pub fn len(&self) -> usize { self.as_ref().len() }

            #[inline]
            pub fn $a(&self) -> &[$A] { &self.$a }
            #[inline]
            pub fn $b(&self) -> &[$B] { &self.$b }

            #[inline]
            pub fn sqnorm(&self) -> $Real { self.as_ref().sqnorm() }
            #[inline]
            pub fn norm(&self) -> $Real { self.as_ref().norm() }

            #[inline]
            pub fn to_normalized(&self) -> $Ket { self.as_ref().to_normalized() }

            // can't do Index because we can't return a borrow
            #[inline]
            pub fn at(&self, i: usize) -> $Complex { self.as_ref().at(i) }
            #[inline]
            pub fn overlap<K: $AsKetRef>(self, other: &K) -> $Real { self.as_ref().overlap(other) }
            #[inline]
            pub fn iter(&self) -> Iter { self.as_ref().iter() }
        }

        impl $Ket {
            pub fn into_normalized(self) -> $Ket {
                let norm = self.norm();
                self.div_real(norm)
            }
        }

        impl IntoIterator for $Ket {
            type Item = $Complex;
            type IntoIter = IntoIter;

            #[inline]
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

        pub trait $AsKetRef {
            fn as_ket_ref(&self) -> $KetRef;
        }

        impl<'a> $AsKetRef for $KetRef<'a> {
            #[inline]
            fn as_ket_ref(&self) -> $KetRef { *self }
        }
        impl $AsKetRef for $Ket {
            #[inline]
            fn as_ket_ref(&self) -> $KetRef { self.as_ref() }
        }

        /// A not-owned ket.
        #[derive(Debug, Copy, Clone)]
        pub struct $KetRef<'a> {
            pub(crate) $a: &'a [$A],
            pub(crate) $b: &'a [$B],
        }

        impl<'a> $KetRef<'a> {
            #[inline]
            pub fn new($a: &'a [$A], $b: &'a [$B]) -> Self {
                $KetRef { $a, $b }
            }

            // can't do Index because we can't return a borrow
            #[inline]
            pub fn at(&self, i: usize) -> $Complex {
                $Complex {
                    $a: self.$a[i],
                    $b: self.$b[i],
                }
            }

            #[inline]
            pub fn len(&self) -> usize { self.$a.len() }
            #[inline]
            pub fn $a(&self) -> &[$A] { self.$a }
            #[inline]
            pub fn $b(&self) -> &[$B] { self.$b }

            #[inline]
            pub fn sqnorm(&self) -> $Real { self.overlap(self) }
            #[inline]
            pub fn norm(&self) -> $Real { self.sqnorm().sqrt() }

            #[inline]
            pub fn to_owned(&self) -> $Ket {
                $Ket {
                    $a: self.$a.to_owned(),
                    $b: self.$b.to_owned(),
                }
            }

            #[inline]
            pub fn to_normalized(&self) -> $Ket { self.to_owned().into_normalized() }

            #[inline]
            pub fn iter(&self) -> Iter<'a> {
                let &$KetRef { $a, $b } = self;
                Box::new($a.iter().zip($b).map(|(&$a, &$b)| $Complex { $a, $b }))
            }
        }

        impl<'a> IntoIterator for $KetRef<'a> {
            type Item = $Complex;
            type IntoIter = Iter<'a>;
            #[inline]
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
            #[inline]
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
            #[inline]
            pub fn width(&self) -> usize { self.width }
            /// Number of kets
            #[inline]
            pub fn rank(&self) -> usize { self.data.len() / (2 * self.width) }
            #[inline]
            pub fn ket(&self, i: usize) -> KetRef {
                let w = self.width;
                KetRef {
                    real: &self.data[w * (2 * i + 0) .. w * (2 * i + 1)],
                    imag: &self.data[w * (2 * i + 1) .. w * (2 * i + 2)],
                }
            }

            #[inline]
            pub fn iter(&self) -> Iter {
                Box::new((0..self.rank()).map(move |i| self.ket(i)))
            }

            pub fn lossy_compress(&self) -> ::compact::Basis {
                use ::complex::compact::PhaseTable;

                let table = PhaseTable::get();
                let mut abs = vec![];
                let mut phase = vec![];
                for ket in self.iter() {
                    for c in ket.iter() {
                        phase.push(table.nearest_phase(c.imag.atan2(c.real)));
                        abs.push(c.abs() as f32);
                    }
                }
                let width = self.width;
                ::basis::compact::basis::Cereal { width, phase, abs }.validate()
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
            #[inline]
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
            real: [f64]
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

        impl Ket {
            fn div_real(mut self, factor: f64) -> Ket {
                for x in &mut self.real { *x /= factor; }
                for x in &mut self.imag { *x /= factor; }
                self
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
        //  - abs.len() == phase.len()
        //  - abs.len() is divisible by width
        // contracts that aren't strictly protected as invariants:
        //  - eigenvectors SHOULD be orthogonal
        //  - eigenvectors SHOULD be normalized
        #[derive(Debug, Clone)]
        #[derive(PartialEq)]
        pub struct Basis {
            width: usize,
            abs:  Vec<f32>,
            phase: Vec<u8>,
        }

        impl Basis {
            #[inline]
            pub fn new(abs: Vec<f32>, phase: Vec<u8>, width: usize) -> Basis {
                Cereal { abs, phase, width }.validate()
            }

            // takes a tuple to be forward-compatible with
            // a possible overload in the future for KetRef
            pub fn insert(&mut self, (norm, phase): (&[f32], &[u8])) {
                assert_eq!(self.width, norm.len());
                assert_eq!(self.width, phase.len());
                self.abs.extend_from_slice(norm);
                self.phase.extend_from_slice(phase);
            }

            #[inline]
            pub fn rank(&self) -> usize { self.abs.len() / self.width }
            #[inline]
            pub fn width(&self) -> usize { self.width }
            #[inline]
            pub fn ket(&self, i: usize) -> KetRef {
                let w = self.width;
                KetRef {
                    abs:   &self.abs  [w * i .. w * (i + 1)],
                    phase: &self.phase[w * i .. w * (i + 1)],
                }
            }

            #[inline]
            pub fn iter(&self) -> Iter {
                Box::new((0..self.rank()).map(move |i| self.ket(i)))
            }
        }

        /// Raw data type with no invariants, for serialization
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
        #[derive(Debug, Clone, PartialEq)]
        pub(crate) struct Cereal {
            pub width: usize,
            pub abs:   Vec<f32>,
            pub phase: Vec<u8>,
        }

        impl Cereal {
            #[inline]
            pub fn validate(self) -> Basis {
                let Cereal { width, abs, phase } = self;
                assert_eq!(abs.len(), phase.len());
                assert_eq!(abs.len() % width, 0);
                Basis { width, abs, phase }
            }
        }

        impl Basis {
            #[cfg(feature = "serde")]
            pub(crate) fn cereal(self) -> Cereal {
                let Basis { width, abs, phase } = self;
                Cereal { width, abs, phase }
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
            elements: [Polar { abs: f32, phase: u8 }]
            real: [f32]
        }

        impl<'a> KetRef<'a> {
            pub fn overlap<K: AsKetRef>(self, other: &K) -> f32 {
                let other = other.as_ket_ref();
                assert_eq!(self.abs.len(), other.abs.len());
                let table = PhaseTable::get();
                (0..self.abs.len())
                    .map(|i| (self.at(i).conj() * other.at(i)).to_rect(table))
                    .fold(Rect::zero(), |a, b| a + b)
                    .sqnorm()
            }
        }

        impl Ket {
            fn div_real(mut self, factor: f32) -> Ket {
                for x in &mut self.abs { *x /= factor; }
                self
            }
        }
    }
}
