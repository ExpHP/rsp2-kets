macro_rules! impl_common_trash {
    (
        types: [$Basis:ident, $Ket:ident, $KetRef:ident]
        traits: [$AsKetRef:ident]
        // what the datatype conceptually contains
        elements: [$Complex:ident { $a:ident : $A:path, $b:ident : $B:path }]
        // even polar datatypes use Rect where summation is required
        rect: [$Rect:ident]
        // type of probability
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

            /// Computes `<self|self>`
            #[inline]
            pub fn sqnorm(&self) -> $Real { self.as_ref().sqnorm() }
            #[inline]
            pub fn norm(&self) -> $Real { self.as_ref().norm() }

            /// Computes `<self|other>` (i.e. `self` becomes the bra)
            #[inline]
            pub fn dot<K: $AsKetRef>(&self, other: K) -> $Rect { self.as_ref().dot(other) }

            #[inline]
            pub fn to_normalized(&self) -> $Ket { self.as_ref().to_normalized() }

            // can't do Index because we can't return a borrow
            #[inline]
            pub fn at(&self, i: usize) -> $Complex { self.as_ref().at(i) }
            /// Computes `<self|other><other|self>`
            #[inline]
            pub fn overlap<K: $AsKetRef>(&self, other: K) -> $Real { self.as_ref().overlap(other) }
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

        impl<'a> IntoIterator for &'a super::basis::$Basis {
            type Item = $KetRef<'a>;
            type IntoIter = super::basis::Iter<'a>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter { self.iter() }
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
        impl<'a, K: $AsKetRef> $AsKetRef for &'a K {
            #[inline]
            fn as_ket_ref(&self) -> $KetRef { (**self).as_ket_ref() }
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
            pub fn norm(&self) -> $Real { self.sqnorm().sqrt() }
            /// Computes `<self|other><other|self>`
            #[inline]
            pub fn overlap<K: $AsKetRef>(&self, other: K) -> $Real { self.dot(other).sqnorm() }

            #[inline]
            pub fn to_owned(&self) -> $Ket {
                $Ket {
                    $a: self.$a.to_owned(),
                    $b: self.$b.to_owned(),
                }
            }

            #[inline]
            pub fn scale(&self, c: $Complex) -> $Ket {
                self.iter().map(|x| x * c).collect()
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

            // Reminder to self:
            // Many sources misrepresent the Modified Gram Schmidt method by implying
            // that its key difference from Classical Gram Schmidt is in its order of
            // iteration, describing MGS as a method that "looks ahead" and uses each vector
            // to orthogonalize vectors after it.
            //
            // If you work out the data dependencies, however, you will find that this difference
            // is irrelevant, and that the change in order merely makes the (already present) opportunities
            // for parallelism more obvious.
            //
            // (I won't bother with parallelizing this yet, since we can do *even better* than a
            //  parallel inner loop by explicitly using rayon::join)
            //
            /// Orthonormalize a basis using the Modified Gram Schmidt method.
            pub fn orthonormalize(&self) -> Basis {
                let mut out = Basis::new(vec![], self.width);
                for ket in self {
                    let mut ket = ket.to_owned();
                    for bra in &out {
                        let projected = ket.projected_onto(bra);

                        // ket -= projected
                        ket.real.iter_mut().zip(projected.real).for_each(|(a, b)| { *a -= b; });
                        ket.imag.iter_mut().zip(projected.imag).for_each(|(a, b)| { *a -= b; });
                    }
                    ket = ket.into_normalized();
                    out.insert((ket.real(), ket.imag()));
                }
                out
            }
        }

        #[test]
        fn test_orthonormalize() {
            let dim = 200;
            let num_kets = 30;

            let data = (0..dim * num_kets * 2).map(|_| 0.5 - ::rand::random::<f64>()).collect();
            let basis = Basis::new(data, dim).orthonormalize();
            for (i, ket) in basis.iter().enumerate() {
                for (j, bra) in basis.iter().enumerate() {
                    if i == j {
                        assert!(f64::abs(ket.overlap(bra) - 1.0) < 1e-12, "({},{}): {}", i, j, ket.overlap(bra));
                    } else {
                        assert!(f64::abs(ket.overlap(bra)) < 1e-12, "({},{}): {}", i, j, ket.overlap(bra));
                    }
                }
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
            rect: [Rect]
            real: [f64]
        }

        impl<'a> KetRef<'a> {
            /// Computes `<self|other>` (i.e. `self` becomes the bra)
            #[cfg(not(feature = "faster"))]
            pub fn dot<K: AsKetRef>(self, other: K) -> Rect {
                let other = other.as_ket_ref();
                assert_eq!(self.real.len(), other.real.len());
                // NOTE: I'm not sure why, but iterating through all four lists together
                //       seems to be quite a bit faster (factor of 2) compared to taking
                //       two lists at a time and summing their elementwise products, especially
                //       on larger inputs.  I tried measuring cache events with `perf stat` but
                //       saw nothing extraordinary/didn't know what to look for...
                (0..self.real.len())
                    .map(|i| (self.at(i).conj() * other.at(i)))
                    .fold(Rect::zero(), |a,b| a + b)
            }

            #[cfg(feature = "faster")]
            pub fn dot<K: AsKetRef>(self, other: K) -> Rect {
                use ::faster::prelude::*;

                let other = other.as_ket_ref();
                assert_eq!(self.real.len(), other.real.len());
                assert_eq!(self.real.len(), other.imag.len());
                assert_eq!(self.real.len(), self.imag.len());

                let zero = (f64s(0.0), f64s(0.0));
                let add = |(ar, ai): (f64s, f64s), (br, bi): (f64s, f64s)| {
                    (ar + br, ai + bi)
                };
                let simd_map_func = |(ar, ai, br, bi): (f64s, f64s, f64s, f64s)| {
                    let real = ar * br + ai * bi;
                    let imag = ar * bi - ai * br;
                    (real, imag)
                };

                let mut iter = (
                    // (the defaults of zero here may show up in the unaligned remainder,
                    //  where they will harmlessly "contribute" to the sum)
                    self.real.simd_iter(f64s(0.0)),
                    self.imag.simd_iter(f64s(0.0)),
                    other.real.simd_iter(f64s(0.0)),
                    other.imag.simd_iter(f64s(0.0)),
                ).zip();

                let partial_sums = {
                    // deal with aligned elements
                    iter.by_ref()
                        .map(simd_map_func)
                        .fold(zero, add)
                };
                // deal with unaligned remainder
                let leftovers = iter.end().map(|(vs, _)| vs).map(simd_map_func).unwrap_or(zero);
                let (real, imag) = add(partial_sums, leftovers);
                let (real, imag) = (real.sum(), imag.sum());

                Rect { real, imag }
            }

            /// Computes `<self|self>`
            pub fn sqnorm(self) -> f64 {
                self.real.iter().chain(self.imag).map(|x| x*x).sum::<f64>()
            }

            /// Computes `|other><other|self>`
            pub fn projected_onto<K: AsKetRef>(self, other: K) -> Ket {
                let other = other.as_ket_ref();
                other.scale(other.dot(self))
            }
        }

        impl Ket {
            fn div_real(mut self, factor: f64) -> Ket {
                for x in &mut self.real { *x /= factor; }
                for x in &mut self.imag { *x /= factor; }
                self
            }

            #[inline]
            pub fn projected_onto<K: AsKetRef>(&self, other: K) -> Ket {
                self.as_ref().projected_onto(other)
            }
        }

        #[cfg(test)]
        fn random_vector(n: usize) -> Vec<f64> {
            (0..n).map(|_| 0.5 - ::rand::random::<f64>()).collect()
        }

        #[test]
        fn test_normalize() {
            let dim = 200;
            let ket = Ket::new(random_vector(dim), random_vector(dim));
            assert!(1.0 - f64::abs(ket.to_normalized().norm()) < 1e-12);
            assert!(1.0 - f64::abs(ket.as_ref().to_normalized().norm()) < 1e-12);
            assert!(1.0 - f64::abs(ket.into_normalized().norm()) < 1e-12);
        }

        #[test]
        fn test_dot() {
            let a = KetRef::new(&[1.0, 1.0, 3.0], &[0.0, 0.0, 0.0]);
            let b = KetRef::new(&[1.0, 1.0, 2.0], &[0.0, 0.0, 0.0]);
            assert_eq!(a.dot(a), Rect { real: 11.0, imag: 0.0 });
            assert_eq!(a.sqnorm(), 11.0);
            assert_eq!(a.dot(b), Rect { real: 8.0, imag: 0.0 });

            let a = KetRef::new(&[0.0, 0.0, 0.0], &[1.0, 1.0, 3.0]);
            let b = KetRef::new(&[0.0, 0.0, 0.0], &[1.0, 1.0, 2.0]);
            assert_eq!(a.dot(a), Rect { real: 11.0, imag: 0.0 });
            assert_eq!(a.sqnorm(), 11.0);
            assert_eq!(a.dot(b), Rect { real: 8.0, imag: 0.0 });

            let a = KetRef::new(&[0.0, 1.0, 0.0], &[0.0, 0.0, 0.0]);
            let b = KetRef::new(&[0.0, 0.0, 0.0], &[0.0, 1.0, 0.0]);
            assert_eq!(a.dot(b), Rect { real: 0.0, imag: 1.0 });
            assert_eq!(b.dot(a), Rect { real: 0.0, imag: -1.0 });

            // if the above tests failed to ever exercise one of the two code paths
            // in the SIMD code, this will hopefully exercise both
            let a = KetRef::new(&[1.0, 0.0, 0.0, 0.0, 1.0], &[0.0, 0.0, 0.0, 0.0, 0.0]);
            let b = KetRef::new(&[0.0, 0.0, 0.0, 0.0, 0.0], &[1.0, 0.0, 0.0, 0.0, 1.0]);
            assert_eq!(a.dot(b), Rect { real: 0.0, imag: 2.0 });
            assert_eq!(b.dot(a), Rect { real: 0.0, imag: -2.0 });
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
            rect: [Rect]
            real: [f32]
        }

        impl<'a> KetRef<'a> {
            pub fn dot<K: AsKetRef>(self, other: K) -> Rect {
                let other = other.as_ket_ref();
                assert_eq!(self.abs.len(), other.abs.len());
                let table = PhaseTable::get();
                (0..self.abs.len())
                    .map(|i| (self.at(i).conj() * other.at(i)).to_rect(table))
                    .fold(Rect::zero(), |a, b| a + b)
            }

            pub fn sqnorm(&self) -> f32 {
                self.abs.iter().map(|x| x*x).sum()
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
