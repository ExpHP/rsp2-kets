#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[macro_use]
extern crate lazy_static;
#[cfg(test)]
extern crate rand;
#[cfg(feature = "faster")]
extern crate faster;

pub(crate) mod basis;
pub(crate) mod complex;

// actual public API, employing the "pick your parallel namespace" model
// frequently employed by e.g. parsing libraries when they want to
//  support either bytes or str.

pub use lossless::*;
pub mod lossless {
    //! The items in this module are identical to those
    //! already available in the root module.
    //! This module merely exists to give them a more explicit path.

    // (actually it only exists so I have an excuse to indent these
    //  items similarly to the compact module below >_>)
    pub use ::complex::lossless::Rect;

    pub use ::basis::lossless::basis::Basis;
    pub use ::basis::lossless::basis::Iter as BasisIter;

    pub use ::basis::lossless::ket::Ket;
    pub use ::basis::lossless::ket::KetRef;
    pub use ::basis::lossless::ket::AsKetRef;
    pub use ::basis::lossless::ket::Iter as KetIter;
    pub use ::basis::lossless::ket::IntoIter as KetIntoIter;
}

pub mod compact {
    pub use ::complex::compact::Rect;
    pub use ::complex::compact::Polar;

    pub use ::basis::compact::basis::Basis;
    pub use ::basis::compact::basis::Iter as BasisIter;

    pub use ::basis::compact::ket::Ket;
    pub use ::basis::compact::ket::KetRef;
    pub use ::basis::compact::ket::AsKetRef;
    pub use ::basis::compact::ket::Iter as KetIter;
    pub use ::basis::compact::ket::IntoIter as KetIntoIter;
}
