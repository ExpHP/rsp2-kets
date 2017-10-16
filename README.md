# `rsp2-kets`

```toml
[dependencies.rsp2-kets]
git = "https://github.com/ExpHP/rsp2-kets"
version = "0.3"
features = ["serde"]
```

Datatypes representing kets and collections thereof, stored as blobs of `f64`, used in my personal projects.

No tests.  Little documentation.  Not really intended for public consumption.  This would never even have seen the light of day were it not for the fact that I have more than one separate project that uses it.

The `compressed` stuff is optimized for use with [`bincode`](https://github.com/TyOverby/bincode) and [`flate2`](https://github.com/alexcrichton/flate2-rs).
