# `rsp2-kets`

```toml
[dependencies.rsp2-kets]
git = "https://github.com/ExpHP/rsp2-kets"
version = "0.1"
features = ["serde"]
```

Datatypes representing kets and collections thereof, stored as blobs of `f64`, used in my personal projects.

No tests.  Unusual public module layout.  Little documentation.  Not really intended for public consumption.

### Tips

* You can't deserialize `Basis` directly. Deserialize a `Raw` type and use `.validate()`.
* The `compressed` stuff is optimized for use with [`bincode`](https://github.com/TyOverby/bincode) and [`flate2`](https://github.com/alexcrichton/flate2-rs).
