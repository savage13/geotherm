
use std::env;
use std::fs::read_to_string;
use toml;
use dimensioned as dim;
use geotherm2::GeothermBuilder as Builder;

fn arange(x0: f64, dx: f64, x1: f64) -> Vec<f64> {
    let mut x = vec![];
    let mut xi = x0;
    while xi <= x1 {
        x.push(xi);
        xi += dx;
    }
    x
}

pub trait PressureAtDepth {
    fn pressure(&self) -> dim::si::Pascal<f64>;
}
impl PressureAtDepth for dim::si::Meter<f64> {
    fn pressure(&self) -> dim::si::Pascal<f64> {
        use dim::si::{MPS2, KGPM3};
        *self * 9.81 * MPS2 * 3300.0 * KGPM3
    }
}

macro_rules! convert {
    ($name:ident, $abbrev:expr, $base:ty, $conv: block) => {
        pub struct $name { v: f64 }
        impl From<$base> for $name {
            fn from(thing: $base) -> Self {
                Self { v: $conv(thing.value_unsafe) }
            }
        }
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{:.4} {}", self.v, $abbrev)
            }
        }
        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{:.4} {}", self.v, $abbrev)
            }
        }
        impl std::ops::Deref for $name {
            type Target = f64;
            fn deref(&self) -> &f64 {
                &self.v
            }
        }
    };
}

convert!(Bar,       "bar",  dim::si::Pascal<f64>, { |x| x / 1e5 } );
convert!(GPa,       "GPa",  dim::si::Pascal<f64>, { |x| x / 1e9 } );
convert!(Kbar,      "kbar", dim::si::Pascal<f64>, { |x| x / 1e8 } );
convert!(Celsius,   "C",    dim::si::Kelvin<f64>, { |x| x - 273.15 } );
convert!(Kilometer, "km",   dim::si::Meter<f64>,  { |x| x / 1e3 } );

pub const KM: dim::si::Meter<f64> = dim::si::SI { value_unsafe: dim::si::M.value_unsafe * 1000.0,
                                                  _marker: std::marker::PhantomData, };

fn main() {
    let args: Vec<String> = env::args().collect();
    let par = read_to_string(&args[1]).unwrap();
    let g : Builder = toml::from_str(&par).unwrap();
    let g = g.build();

    let z : Vec<f64> = args[2].split("/").map(|x| x.parse().unwrap()).collect();

    for zi in arange(z[0],z[1],z[2]).iter() {
        let zi = *zi * KM;
        let t = g.t(zi).unwrap();
        let c = Celsius::from(t);
        let p = zi.pressure();
        let bars = Bar::from(p);
        let gpa = GPa::from(p);
        let km = Kilometer::from(zi);
        println!("{:10.3} {:10.3}    {} {} {}", *bars, t.value_unsafe, km, c, gpa );
    }

}
