
use std::env;
use dimensioned as dim;
use geotherm::Geotherm;
use geotherm::convert;

convert!(Bar,       "bar",  dim::si::Pascal<f64>, { |x| x / 1e5 } );
convert!(GPa,       "GPa",  dim::si::Pascal<f64>, { |x| x / 1e9 } );
convert!(Kbar,      "kbar", dim::si::Pascal<f64>, { |x| x / 1e8 } );
convert!(Celsius,   "C",    dim::si::Kelvin<f64>, { |x| x - 273.15 } );
convert!(Kilometer, "km",   dim::si::Meter<f64>,  { |x| x / 1e3 } );

fn arange(x0: f64, dx: f64, x1: f64) -> Vec<f64> {
    let mut x = vec![];
    let mut xi = x0;
    while xi <= x1 {
        x.push(xi);
        xi += dx;
    }
    x
}

trait PressureAtDepth {
    fn pressure(&self) -> dim::si::Pascal<f64>;
}
impl PressureAtDepth for dim::si::Meter<f64> {
    fn pressure(&self) -> dim::si::Pascal<f64> {
        use dim::si::{MPS2, KGPM3};
        *self * 9.81 * MPS2 * 3300.0 * KGPM3
    }
}

fn main() -> Result<(), geotherm::GeothermError> {
    let args: Vec<String> = env::args().collect();
    let g = Geotherm::from_file(&args[1])?;

    let z : Vec<f64> = args[2].split("/").map(|x| x.parse().unwrap()).collect();

    for &zi in arange(z[0],z[1],z[2]).iter() {
        let zi = zi * 1e3 * dim::si::M;
        let t = g.t(zi)?;

        let p    = zi.pressure();
        let bars = Bar::from(p);
        let gpa  = GPa::from(p);
        let c    = Celsius::from(t);
        let km   = Kilometer::from(zi);
        println!("{:10.3} {:10.3}    {:.2} {:.2} {:.3}", bars.val(), t.value_unsafe, km, c, gpa );
    }
    Ok(())
}
