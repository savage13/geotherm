//! Geotherm Calculator for the Crust and Mantle
//!
//! ```rust
//! use dimensioned::si::{M,K};
//!
//! let g = geotherm::Geotherm::from_file("test/oceanic_100Ma.toml").unwrap();
//! assert_eq!(g.t(0.0 * M).unwrap(), 273.15 * K);
//! assert_eq!(g.t(100e3 * M).unwrap(), 1328.8354705757051 * K);
//!
//! let ad = geotherm::Adiabat::from_tp((1300.0 + 273.15) * K);
//! assert_eq!(ad.t(0.0 * M).unwrap(),   (1300.0 + 273.15) * K);
//! assert_eq!(ad.t(100e3 * M).unwrap(),  1606.290193479815 * K);
//! ```
//!
//! ## Refernces
//!  - Chapman, D. (1986), Thermal gradients in the continental crust, Geological Society, Lon- don, Special Publications, 24(1), 63–70.
//!  - Faul, U. H., and I. Jackson (2005), The seismological signature of temperature and grain size variations in the upper mantle, Earth Planet. Sci. Lett., 234(1-2), 119–134.
//!  - Jaupart, C., and J. Mareschal (1999), The thermal structure and thickness of continental roots, in Developments in Geotectonics, vol. 24, pp. 93–114, Elsevier.
//!  - Schatz, J. F., and G. Simmons (1972), Thermal conductivity of earth materials at high temperatures, Journal of Geophysical Research, 77(35), 6966–6983.


#![allow(non_snake_case)]
#![allow(unused_imports)]

#[macro_use]
extern crate dimensioned as dim;

use serde::Deserialize;

pub mod erf;
use erf::ErrorFunction;

/// Create a simple Unit Newtype for doing Unit Conversions
///
/// Allows simple unit conversion of units of the same type.
///   This is really only useful for displaying a value in alternative
///   units.  The underlying base unit/type is wrapped in a NewType and only converted
///   when the value is requested either using Display, Debug, or directly using
///   `val`.  The underlying base unit can be accessed through `deref`.
///
/// ## Arguments
///  - New Unit Type
///  - Unit Abbreviation
///  - Underlying Base Unit
///  - Converson from base unit to New Unit Type
///
/// ## Example
///
/// ```rust
/// use dimensioned as dim;
/// use dim::si::{K,Kelvin};
/// use geotherm::convert;
///
/// convert!(Celsius, "C", Kelvin<f64>, { |x| x - 273.15 } );
///
/// let t_in_k = 273.15 * K;
/// let t_in_c = Celsius::from(t_in_k);
/// assert_eq!(t_in_c.val(), 0.0); // Get value
/// assert_eq!(t_in_k, *t_in_c);   // Get base type using deref
/// ```
#[macro_export]
macro_rules! convert {
    ($name:ident, $abbrev:expr, $base:ty, $conv:block) => {
        pub struct $name ( $base );

        impl $name {
            #[inline]
            pub fn val(&self) -> f64 {
                let v = $conv(self.0.value_unsafe);
                v
            }
        }
        impl From<$base> for $name {
            #[inline]
            fn from(thing: $base) -> Self {
                $name ( thing )
            }
        }
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                self.val().fmt(f)?;
                write!(f, " {}", $abbrev)
            }
        }
        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                self.val().fmt(f)?;
                write!(f, " {}", $abbrev)
            }
        }
        impl std::ops::Deref for $name {
            type Target = $base;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }
    };
}

macro_rules! C {
    ($v: expr) => {{
        Celsius::new($v)
    }};
}
/// A Year in Seconds
pub const YEAR: Second<f64> = SI { value_unsafe: dim::si::DAY.value_unsafe * 365.24,
                                   _marker: std::marker::PhantomData, };
/// A Billon Years in Seconds
pub const GA: Second<f64> = SI { value_unsafe: YEAR.value_unsafe * 1e9,
                                   _marker: std::marker::PhantomData, };
/// A Millon Years in Seconds
pub const MA: Second<f64> = SI { value_unsafe: YEAR.value_unsafe * 1e6,
                                   _marker: std::marker::PhantomData, };


use dim::si::{self,SI,Meter,Kelvin,Second};

derived!(si, SI: HeatGeneration = Watt / Meter3);
derived!(si, SI: Conductivity = Watt / Meter / Kelvin);
derived!(si, SI: HeatFlux = Watt / Meter2 );
derived!(si, SI: ThermalDiffusivity = Meter2 / Second );
derived!(si, SI: ThermalGradient = Kelvin / Meter );
derived!(si, SI: Acceleration = Meter / Second2 );
derived!(si, SI: CoefficientThermalExpansion = Unitless / Kelvin );
derived!(si, SI: SpecificHeat = Joule / Kilogram / Kelvin );

#[derive(Debug)]
enum KType {
    Constant(Conductivity<f64>),
    Function(fn(Kelvin<f64>, Meter<f64>) -> Conductivity<f64>),
}

#[derive(Debug)]
struct ContinentalLayer {
    pub zt : Meter<f64>,
    pub zb : Meter<f64>,
    pub rh : HeatGeneration<f64>,
    //pub k  : Conductivity<f64>,
    pub k  : KType,

    pub qt : HeatFlux<f64>,
    pub qb : HeatFlux<f64>,
    pub tt : Kelvin<f64>,
    pub tb : Kelvin<f64>,
}

/// Continental Geotherm Calculator
///
/// ## Example
///
/// ```rust
/// use dimensioned::si::{M,W,K,M2,M3};
/// let km = 1e3 * M;
///
/// let kcrust    = 2.5 * W/(M*K);
/// let kmantle   = 3.0 * W/(M*K);
/// let t_surface = 273.15 * K;
/// let q_surface = 79e-3 * W/M2;
///
/// let mut g = geotherm::ContinentalGeotherm::new(t_surface, q_surface);
/// g.add_layer(11.2  * km, kcrust,  3.00e-6 * W/M3);
/// g.add_layer(27.8  * km, kcrust,  1.00e-6 * W/M3);
/// g.add_layer(211.0 * km, kmantle, 0.03e-6 * W/M3);
///
/// assert_eq!(g.t(0.0   * M).unwrap(), 273.15 * K);
/// assert_eq!(g.t(100e3 * M).unwrap(), 1241.3476666666668 * K);
///```
///
/// ## Method
///
/// Construction of a continental geother with depth based on
///  Faul and Jackson (2005) and Chapman (1986)
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><msub><mi>T</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>=</mo><msub><mi>T</mi><mi>i</mi></msub><mo>+</mo><mfrac><msub><mi>q</mi><mi>i</mi></msub><msub><mi>k</mi><mi>i</mi></msub></mfrac><mi>&#x394;</mi><mi>z</mi><mo>-</mo><mfrac><mi>A</mi><mrow><mn>2</mn><msub><mi>k</mi><mi>i</mi></msub></mrow></mfrac><mi>&#x394;</mi><msup><mi>z</mi><mn>2</mn></msup></math>
///
/// and
///  <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><msub><mi>q</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>=</mo><msub><mi>q</mi><mi>i</mi></msub><mo>-</mo><mi>A</mi><mi>&#x394;</mi><mi>z</mi></math>
///
/// where
///
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>T</mi><mi>i</mi></msub></math> is the Temperature at the top a layer
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>T</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub></math> is the Temperature at the bottom of the layer
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>q</mi><mi>i</mi></msub></math> is the Heat flow at the top of the layer
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>q</mi><mrow><mi>i</mi><mo>+</mo><mn>1</mn></mrow></msub></math> is the heat flow at the bottom of the layer
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>k</mi><mi>i</mi></msub></math> is the Thermal Conductivity within the layer and is assumed constant
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>A</mi></math> is the Heat Generation within the layer
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>&#x394;</mi><mi>z</mi></math> is the thickness of the layer
/// ## Refernces
///   Chapman, D. (1986), Thermal gradients in the continental crust, Geological Society, Lon- don, Special Publications, 24(1), 63–70.
///
///   Faul, U. H., and I. Jackson (2005), The seismological signature of temperature and grain size variations in the upper mantle, Earth Planet. Sci. Lett., 234(1-2), 119–134.
///
#[derive(Debug)]
pub struct ContinentalGeotherm {
    /// Continental Layers to compute the geotherm
    layers: Vec<ContinentalLayer>,
    /// Surface Heat floe
    heat_flow: HeatFlux<f64>,
    /// Surface Temperature
    ts: Kelvin<f64>,
    /// Mantle Potential Temperature
    tp: Kelvin<f64>,
}

//pub fn q_in_layer(z: Length, zt: Length, qt: HeatFlux, rh: HeatGeneration) -> HeatFlux {
//    qt - rh * (z-zt)
//}

impl ContinentalLayer {
    /// Create a new Continental Layer between the Top and bottom depths `zt`
    /// and `zb`, with Conductivity `k` and Internal Heat Generation `rh`
    ///
    pub fn new<KT: Into<Kelvin<f64>>>(zt: Meter<f64>, zb: Meter<f64>,
                                      k: KType, rh: HeatGeneration<f64>,
                                      tt: KT, qt: HeatFlux<f64> ) -> Self {
        use dim::si::{WPM2,K};
        let mut s = Self {
            zt, zb, k, rh,
            qt, qb: 0.0 * WPM2,
            tt: tt.into(), tb: 0.0 * K,
        };
        s.tb = s.t(s.zb).unwrap();
        s.qb = s.q(s.zb).unwrap();
        s
    }
    /// Compnute the Temperature at a Depth `z`
    ///
    /// Following Faul and Jackson, 2005, Eqn 18.
    ///
    /// <!-- begin MathToWeb --><!-- (your LaTeX) $T_b = T_t + \dfrac{q_t}{k}\Delta z - \dfrac{\rho H}{2k}\Delta z^2$ --><math xmlns="http://www.w3.org/1998/Math/MathML"><mrow><msub><mi>T</mi><mi>b</mi></msub><mo>=</mo><msub><mi>T</mi><mi>t</mi></msub><mo>+</mo><mstyle scriptlevel="0" displaystyle="true"><mrow><mfrac linethickness="1"><mrow><msub><mi>q</mi><mi>t</mi></msub></mrow><mi>k</mi></mfrac></mrow></mstyle><mi>&#x00394;</mi><mi>z</mi><mo>-</mo><mstyle scriptlevel="0" displaystyle="true"><mrow><mfrac linethickness="1"><mrow><mi>&#x003C1;</mi><mi>H</mi></mrow><mrow><mn>2</mn><mi>k</mi></mrow></mfrac></mrow></mstyle><mi>&#x00394;</mi><msup><mi>z</mi><mn>2</mn></msup></mrow></math><!-- end MathToWeb -->
    ///
    pub fn t(&self, z: Meter<f64>) -> Result<Kelvin<f64>,GeothermError> {
        use dim::si::WPMK;
        if ! self.within(z) {
            return Err(GeothermError::InvalidDepth);
        }
        let tt = self.tt;
        let qt = self.qt;
        let z0 = z - self.zt;
        let kv = match self.k {
            KType::Constant(v) => v,
            KType::Function(fun) => fun(tt, z),
        };
        Ok( tt + (qt * z0) / kv - self.rh * z0 * z0/ (2.0*kv) )
    }
    /// Determine with depth `z` is within the layer
    fn within(&self, z: Meter<f64>) -> bool {
        self.zt <= z && z <= self.zb
    }
    /// Computes the Heat flow at the bottom of the layer
    ///
    /// After Faul and Jackson, 2005, Eqn. 18 (part 2)
    ///
    /// <!-- begin MathToWeb --><!-- (your LaTeX) $q_b = q_t - \rho H \Delta z$ --><math xmlns="http://www.w3.org/1998/Math/MathML"><mrow><msub><mi>q</mi><mi>b</mi></msub><mo>=</mo><msub><mi>q</mi><mi>t</mi></msub><mo>-</mo><mi>&#x003C1;</mi><mi>H</mi><mi>&#x00394;</mi><mi>z</mi></mrow></math><!-- end MathToWeb -->
    ///
    pub fn q(&self, z: Meter<f64>) -> Result<HeatFlux<f64>, GeothermError> {
        if ! self.within(z) {
            return Err(GeothermError::InvalidDepth);
        }
        Ok( self.qt - self.rh * (z - self.zt) )
    }
    /// Set the Top Temperature `tt` and Heat Flux, `qt`
    ///
    ///
    ///
    fn _set_top_t_q<KT: Into<Kelvin<f64>>>(&mut self, tt: KT, qt: HeatFlux<f64>) {
        self.tt = tt.into();
        self.qt = qt;
        self.tb = self.t(self.zb).unwrap();
        self.qb = self.q(self.zb).unwrap();
    }
}

impl ContinentalGeotherm {
    /// Create an empty Continental Geotherm with surface temperature `ts`
    ///  and a surface `heat_flow`
    ///
    /// To properly compute the temperature within the geotherm
    /// layers need to be added using either `add_layer` or `add_layer_kfun`
    ///
    pub fn new<KT: Into<Kelvin<f64>>>(ts: KT, heat_flow: HeatFlux<f64>) -> Self {
        Self { layers: vec![], heat_flow, ts: ts.into(), tp: C!(1300.0).into(), }
    }
    /// Set the Mantle Potential Temperature associated with the Geotherm
    ///
    /// Temperatures exceeding the `Adiabat`, set by the Mantle Potential Temperature
    ///    will be set to the Adiabatic Temperature at that depth, or
    ///
    /// `T = min(T_geotherm, T_adiabat)`
    pub fn set_tp<KT: Into<Kelvin<f64>>>(&mut self, tp: KT) {
        self.tp = tp.into();
    }
    /// Get a reference to the ContinentalLayer which contains depth `z`
    fn z_to_layer<'a> (&'_ self, z: Meter<f64>) -> Result<&'_ ContinentalLayer, GeothermError> {
        for layer in self.layers.iter() {
            if layer.within(z) {
                return Ok(layer)
            }
        }
        Err(GeothermError::InvalidDepth)
    }
    /// Add a Layer with a `thickness`, conductivity `k` and Internal Heat
    ///   Generation `rH`.
    ///
    /// Conductivity within the layer is constant.
    pub fn add_layer(&mut self, thickness: Meter<f64>, k: Conductivity<f64>, rH: HeatGeneration<f64>) {
        let zb = if let Some(bot) = self.layers.last() { bot.zb } else { 0.0 * dim::si::M  };
        let layer =
        if let Some(last) = self.layers.last() {
            ContinentalLayer::new(zb, zb + thickness, KType::Constant(k), rH, last.tb, last.qb)
        } else {
            ContinentalLayer::new(zb, zb + thickness, KType::Constant(k), rH, self.ts, self.heat_flow)
        };
        self.push(layer);
    }
    /// Add a Layer with a `thickness`,  Internal Heat Generation `rH`,
    ///   and Conductivity defined by a function
    ///
    /// Conductivity within the layer is constant. To produce a gradient in
    ///   conductivity, multiple thin layers can be added
    pub fn add_layer_kfun(&mut self, thickness: Meter<f64>,
                          rH: HeatGeneration<f64>,
                          kfun: fn(Kelvin<f64>, Meter<f64>) -> Conductivity<f64>) {
        let zb = if let Some(bot) = self.layers.last() { bot.zb } else { 0.0 * dim::si::M  };
        let layer =
        if let Some(last) = self.layers.last() {
            ContinentalLayer::new(zb, zb + thickness, KType::Function(kfun), rH, last.tb, last.qb)
        } else {
            ContinentalLayer::new(zb, zb + thickness, KType::Function(kfun), rH, self.ts, self.heat_flow)
        };
        self.push(layer);
    }

    /// Add a Layer to the bottom of the ContinentalGeotherm
    ///
    /// Typically layers are added using `add_layer` or `add_layer_kfun`
    ///
    /// If layers exist and the bottom layer has the bottom temperature and
    /// heat flow defined, set the top temperature and heat flow of the
    /// added layer to match the layer above (boundary condition) and compute
    /// the new layers bottom temperature and heat flow
    ///
    fn push(&mut self, layer: ContinentalLayer) {
        self.layers.push(layer);
    }
    /// Compute the temperature at of the geotherm at a depth `z`
    ///
    /// Depth, `z`, values outside the defined range return `None`
    ///  otherwise the Temperature is returned in <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>K</mi></math>
    ///
    pub fn t(&self, z: Meter<f64>) -> Result<Kelvin<f64>, GeothermError> {
        let layer = self.z_to_layer(z)?;
        let abat = Adiabat::from_tp(self.tp);
        let ta = abat.t(z)?;
        let tg = layer.t(z)?;
        let tg = if tg > ta { ta } else { tg };
        Ok(tg)
    }
    /// Compute the Heat Flow at of the geotherm at a depth `z`
    ///
    /// Depth, `z`, values outside the defined range return `None`
    ///  otherwise the Heat Flow is returned in <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>W</mi><msup><mi>m</mi><mrow><mo>-</mo><mn>2</mn></mrow></msup></math>
    ///
    pub fn q(&self, z: Meter<f64>) -> Result<HeatFlux<f64>, GeothermError> {
        let layer = self.z_to_layer(z)?;
        layer.q(z)
    }
    /// Compute a full geotherm at a depth spacing of `dz`
    pub fn geotherm(&self, dz: Meter<f64>) -> Option<(Vec<Meter<f64>>,Vec<Kelvin<f64>>)> {
        if let Some(layer) = self.layers.last() {
            let z1 = layer.zb;
            let z0 = 0.0 * dim::si::M;
            let z : Vec<_> = (0..)
                .map(|n| z0 + dz*n as f64)
                .take_while(|&z| z <= z1)
                .collect();
            let t : Vec<_> = z.iter()
                .map(|&d| self.t(d).unwrap())
                .collect();
            return Some((z,t));
        }
        None
    }
}


/// Oceanic Geotherm Calculator
///
/// ## Example
///
/// ```rust
/// use dimensioned::si::{M,K};
///
/// let g = geotherm::Geotherm::from_file("test/oceanic_100Ma.toml").unwrap();
/// assert_eq!(g.t(0.0 * M).unwrap(), 273.15 * K);
/// assert_eq!(g.t(100e3 * M).unwrap(), 1328.8354705757051 * K);
/// ```
///
/// ## Method
///
/// Construction of this geotherm with depth is based on
/// Faul and Jackson, 2005, Eqn. 17
///
/// <!--  $T(z) = T_s + (T_{ad} - T_{s}) \textrm{erf} \left( \dfrac{z}{2\sqrt{\kappa t}}\right)$ -->
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><mi>T</mi><mrow><mo form="prefix">(</mo><mi>z</mi><mo form="postfix">)</mo></mrow><mo>=</mo><msub><mi>T</mi><mi>s</mi></msub><mo>+</mo><mrow><mo form="prefix">(</mo><msub><mi>T</mi><mrow><mi>a</mi><mi>d</mi></mrow></msub><mo>-</mo><msub><mi>T</mi><mi>s</mi></msub><mo form="postfix">)</mo></mrow><mstyle mathvariant="normal"><mi>e</mi><mi>r</mi><mi>f</mi></mstyle><mrow><mo rspace="0.3em" lspace="0em" stretchy="true" fence="true" form="prefix">(</mo><mstyle scriptlevel="0" displaystyle="true"><mrow><mfrac linethickness="1"><mi>z</mi><mrow><mn>2</mn><msqrt><mrow><mi>&#x003BA;</mi><mi>t</mi></mrow></msqrt></mrow></mfrac></mrow></mstyle><mo rspace="0em" lspace="0.3em" stretchy="true" fence="true" form="postfix">)</mo></mrow></mrow></math>
///
/// where
///
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>T</mi><mi>s</mi></msub></math> is the surface Temperature
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>T</mi><mrow><mi>a</mi><mi>d</mi></mrow></msub></math> is the Mantle Adiabatic Temperature
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>z</mi></math> is depth in meters
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>&#x3BA;</mi></math> is the Thermal Diffusivity
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>t</mi></math> is the age in seconds
///
///
/// ## References
///   Faul, U. H., and I. Jackson (2005), The seismological signature of temperature and grain size variations in the upper mantle, Earth Planet. Sci. Lett., 234(1-2), 119–134.
///
#[derive(Debug)]
pub struct OceanicGeotherm {
    /// Thermal Diffusivity m^2/s
    kappa: ThermalDiffusivity<f64>,
    /// Surface Temperature
    ts: Kelvin<f64>,
    /// Mantle Adiabat
    ad: Adiabat,
    /// Age of the Oceanic Lithosphere
    age: Second<f64>,
}

impl OceanicGeotherm {
    /// Create a new Oceanic Geotherm from the Thermal Diffusivity, `kappa`,
    /// `age`, Surface Temperature, `ts` and Adibat `ad`
    pub fn new<KT: Into<Kelvin<f64>>>(ad: Adiabat, kappa: ThermalDiffusivity<f64>, ts: KT, age: Second<f64>) -> Self {
        OceanicGeotherm { ad, kappa, ts: ts.into(), age }
    }
    /// Create a new Oceanic Geotherm from a Mantle Potential Temperature `tp`
    /// and an `age` of the lithosphere
    ///
    /// Surface Temperature, <math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>T</mi><mi>s</mi></msub></math>, is set to <math xmlns="http://www.w3.org/1998/Math/MathML"><mn>0</mn><mo>&#xB0;</mo><mi>C</mi></math>
    ///
    /// Thermal Diffusivity, <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>&#x3BA;</mi></math>, is set at <math xmlns="http://www.w3.org/1998/Math/MathML"><mn>10</mn><mo>.</mo><mn>0</mn><mo>&#xD7;</mo><msup><mn>10</mn><mrow><mo>-</mo><mn>7</mn></mrow></msup><msup><mi>m</mi><mn>2</mn></msup><msup><mi>s</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup></math>.
    ///
    /// Diffusivity value is derived from the lower temperature estimate of Gibert et al., 2003, Fig. 6.
    ///
    pub fn from_tp_and_age<KT: Into<Kelvin<f64>>>(tp: KT, age: Second<f64>) -> Self {
        use dim::si::M2PS;
        let ad = Adiabat::from_tp(tp.into());
        let kappa = 10.0e-7 * M2PS; // Thermal diffusivity
        let ts = C!(0.0).into();
        OceanicGeotherm { ad, kappa, age, ts }
    }
    /// Compute the Temperature at a depth for an Oceanic Geotherm
    ///
    pub fn t(&self, z: Meter<f64>) -> Result<Kelvin<f64>, GeothermError> {
        use dim::Sqrt;
        if z < 0.0 * dim::si::M {
            return Err(GeothermError::InvalidDepth);
        }
        let v = z / (2.0 * (self.kappa * self.age).sqrt() );
        let vo = if v < 26.0 * dim::si::ONE { v.value_unsafe.erf() } else { 1.0 };
        let dt = (self.ad.t(z)? - self.ts) * vo;
        let t = self.ts + dt;
        Ok(t)
    }
}

/// Adiabatic Temperature Gradient for Earth's Mantle
///
/// ## Example
///
/// ```rust
/// use dimensioned::si::{M,K};
///
/// let ad = geotherm::Adiabat::from_tp((1300.0 + 273.15) * K);
/// assert_eq!(ad.t(0.0 * M).unwrap(),   (1300.0 + 273.15) * K);
/// assert_eq!(ad.t(100e3 * M).unwrap(),  1606.290193479815 * K);
/// ```
///
/// ## Method
///
/// This follows from Faul and Jackson, 2005 Eqn. 15
///
/// <!-- $\dfrac{dT}{dz} = \dfrac{T_p * \alpha * g}{c_p}$ -->
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><mstyle scriptlevel="0" displaystyle="true" display="block"><mrow><mfrac linethickness="1"><mrow><mi>dT</mi></mrow><mrow><mi>dz</mi></mrow></mfrac></mrow></mstyle><mo>=</mo><mstyle scriptlevel="0" displaystyle="true"><mrow><mfrac linethickness="1"><mrow><msub><mi>T</mi><mi>p</mi></msub><mi>&#x003B1;</mi><mi>g</mi></mrow><mrow><msub><mi>c</mi><mi>p</mi></msub></mrow></mfrac></mrow></mstyle></mrow></math>
///
/// ## References
///  - Faul, U. H., and I. Jackson (2005), The seismological signature of temperature and grain size variations in the upper mantle, Earth Planet. Sci. Lett., 234(1-2), 119–134.
#[derive(Debug,Copy,Clone)]
pub struct Adiabat {
    /// Mantle Potential Temperature (K)
    tp: Kelvin<f64>,
    /// Thermal Gradient with Depth (K/m)
    dtdz: ThermalGradient<f64>,
}

impl Adiabat {
    /// Create new Adiabat from the Mantle Potential Temperature `tp` in Kelvin and
    /// the `dtdz` thermal gradient in Kelvin/meter
    pub fn new<KT: Into<Kelvin<f64>>>(tp: KT, dtdz: ThermalGradient<f64>) -> Adiabat {
        Adiabat { tp: tp.into(), dtdz }
    }
    /// Compute Adiabat from the Mantle Potential Temperature `tp`, gravity `g`
    /// Coefficient of Thermal Expansion `alpha` and Specific Heat `cp`
    ///
    pub fn from_values<KT: Into<Kelvin<f64>>>(tp: KT,
                                              g: Acceleration<f64>,
                                              alpha: CoefficientThermalExpansion<f64>,
                                              cp: SpecificHeat<f64>) -> Adiabat {
        let tp = tp.into();
        let dtdz = tp * alpha * g / cp;
        Adiabat::new(tp, dtdz)
    }
    /// Compute an Adiabat from a Potential Temperautre `tp` with "reasonable"
    /// assumptions to compute the thermal gradient
    ///
    /// - g: 9.8 m/s^2 (Gravity)
    /// - &#x003B1;: 2.9e-5 1/K (Coefficient of Thermal Expansion)
    /// - cp: 1350 J/ (kg K) (Specific Heat)
    ///
    pub fn from_tp<KT: Into<Kelvin<f64>>>(tp: KT) -> Adiabat {
        use dim::si::{M,S2,K,J,KG};
        Adiabat::from_values(tp.into(),           // Potential Temperature
                             9.80665 * M/S2,      // Gravity
                             2.9e-5 / K,          // Coeff. Thermal Expansion
                             1350.0 * J / KG / K) // Specific Heat
    }
    /// Compute the Adiabatic Temperature at a depth `z` in meters
    pub fn t(&self, z: Meter<f64>) -> Result<Kelvin<f64>,GeothermError> {
        if z < 0.0 * dim::si::M {
            Err(GeothermError::InvalidDepth)
        } else {
            Ok( self.tp + self.dtdz * z )
        }
    }
}

/// Temperature dependent thermal conductivity
///
/// ## Arguments
/// - t : Temperature in Kelvin
/// - z : Depth in meters
///
/// ## Method
/// From Jaupart and Mareschal 1999, Eqn 3
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>k</mi><mo>(</mo><mi>T</mi><mo>)</mo><mo>=</mo><mfrac><mn>1</mn><mrow><mi>a</mi><mo>+</mo><mi>b</mi><mi>T</mi></mrow></mfrac><mo>+</mo><mi>c</mi><msup><mi>T</mi><mn>3</mn></msup></math>
///
/// where
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>a</mi></math> = 0.174
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>b</mi></math> = 0.000265
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>c</mi></math> = <math xmlns="http://www.w3.org/1998/Math/MathML"><mn>3</mn><mo>.</mo><mn>68</mn><mo>&#xD7;</mo><msup><mn>10</mn><mrow><mo>-</mo><mn>10</mn></mrow></msup></math>
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>T</mi></math> is the temperature measured in Kelvin
///
/// The depth value `z` is unused in this formulation
///
/// ## References
///   Jaupart, C., and J. Mareschal (1999), The thermal structure and thickness of continental roots, in Developments in Geotectonics, vol. 24, pp. 93–114, Elsevier.
pub fn kappa_mantle_jm99(t: Kelvin<f64>, _z: Meter<f64>) -> Conductivity<f64> {
    let t = t.value_unsafe;
    let k = 1.0 / (0.174 + 0.000265 * t) + 3.68e-10 * t.powi(3);
    k * dim::si::WPMK
}
/// Thermal Conductivity formulation from Chapman (1986)
///
/// ## Arguments
/// - t : Temperature in Kelvin
/// - z : Depth in meters
///
/// ## Reference
/// Chapman, D. Thermal gradients in the continental crust. Geological Society,
///    London, Special Publications, 24(1):63–70, 1986.
///
fn kcrust(t: Kelvin<f64>, z: Meter<f64>, k0: f64, b: f64) -> Conductivity<f64> {
    use dim::si::WPMK;
    let z = z.value_unsafe / 1e3;
    let t = t.value_unsafe - 273.15;
    let k = k0.powf(1.0 + 1.5e-3 * z) / (1.0 * b * t);
    k * WPMK
}

/// Thermal conductivity calculation, relevant to the upper crust
///
/// ## Arguments
/// - t : Temperature in Kelvin
/// - z : Depth in meters
///
/// ## Method
/// Conductivity is calculated using the equation below where `t` Temperature is in Celcius and
/// `z` is in meters.  Values are automatically converted to the correct units
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" alttext="k(T,z)=k_{0}\dfrac{1+cz}{1+bz}" display="block">
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="k(T,z)=k_{0}\frac{1+cz}{1+bz}" display="block">
/// <mrow> <mrow> <mi>k</mi> <mo>⁢</mo> <mrow> <mo stretchy="false">(</mo> <mi>T</mi> <mo>,</mo> <mi>z</mi> <mo stretchy="false">)</mo> </mrow> </mrow> <mo>=</mo> <mrow> <msub> <mi>k</mi> <mn>0</mn> </msub> <mo>⁢</mo> <mfrac> <mrow> <mn>1</mn> <mo>+</mo> <mrow> <mi>c</mi> <mo>⁢</mo> <mi>z</mi> </mrow> </mrow> <mrow> <mn>1</mn> <mo>+</mo> <mrow> <mi>b</mi> <mo>⁢</mo> <mi>T</mi> </mrow> </mrow> </mfrac> </mrow> </mrow></math></math>
///
/// where
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="k_{0}=3.0\;Wm^{-1}K^{-1}" display="inline"> <mrow> <msub> <mi>k</mi> <mn>0</mn> </msub> <mo>=</mo> <mrow> <mpadded width="+2.8pt"> <mn>3.0</mn> </mpadded> <mo>⁢</mo> <mi>W</mi> <mo>⁢</mo> <msup> <mi>m</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> <mo>⁢</mo> <msup> <mi>K</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> </mrow> </mrow></math></math><br/>
///
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="b=1.5\times 10^{-3}\;K^{-1}" display="inline"> <mrow> <mi>b</mi> <mo>=</mo> <mrow> <mrow> <mn>1.5</mn> <mo>×</mo> <mpadded width="+2.8pt"> <msup> <mn>10</mn> <mrow> <mo>-</mo> <mn>3</mn> </mrow> </msup> </mpadded> </mrow> <mo>⁢</mo> <msup> <mi>K</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> </mrow> </mrow></math></math><br/>
///
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="c=1.5\times 10^{-3}\;km^{-1}" display="inline"> <mrow> <mi>c</mi> <mo>=</mo> <mrow> <mrow> <mn>1.5</mn> <mo>×</mo> <mpadded width="+2.8pt"> <msup> <mn>10</mn> <mrow> <mo>-</mo> <mn>3</mn> </mrow> </msup> </mpadded> </mrow> <mo>⁢</mo> <mi>k</mi> <mo>⁢</mo> <msup> <mi>m</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> </mrow> </mrow></math></math><br/>
///
///
/// ## Reference
/// Chapman, D. Thermal gradients in the continental crust. Geological Society,
///    London, Special Publications, 24(1):63–70, 1986.
///
pub fn kappa_upper_crust_c86(t: Kelvin<f64>, z: Meter<f64>) -> Conductivity<f64> {
    kcrust(t, z, 3.0, 1.5e-3)
}

/// Thermal conductivity calculation, relevant to the lower crust
///
/// ## Arguments
/// - t : Temperature in Kelvin
/// - z : Depth in meters
///
/// ## Method
/// Conductivity is calculated using the equation below where `t` Temperature is in Celcius and
/// `z` is in meters.  Values are automatically converted to the correct units
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" alttext="k(T,z)=k_{0}\dfrac{1+cz}{1+bz}" display="block">
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="k(T,z)=k_{0}\frac{1+cz}{1+bz}" display="block">
/// <mrow> <mrow> <mi>k</mi> <mo>⁢</mo> <mrow> <mo stretchy="false">(</mo> <mi>T</mi> <mo>,</mo> <mi>z</mi> <mo stretchy="false">)</mo> </mrow> </mrow> <mo>=</mo> <mrow> <msub> <mi>k</mi> <mn>0</mn> </msub> <mo>⁢</mo> <mfrac> <mrow> <mn>1</mn> <mo>+</mo> <mrow> <mi>c</mi> <mo>⁢</mo> <mi>z</mi> </mrow> </mrow> <mrow> <mn>1</mn> <mo>+</mo> <mrow> <mi>b</mi> <mo>⁢</mo> <mi>T</mi> </mrow> </mrow> </mfrac> </mrow> </mrow></math></math>
///
/// where
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="k_{0}=2.6\;Wm^{-1}K^{-1}" display="inline"> <mrow> <msub> <mi>k</mi> <mn>0</mn> </msub> <mo>=</mo> <mrow> <mpadded width="+2.8pt"> <mn>2.6</mn> </mpadded> <mo>⁢</mo> <mi>W</mi> <mo>⁢</mo> <msup> <mi>m</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> <mo>⁢</mo> <msup> <mi>K</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> </mrow> </mrow></math></math><br/>
///
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="b=1.0\times 10^{-4}\;K^{-1}" display="inline"> <mrow> <mi>b</mi> <mo>=</mo> <mrow> <mrow> <mn>1.5</mn> <mo>×</mo> <mpadded width="+2.8pt"> <msup> <mn>10</mn> <mrow> <mo>-</mo> <mn>3</mn> </mrow> </msup> </mpadded> </mrow> <mo>⁢</mo> <msup> <mi>K</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> </mrow> </mrow></math></math><br/>
///
/// <?xml version="1.0" encoding="UTF-8"?><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="c=1.5\times 10^{-3}\;km^{-1}" display="inline"> <mrow> <mi>c</mi> <mo>=</mo> <mrow> <mrow> <mn>1.5</mn> <mo>×</mo> <mpadded width="+2.8pt"> <msup> <mn>10</mn> <mrow> <mo>-</mo> <mn>3</mn> </mrow> </msup> </mpadded> </mrow> <mo>⁢</mo> <mi>k</mi> <mo>⁢</mo> <msup> <mi>m</mi> <mrow> <mo>-</mo> <mn>1</mn> </mrow> </msup> </mrow> </mrow></math></math><br/>
///
/// 
/// ## Reference
/// Chapman, D. Thermal gradients in the continental crust. Geological Society,
///    London, Special Publications, 24(1):63–70, 1986.
///
pub fn kappa_lower_crust_c86(t: Kelvin<f64>, z: Meter<f64>) -> Conductivity<f64> {
    kcrust(t, z, 2.6, 1.0e-4)
}
///
/// Calculation of Thermal Conductivity of the Mantle determined from T and z
///
/// ## Arguments
/// - t : Temperature in Kelvin
/// - z : Depth in meters
///
/// ## Method
///
/// Uses Schatz and Simmons (1972) Eqs 9, 10, and 11
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><msub><mi>K</mi><mi>L</mi></msub><mo>=</mo><mfrac><mn>1</mn><mrow><mo>(</mo><mi>a</mi><mo>+</mo><mi>b</mi><mi>T</mi><mo>)</mo></mrow></mfrac></math>
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><msub><mi>K</mi><mi>R</mi></msub><mo>&#xA0;</mo><mo>=</mo><mo>&#xA0;</mo><mfenced open="{" close=""><mtable columnspacing="1.4ex" columnalign="left"><mtr><mtd><mn>0</mn></mtd><mtd><mi>T</mi><mo>&lt;</mo><mn>500</mn><mi>K</mi></mtd></mtr><mtr><mtd><mi>c</mi><mo>(</mo><mn>500</mn><mo>-</mo><mi>T</mi><mo>)</mo></mtd><mtd><mi>T</mi><mo>&#x2265;</mo><mn>500</mn><mi>K</mi></mtd></mtr></mtable></mfenced></math>
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><msub><msub><mi>K</mi><mi>L</mi></msub><mrow><mi>m</mi><mi>i</mi><mi>n</mi></mrow></msub><mo>=</mo><mi>d</mi><mo>(</mo><mn>1</mn><mo>+</mo><mi>z</mi><mo>)</mo></math>
///
///<math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>k</mi><mo>(</mo><mi>T</mi><mo>,</mo><mi>z</mi><mo>)</mo><mo>=</mo><msub><mi>K</mi><mi>R</mi></msub><mo>+</mo><mtext>max</mtext><mo>(</mo><msub><mi>K</mi><mi>L</mi></msub><mo>,</mo><msub><msub><mi>K</mi><mi>L</mi></msub><mrow><mi>m</mi><mi>i</mi><mi>n</mi></mrow></msub><mo>)</mo></math>
///
/// where
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>a</mi><mo>=</mo><mn>7</mn><mo>.</mo><mn>409177</mn><mo>&#xD7;</mo><msup><mn>10</mn><mrow><mo>-</mo><mn>2</mn></mrow></msup><mo>&#xA0;</mo><mi>m</mi><mi>K</mi><mo>/</mo><mi>W</mi></math>
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>b</mi><mo>=</mo><mn>5</mn><mo>.</mo><mn>0191204</mn><mo>&#xD7;</mo><msup><mn>10</mn><mrow><mo>-</mo><mn>4</mn></mrow></msup><mo>&#xA0;</mo><mi>m</mi><mo>/</mo><mi>W</mi></math>
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>c</mi><mo>=</mo><mn>2</mn><mo>.</mo><mn>3</mn><mo>&#xD7;</mo><msup><mn>10</mn><mrow><mo>-</mo><mn>3</mn></mrow></msup><mo>&#xA0;</mo><mi>W</mi><mo>/</mo><mi>m</mi><msup><mi>K</mi><mn>2</mn></msup></math>
/// - <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>d</mi><mo>=</mo><mn>1</mn><mo>.</mo><mn>255199</mn><mi>W</mi><mo>/</mo><mi>m</mi><mi>K</mi><mo>/</mo><mo>(</mo><mn>1000</mn><mi>k</mi><mi>m</mi><mo>)</mo></math>
///
/// ## References
/// Schatz, J. F., and G. Simmons (1972), Thermal conductivity of earth 
///    materials at high temperatures, Journal of Geophysical Research, 
///    77(35), 6966–6983.
///
pub fn kappa_mantle_ss72(t: Kelvin<f64>, z: Meter<f64>) -> Conductivity<f64> {
    use dim::si::WPMK;
    let t = t.value_unsafe;
    let z = z.value_unsafe;
    // Schatz and Simmons, 1972, Eqn (9) pg 6975
    //    Temperature in Kelvin
    //    31 (calories / s C))^-1 = 7.409177e-2 (W/(m K))^-1
    //    0.21 (calories / s ))^-1 = 5.0191204e-4 (W/m )^-1
    let a = 7.4091778e-2;
    let b = 5.0191204e-4;
    let kl = 1.0/(a + b * t);
    // Schatz and Simmons, 1972, Eqn (10)
    //    Temperature in Kelvin
    //    (W / (m K)) / (calories / s C) = 418.4
    //    5.5e-6 (calories / s C^2 ) = 2.3e-3 (W / m K^2)
    let kr = if t > 500.0 {
        let c = 2.3e-3;
        c * (t - 500.0)
    } else {
        0.0
    };
    // Schatz and Simmons, 1972, Eqn (11)
    //    z in 1000's of km, hence 1e3 for m -> km and 1e3 for 1km -> 1000km
    //    0.003 (calories / s C) = 1.255199 (W / (m K))
    let d = 1.255199;
    let klmin =  d * (1.0 + z /1e6);
    let k = klmin + if kl > kr { kl } else { kr };
    k * WPMK
}

/// Temperature in Celsius
#[derive(Copy,Clone,PartialEq)]
struct Celsius {
    /// Temperature in Celsius
    v: f64
}

impl Celsius {
    pub fn new(v: f64) -> Self {
        Self { v }
    }
}

impl std::ops::Deref for Celsius {
    type Target = f64;
    fn deref(&self) -> &f64 {
        &self.v
    }
}

impl From<Kelvin<f64>> for Celsius {
    fn from(t: Kelvin<f64>) -> Self {
        Celsius { v: t.value_unsafe - 273.15 }
    }
}
impl From<Celsius> for Kelvin<f64> {
    fn from(t: Celsius) -> Self {
        (*t + 273.15) * dim::si::K
    }
}
impl std::fmt::Display for Celsius {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:.4} C", self.v)
    }
}
impl std::fmt::Debug for Celsius {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:.4} C", self.v)
    }
}



#[cfg(test)]
mod tests {
    //use super::*;
    use crate::Celsius;
    use std::fs::read_to_string;
    use dim::Abs;
    #[test]
    fn test2() {

        use dim::si::{M,M2,K,W};
        let km = 1000.0 * M;
        let a0 = 1.0 * W/M2;
        let k  = 1.0 * W / M / K;
        let z  = 1.0 * km;
        println!("k: {:?} {:?} a0: {:?}", k, z,  a0);
    }
    #[test]
    fn layer_test() -> Result<(), crate::GeothermError> {
        // Test using South Eastern Australia from Faul and Jackson, 2005, Table 2
        use crate::Celsius;
        use crate::ContinentalLayer as Layer;
        use dim::si::{K,M,W,M3,M2};
        let km = 1000.0 * M;

        let mut stack = crate::ContinentalGeotherm::new(C!(273.15), 79e-3 * W/M2);
        let kcrust  = 2.5 * W/(M*K);
        let kmantle = 3.0 * W/(M*K);
        stack.add_layer(11.2  * km, kcrust,  3.00e-6 * W/M3);
        stack.add_layer(27.8  * km, kcrust,  1.00e-6 * W/M3);
        stack.add_layer(211.0 * km, kmantle, 0.03e-6 * W/M3);

        let (z,t) = stack.geotherm(5.0 * km).unwrap();
        let ma = crate::Adiabat::from_tp(C!(1300.0));

        for (&zi,&ti) in z.iter().zip(t.iter()) {
            let tp = ma.t(zi)?;
            let t1 = if tp < ti { tp } else { ti };
            println!("{:?} {:?} -- {:?}", zi, Celsius::from(t1), Celsius::from(ti));
        }
        Ok(())
    }
    #[test]
    fn dalton_and_faul_2010_fig8a() {
        use std::fs::write;
        use dim::si::{K,W,M,M3,M2};
        let km = 1000.0 * M;

        let mut gt = crate::ContinentalGeotherm::new(C!(0.0),
                                                     45e-3 * W/M2);
        gt.add_layer( 10. * km, 2.5 * W/(M*K), 1.480e-6 * W/M3);
        gt.add_layer( 30. * km, 2.5 * W/(M*K), 0.250e-6 * W/M3);
        gt.add_layer( 80. * km, 3.0 * W/(M*K), 0.010e-6 * W/M3);
        gt.add_layer(280. * km, 3.0 * W/(M*K), 0.084e-6 * W/M3);

        let (z0,t0) = read_tz("test/Q45.csv");
        //let mut out = String::new();
        for (z,t) in z0.iter().zip(t0.iter()) {
            let d = *z * km;
            let t2 = Celsius::from( gt.t(d).unwrap() );
            assert!((t - *t2).abs() < 45.0, "{} != {} depth: {}", t, t2, z);
            //out += &format!("{} {} {}\n", z, t, *t2);
        }
        //write("ad.txt", out).unwrap();
    }
    #[test]
    fn test_kT() {
        let kt = [(   0.0, 4.066192586184664),
                  ( 200.0, 3.3791637189091492),
                  ( 400.0, 2.950056484329319),
                  ( 600.0, 2.711763191041848),
                  ( 800.0, 2.6363825099351446),
                  (1000.0, 2.714902704854108),
                  (1200.0, 2.94833225380584),
                  (1300.0, 3.125087175451463),
                  (1400.0, 3.343398873347626),
                  (1600.0, 3.910288729963982),
        ];
        let z = 0.0 * dim::si::M;
        for (t,k) in kt.iter() {
            let t = C!(*t);
            assert_eq!(crate::kappa_mantle_jm99(t.into(), z), *k * dim::si::WPMK);
        }

    }

    fn read_tz(file: &str) -> (Vec<f64>, Vec<f64>) {
        let txt = read_to_string(file).unwrap();
        let mut t0 = vec![];
        let mut z0 = vec![];
        for line in txt.lines() {
            let mut v = line.split(",");
            t0.push( v.next().unwrap().trim().parse::<f64>().unwrap() );
            z0.push( v.next().unwrap().trim().parse::<f64>().unwrap() );
        }
        (z0,t0)
    }

    #[test]
    fn faul_2005_fig5a() -> Result<(), crate::GeothermError> {
        use crate::Celsius;
        use crate::OceanicGeotherm as OG;
        let km = 1000.0 * dim::si::M;
        let tp = C!(1350.0);
        let ages : Vec<_> = [2.,12.,35., 80., 110.].iter()
            .map(|x : &f64| x * 1e6 * crate::YEAR).collect();
        let abat = crate::Adiabat::from_tp(tp);
        let ogt : Vec<_> = ages.iter()
            .map(|x| OG::from_tp_and_age(tp, x.clone())).collect();
        for z in (0 ..= 400_i32).step_by(25) {
            let d = z as f64 * km;
            let t = abat.t(d)?;
            let ts : Result<Vec<_>,_> = ogt.iter().map(|g| g.t(d) ).collect();
            let ts : Vec<_> = ts?.iter().map(|x| Celsius::from(*x)).collect();
            println!("z,t: {:.4?} {:.4?} {:.4} ", d, ts, Celsius::from(t));
        }
        Ok(())
    }
    #[test]
    fn faul_2005_fig7a_ocean() {
        use crate::OceanicGeotherm as OG;
        let km = 1000.0 * dim::si::M;
        let tp = C!(1300.0);
        let age = 35.0 * crate::MA;
        let ogt = OG::from_tp_and_age(tp, age);
        let (z0,t0) = read_tz("test/Ocean.csv");

        for (z,t) in z0.iter().zip(t0.iter()) {
            let d = *z * km;
            let t2 = Celsius::from( ogt.t(d).unwrap() );
            assert!((t - *t2).abs() < 10.0, "{} != {} depth: {}", t, t2, z);
        }
    }
    #[test]
    fn faul_2005_fig7a_yilgarn() {
        use dim::si::{WPMK,WPM3,M};
        use crate::ContinentalGeotherm as CG;
        let mut cgt = CG::new(C!(25.0), 51e-3 * dim::si::WPM2);
        let km = 1000.0 * M;
        let hm = 37.0 * km;
        let huc = 6.6 * km;
        let rHm = 0.03e-6 * WPM3;
        cgt.add_layer(huc,     2.5 * WPMK, 2.200e-6 * WPM3);
        cgt.add_layer(hm-huc,  2.5 * WPMK, 0.800e-6 * WPM3);

        let n = 100;
        let thick = (400.0 * km - hm) / n as f64;
        dbg!(thick);
        for _ in 0 .. n {
            cgt.add_layer_kfun(thick, rHm, crate::kappa_mantle_jm99);
        }

        // Check Mantle Values
        let t0 = 436.0;
        let t = Celsius::from(cgt.t(hm).unwrap());
        assert!((*t-t0).abs() < 1.0, "{} != {}", *t, t0);
        let q0 = 12.16e-3 * dim::si::WPM2;
        let q = cgt.q(hm).unwrap();
        assert!((q-q0).abs() < 1e-5 * dim::si::WPM2, "{} != {}", q, q0);

        let zl = 398.0 * km;
        let q0 = 1.32999999e-3 * dim::si::WPM2;
        let q = cgt.q(zl).unwrap();
        assert!((q-q0).abs() < 1e-5 * dim::si::WPM2, "{} != {}", q, q0);

        let (z0,t0) = read_tz("test/Yilgarn.csv");

        //let mut out = String::new();
        for (z,t) in z0.iter().zip(t0.iter()) {
            let d = *z * km;
            let t2 = Celsius::from( cgt.t(d).unwrap() );
            assert!((t - *t2).abs() < 37.0, "{} != {} depth: {} diff: {}", t, t2, z, (t-*t2).abs());
            //out += &format!("{} {} {}\n", z, t, *t2);
        }
        //std::fs::write("ad.txt", out).unwrap();
    }
    #[test]
    fn faul_2005_fig7a_se_aus() {
        use dim::si::{WPM3,WPMK,M,WPM2};
        use crate::ContinentalGeotherm as CG;
        let km = 1000.0 * M;
        let hm = 39.0 * km;
        let huc = 11.2 * km;
        let rHuc = 3.0e-6 * WPM3;
        let rHlc = 1.0e-6 * WPM3;
        let rHm = 0.03e-6 * WPM3;
        let mut cgt = CG::new(C!(25.0), 79e-3 * WPM2);
        cgt.add_layer(huc,      2.5 * WPMK, rHuc);
        cgt.add_layer(hm-huc,   2.5 * WPMK, rHlc);
        cgt.add_layer(400.0 * km-hm, 3.0 * WPMK, rHm );

        // Check Mantle Values
        let t0 = 653.936;
        let t = Celsius::from(cgt.t(hm).unwrap());
        assert!((*t-t0).abs() < 1.0, "{} != {}", *t, t0);
        let q0 = 17.6e-3 * dim::si::WPM2;
        let q = cgt.q(hm).unwrap();
        assert!((q-q0).abs() < 1e-5 * dim::si::WPM2, "{} != {}", q, q0);

        let zl = 398.0 * km;
        let q0 = 6.83e-3 * dim::si::WPM2;
        let q = cgt.q(zl).unwrap();
        assert!((q-q0).abs() < 1e-5 * dim::si::WPM2, "{} != {}", q, q0);

        let (z0,t0) = read_tz("test/SEAus.csv");

        //let mut out = String::new();
        for (z,t) in z0.iter().zip(t0.iter()) {
            let d = *z * km;
            let t2 = Celsius::from( cgt.t(d).unwrap() );
            assert!((t - *t2).abs() < 25.0, "{} != {} depth: {} diff: {}", t, t2, z, (t-*t2).abs());
            //out += &format!("{} {} {}\n", z, t, *t2);
        }
        //std::fs::write("ad.txt", out).unwrap();
    }
        #[test]
    fn faul_2005_fig7a_n_aus() {
        use dim::si::{WPMK,WPM3,WPM2, M};
        use crate::ContinentalGeotherm as CG;
        let km = 1000.0 * M;
        let qs = 83e-3 * WPM2;
        let hm = 42.0 * km;
        let huc = 9.8 * km;
        let rHuc = 4.0e-6 * WPM3;
        let rHlc = 1.0e-6 * WPM3;
        let rHm  = 0.03e-6 * WPM3;
        let mut cgt = CG::new(C!(25.0), qs);
        cgt.add_layer(huc,      2.5 * WPMK, rHuc);
        cgt.add_layer(hm-huc,   2.5 * WPMK, rHlc);

        let n = 120;
        let thick = (500.0 * km - hm) / n as f64;
        dbg!(thick);
        for _ in 0 .. n {
            cgt.add_layer_kfun(thick, rHm, crate::kappa_mantle_jm99);
        }

        assert!( (*Celsius::from(cgt.t(huc).unwrap())- 273.528).abs() < 0.2);
        assert!( (*Celsius::from(cgt.t(hm).unwrap())- 630.304).abs() < 0.2);

        let q = cgt.q(hm).unwrap();
        let q0 = 11.6e-3 * dim::si::WPM2;
        assert!( (q - q0).abs() < 1e-3 * dim::si::WPM2, "{} != {} diff {}", q, q0, (q-q0).abs());

        let hl = 398.0 * km;
        let q = cgt.q(hl).unwrap();
        let q0 = 0.9e-3 * dim::si::WPM2;
        assert!( (q - q0).abs() < 1e-3 * dim::si::WPM2 , "{} != {} diff {}", q, q0, (q-q0).abs());

        let (z0,t0) = read_tz("test/NAus.csv");

        //let mut out = String::new();
        for (z,t) in z0.iter().zip(t0.iter()) {
            let d = *z * km;
            let t2 = Celsius::from( cgt.t(d).unwrap() );
            assert!((t - *t2).abs() < 40.0, "{} != {} depth: {} diff {}", t, t2, z, (t-*t2).abs());
            //out += &format!("{} {} {}\n", z, t, *t2);
        }
        //std::fs::write("ad.txt", out).unwrap();
    }

    #[test]
    fn faul_2005_fig4a_adiabat() -> Result<(), crate::GeothermError> {
        use std::fs::read_to_string;
        use crate::Celsius;
        use dim::si::M;
        let km = 1000.0 * M;
        let tp = C!(1300.0);
        let abat = crate::Adiabat::from_tp(tp);
        let (z0,t0) = read_tz("test/Adiabat.csv");
        for (z,t) in z0.iter().zip(t0.iter()) {
            let d = *z * km;
            let ta = Celsius::from( abat.t(d)? );
            assert!((t-*ta).abs() < 7.0, "{} != {} at {}", t, ta, z);
        }
        Ok(())
    }
    #[test]
    fn parse() {
        use crate::{KelvinB,MeterB};

        let t : KelvinB = "0 K".parse().unwrap();
        assert_eq!(*t, 0.0 * dim::si::K);
        let t : KelvinB = "0 C".parse().unwrap();
        assert_eq!(*t, 273.15 * dim::si::K);
        let t : MeterB = "0 m".parse().unwrap();
        assert_eq!(*t, 0.0 * dim::si::M);
        let t : MeterB = "1 m".parse().unwrap();
        assert_eq!(*t, 1.0 * dim::si::M);
        let t : MeterB = "1 km".parse().unwrap();
        assert_eq!(*t, 1000.0 * dim::si::M);
    }
    #[test]
    fn geotherm_either() {
        use crate::GeothermBuilder;
        let txt = r#"
type = "Oceanic"
tp  = "1300 C"
age = "100.0 Ma"
"#;
        let _b : GeothermBuilder = toml::from_str(&txt).unwrap();
    }
    
    #[test]
    fn geotherm_ocean() {
        use toml;
        use crate::OceanicGeothermBuilder;
        let txt = r#"
tp  = "1300 C"
age = "100.0 Ma"
"#;
        let b : OceanicGeothermBuilder = toml::from_str(&txt).unwrap();
        let _g = b.build();
    }
    #[test]
    fn geotherm_toml() {
        use toml;
        use crate::ContinentalGeothermBuilder;
        let txt = r#"
heat_flow = "51.0e-3 W/m^2"
ts       =  "25.0 C"

[[layers]]
top    = "0.0 km"
bottom = "10.0 km"
rh     = "2.2e-6 W/m^3"
k      = "2.5 W/mK"

[[layers]]
top    = "10.0 km"
bottom = "30.0 km"
rh     = "0.8e-6 W/m^3"
k      = "2.5 W/mK"

[[layers]]
top    = "30.0 km"
bottom = "500.0 km"
rh     = "0.03e-6 W/m^3"
k      = "Jaupart_Mareschal_1999_Eq3"
"#;
        let b : ContinentalGeothermBuilder = toml::from_str(&txt).unwrap();
        let _g = b.build();
    }
}
#[derive(Debug,Deserialize)]
#[serde(untagged)]
enum KTypeB {
    Constant(ConductivityB),
    Func(String),
}

#[derive(Debug, Deserialize)]
struct GeothermLayerBuilder {
    top: MeterB,
    bottom: MeterB,
    rh: HeatGenerationB,
    k: KTypeB,
}
/// Constructor for a `ContinentalGeotherm`
///
/// This is the preferred method of creating a continental geotherm
///   as the model description is stored in a file and documented
///   with associated units
///
/// ```rust
/// use toml;
/// use dimensioned::si::{M,K};
/// use geotherm::Geotherm;
/// let g = Geotherm::from_file("test/dalton_q45.toml").unwrap();
/// assert_eq!(g.t(0.0 * M).unwrap(),   (25.0 + 273.15) * K);
/// assert_eq!(g.t(100e3 * M).unwrap(), 1213.9499999999998 * K);
/// ```
#[derive(Debug, Deserialize)]
struct ContinentalGeothermBuilder {
    ts: KelvinB,
    heat_flow: HeatFluxB,
    layers: Vec<GeothermLayerBuilder>,
}

/// Constructor for a `OceanicGeoterm`
///
/// This is the preferred method of creating a oceanic geotherm
///   as the model description is stored in a file and documented
///   with associated units
///
#[derive(Debug, Deserialize)]
struct OceanicGeothermBuilder {
    tp: KelvinB,
    age: SecondB,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
enum GeothermBuilder {
    Continental(ContinentalGeothermBuilder),
    Oceanic(OceanicGeothermBuilder),
}

/// Geotherm Error Wrapper
#[derive(Debug)]
pub enum GeothermError {
    /// Error reading file
    Io(std::io::Error),
    /// Error parsing geotherm toml file
    Parse(toml::de::Error),
    /// Depth within geotherm is invalud
    InvalidDepth,
}
impl From<std::io::Error> for GeothermError {
    fn from(err: std::io::Error) -> Self {
        GeothermError::Io(err)
    }
}
impl From<toml::de::Error> for GeothermError {
    fn from(err: toml::de::Error) -> Self {
        GeothermError::Parse(err)
    }
}

/// Oceanic and Continental Geotherms 
///
/// ## Example
/// ```rust
/// use dimensioned::si::{M,K};
///
/// let g = geotherm::Geotherm::from_file("test/oceanic_100Ma.toml").unwrap();
/// assert_eq!(g.t(0.0 * M).unwrap(), (0.0 + 273.15) * K);
/// assert_eq!(g.t(100e3 * M).unwrap(), 1328.8354705757051 * K);
/// ```
#[derive(Debug)]
pub enum Geotherm {
    Continental(ContinentalGeotherm),
    Oceanic(OceanicGeotherm),
}

impl Geotherm {
    pub fn from_file<P: AsRef<Path>>(file: P) -> Result<Geotherm, GeothermError> {
        let b = GeothermBuilder::from_file(file)?;
        Ok(b.build())
    }
    pub fn t(&self, z: Meter<f64>) -> Result<Kelvin<f64>, GeothermError> {
        match self {
            Geotherm::Oceanic(g)     => g.t(z),
            Geotherm::Continental(g) => g.t(z)
        }
    }
}


use std::path::Path;
impl GeothermBuilder {
    /// Read a Geotherm from a toml File
    pub fn from_file<P: AsRef<Path>>(file: P) -> Result<GeothermBuilder,GeothermError> {
        let txt = std::fs::read_to_string(file)?;
        let b = toml::from_str(&txt)?;
        Ok(b)
    }
    pub fn build(self) -> Geotherm {
        match self {
            GeothermBuilder::Oceanic(g) => Geotherm::Oceanic(g.build()),
            GeothermBuilder::Continental(g) => Geotherm::Continental(g.build()),
        }
    }
}

impl OceanicGeothermBuilder {
    pub fn build(self) -> OceanicGeotherm {
        OceanicGeotherm::from_tp_and_age(*self.tp, *self.age)
    }
}

impl ContinentalGeothermBuilder {
    pub fn build(self) -> ContinentalGeotherm {
        use dim::si::{WPM2, WPMK, M, WPM3};
        let km = 1000.0 * M;
        //let mut cgt = ContinentalGeotherm::new(C!(self.ts.0), self.heat_flow * WPM2);
        let mut cgt = ContinentalGeotherm::new(*self.ts, *self.heat_flow);
        for layer in &self.layers {
            let thick = *layer.bottom - *layer.top;
            match &layer.k {
                KTypeB::Constant(v) => cgt.add_layer(thick, **v, *layer.rh),
                KTypeB::Func(name) => {
                    let fun = match name.as_str() {
                        "Jaupart_Mareschal_1999_Eq3" => kappa_mantle_jm99,
                        "Schartz_Simmons_1972" => kappa_mantle_ss72,
                        "Chapman_1986_Upper_Crust" => kappa_upper_crust_c86,
                        "Chapman_1986_Lower_Crust" => kappa_lower_crust_c86,
                        _ => panic!("Unknown Conductivity Function {}", name),
                    };
                    let min_thick = 3.0 * km; // kilometers
                    let n = if min_thick < thick {
                        1
                    } else {
                        1 + (thick / min_thick).floor() as usize
                    };
                    let thick = thick / n as f64;
                    for _ in 0 .. n as usize {
                        cgt.add_layer_kfun(thick, *layer.rh, fun);
                    }
                }
            };
        }
        cgt
    }
}
use serde::de::{self, Visitor};
use serde::de::Deserializer;
use std::num::ParseFloatError;


macro_rules! unit_serde {
    ($name:ident, $visitor:ident, $base:ty, $( $unit: expr, $conv: block ),* ) => {
        // Create newtype for dimensioned value
        #[derive(Debug)]
        struct $name ( $base );
        // Crate visitor type for newtype
        struct $visitor;
        // Deref for type to get to unit easier
        impl std::ops::Deref for $name {
            type Target = $base;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl<'de> Visitor<'de> for $visitor {
            type Value = $name;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("an floating point value with units")
            }
            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E> where E: de::Error {
                match value.parse() {
                    Ok(v) => Ok(v),
                    Err(UnitParseError::Float(x)) => Err(E::custom(x)),
                    Err(UnitParseError::NoValue) => Err(E::custom("no value found")),
                    Err(UnitParseError::NoUnits) => Err(E::custom("no units found")),
                    Err(UnitParseError::UnknownUnit) => Err(E::custom("unknown unit")),
                }
            }
        }

        impl<'de> Deserialize<'de> for $name {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where D: Deserializer<'de>
            {
                deserializer.deserialize_str( $visitor )
            }
        }
        impl std::str::FromStr for $name {
            type Err = UnitParseError;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let items : Vec<_> = s.trim().splitn(2, " ").collect();
                if items.len() == 0 {
                    return Err(UnitParseError::NoValue);
                }
                if items.len() == 1 {
                    return Err(UnitParseError::NoUnits);
                }
                let v : f64 = items[0].trim().parse()?;
                let out = match items[1].trim() {
                    $(
                        $unit  => { Ok( $name ( $conv(v) ) ) }
                    ),*
                        _ => return Err(UnitParseError::UnknownUnit),
                };
                out
            }
        }
    }
}

#[derive(Debug,PartialEq)]
enum UnitParseError {
    Float(ParseFloatError),
    NoValue,
    NoUnits,
    UnknownUnit,
}

impl From<ParseFloatError> for UnitParseError {
    fn from(err: ParseFloatError) -> Self {
        UnitParseError::Float(err)
    }
}


unit_serde!(MeterB, MeterBVisitor, Meter<f64>,
            "m", {|v| v * dim::si::M },
            "km", {|v| 1e3 * v * dim::si::M });

unit_serde!(KelvinB, KelvinBVisitor, Kelvin<f64>,
            "K", {|v| v * dim::si::K },
            "C", {|v| (v + 273.15) * dim::si::K });

unit_serde!(ConductivityB, ConductivityBVisitor, Conductivity<f64>,
            "W/mK", {|v| v * dim::si::WPMK });

unit_serde!(HeatGenerationB, HeatGenerationBVisitor, HeatGeneration<f64>,
            "W/m^3", {|v| v * dim::si::WPM3 });
unit_serde!(HeatFluxB, HeatFluxBVisitor, HeatFlux<f64>,
            "W/m^2", {|v| v * dim::si::WPM2 });
unit_serde!(SecondB, SecondBVisitor, Second<f64>,
            "s", {|v| v * dim::si::S },
            "yr", {|v| v * crate::YEAR },
            "Ma", {|v| v * crate::MA },
            "Ga", {|v| v * crate::GA }
);

