//! This module solves for a time-steady hydrodynamic wind profile, including
//! a free-expansion zone initiated at a supersonic inlet, a stationary shock
//! wave (standing reverse shock) at a specified radius, and subsequent
//! subsonic decelerating flow to an outer boundary radius.
//!
//! Authors: Sagar Adhikari, Jonathan Zrake
//!
//! TODO:
//!
//! [ ] Continue integrating the flow past the shock
//!
//! [ ] Parameterize the problem parameters: inlet state, shock radius, and
//!     outer boundary with a public data structure.
//!
//! [ ] Modify the main function to receive the problem parameter struct and
//!     return a tabulated wind solution.
//use std::fs;
use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::env;

/// The adiabatic index
const GAMMA_LAW_INDEX: f64 = 4.0 / 3.0;

/// Speed of light in cm / s
const SPEED_OF_LIGHT: f64 = 2.99e10;


/// Public data structure for problem parameters (Data given in the yaml file)
#[derive(Clone, Copy, Debug)]
pub struct Parameters {
    /// Inner boundary (cm)
    pub r_in: f64,
    /// Outer boundary (cm)
    pub r_out: f64,
    /// Shock Radius (cm)
    pub r_shock: f64,
    /// Mass loss rate (g / s / Sr)
    pub mdot: f64,
    /// Wind gamma-beta (dimensionless)
    pub u: f64,
}

/// Holds the primitive hydrodynamic variables for a spherically symmetric
/// relativistic wind.
#[derive(Clone, Copy, Debug)]
struct Primitive {
    /// Radius (cm)
    r: f64,
    /// Radial four-velocity (gamma-beta; dimensionless)
    u: f64,
    /// Mass density (g / cm^3)
    d: f64,
    /// Specific enthalpy (cm^2 / s^2)
    h: f64,
}

impl fmt::Display for Primitive {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} {} {} {}", self.r, self.u, self.d, self.h)
    }
}

impl Primitive {

    /// Create hydrodynamic variables given the radius r (cm), gamma-beta u,
    /// mass flux (over c; g) f = rho u, and specific energy flux (over c;
    /// cm^2 / s^2), l.
    fn from_rufl(r: f64, u: f64, f: f64, l: f64) -> Self {
        let d = f / u;
        let h = l / (1.0 + u * u).sqrt();
        Self {r, u, d, h}
    }

    /// Create hydrodynamic variables given the radius r (cm), gamma-beta u,
    /// mass loss rate per steradian mdot (g / s / Sr), and wind luminosity
    /// per steradian edot (erg / s / Sr).
    fn from_ru_mdot_edot(r: f64, u: f64, mdot: f64, edot: f64) -> Self {
        let c = SPEED_OF_LIGHT;
        let d = mdot / r / r / u / c;
        let h = edot / mdot / (1.0 + u * u).sqrt();
        Self {r, u, d, h}
    }

    /// Return the wind luminosity per steradian (erg / s / Sr).
    fn luminosity(&self) -> f64 {
        self.r.powi(2) * self.l() * self.f() * SPEED_OF_LIGHT
    }

    /// Return the wind mass loss rate per steradian (g / s / Sr).
    fn mass_loss_rate(&self) -> f64 {
        self.r.powi(2) * self.f() * SPEED_OF_LIGHT
    }

    /// Return the mass flux (over c; g).
    fn f(&self) -> f64 {
        self.d * self.u
    }

    /// Return the specific momentum flux (cm^2 / s^2).
    fn k(&self) -> f64 {
        let c = SPEED_OF_LIGHT;
        let h = self.h;
        let u = self.u;
        let d = self.d;
        let mu = h - c * c;
        let e = mu / GAMMA_LAW_INDEX;
        let p = d * e * (GAMMA_LAW_INDEX - 1.0);
        h * u + p / d / u
    }

    /// Return the specific luminosity (over c; cm^2 / s^2).
    fn l(&self) -> f64 {
        let h = self.h;
        let u = self.u;
        h * (1.0 + u * u).sqrt()
    }
}

/// Return the spatial derivative du/dr for an energy and momentum-conserving
/// steady-state wind.
fn du_dr(primitive: Primitive) -> f64 {
    let Primitive {r, u, d: _, h} = primitive;
    let c = SPEED_OF_LIGHT;
    let gm = GAMMA_LAW_INDEX;
    let qh = (gm - 1.0) / gm;
    let hg = h;
    let mu = hg - c * c;
    let qd = (hg / c / c - 1.0) * qh;
    let beta = u / (1.0 + u * u).sqrt();
    let beta_f = (c * c / hg * qd / (1.0 - qh)).sqrt();
    -u / r * 2.0 * mu / hg / (qh - 1.0) * qh / (beta.powi(2) - beta_f.powi(2))
}

/// Return the jump condition for a wind with specific momentum flux (k) and
/// specific luminosity (l), and a guess value for the specific internal
/// energy (e). This function returns zero when e is either of the two
/// possible values of the specific internal energy for given k and l values.
fn jump_condition(e: f64, k: f64, l: f64) -> f64 {
    let gm = GAMMA_LAW_INDEX;
    let c = SPEED_OF_LIGHT;
    (l * l - (c * c + e) * (c * c + gm * e) - k * (l * l - (c * c + gm * e).powi(2)).sqrt()) / c.powi(4)
}

/// Solves for the primitive hydrodynamic state on the downstream side of a
/// standing shock wave. The input primitive state is assumed to be
/// supersonic, and the output state is subsonic.
fn solve_jump_condition(primitive: Primitive) -> Primitive {
    let gm = GAMMA_LAW_INDEX;
    let c = SPEED_OF_LIGHT;
    let f = primitive.f();
    let k = primitive.k();
    let l = primitive.l();

    let mut e1 = (l - c * c) / gm;
    let mut e0 = e1 * 0.99999;
    let mut f0 = jump_condition(e0, k, l);
    let mut f1 = jump_condition(e1, k, l);

    while f1.abs() > 1e-10 {
        let e2 = (e0 * f1 - e1 * f0) / (f1 - f0);
        let f2 = jump_condition(e2, k, l);

        e0 = e1;
        e1 = e2;
        f0 = f1;
        f1 = f2;
    }
    let e = e1;
    let h = c * c + gm * e;
    let g = l / h;
    let u = (g * g - 1.0).sqrt();
    Primitive::from_rufl(primitive.r, u, f, l)
}

/// Return a fiducial wind inflow condition. The initial gamma-beta is 10 and
/// the terminal gamma-beta is about 100.

//fn wind_inlet() -> Primitive {
//    let mdot = Parameters.mdot;
//    let r = Parameters.r_in;
//    let u = Parameters.u;
//    let c = SPEED_OF_LIGHT;
//    let d = mdot / r / r / u / c;
//    let h = 10.0 * c * c;
//    Primitive{r, u, d, h}
//}
fn wind_inlet() -> Primitive {
    let args: Vec<String> = env::args().collect();
    let mdot: f64 = args[1].parse().unwrap();
    let r: f64 = args[2].parse().unwrap();
    let u: f64 = args[3].parse().unwrap();
    // let mdot = 1e20; // g / s / Sr
    // let r = 1e8; // cm
    // let u = 10.0; // gamma-beta (dimensionless)
    let c = SPEED_OF_LIGHT;
    let d = mdot / r / r / u / c;
    let h = 10.0 * c * c;
    Primitive{r, u, d, h}
}

fn main() {
    let inlet_prim = wind_inlet();
    let mut r = inlet_prim.r;
    let mut u = inlet_prim.u;
    let edot = inlet_prim.luminosity();
    let mdot = inlet_prim.mass_loss_rate();
    let rmax = 1e10; // cm

    let path = "/home/sagar/Documents/steady-state/data/out.dat";
    let f = File::create(path).expect("unable to create file");
    let mut f = BufWriter::new(f);

    while r < rmax {
        let dr = 1e-4 * r;
        let p = Primitive::from_ru_mdot_edot(r, u, mdot, edot);
        r += dr;
        u += du_dr(p) * dr;
        let data = p.to_string();
        writeln!(f, "{}", data).expect("unable to write file");
    }
    let pre_shock_prim = Primitive::from_ru_mdot_edot(r, u, mdot, edot);
    let pos_shock_prim = solve_jump_condition(pre_shock_prim);
    u = pos_shock_prim.u;
    while (r > rmax) & (r < 100.0 * rmax) {
        let p = Primitive::from_ru_mdot_edot(r, u, mdot, edot);
        let dr = 1e-4 * r;
        r += dr;
        u += du_dr(p) * dr;
        let data = p.to_string();
        writeln!(f, "{}", data).expect("unable to write file");
    }
    let final_prim = Primitive::from_ru_mdot_edot(r, u, mdot, edot);

    println!("upstream:      {:?}", pre_shock_prim);
    println!("downstream:    {:?}", pos_shock_prim);
    println!("end of stream: {:?}", final_prim);
}
