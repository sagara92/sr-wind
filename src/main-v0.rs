const GAMMA_LAW_INDEX: f64 = 4.0 / 3.0;
const SPEED_OF_LIGHT: f64 = 2.99e10;
const EPSILON: f64 = 0.05;

#[derive(Clone, Copy, Debug)]
struct Primitive {
    r: f64,
    u: f64,
    d: f64,
    h: f64,
}

#[derive(Clone, Debug)]
struct Conserved {
    f: f64,
    l: f64,
    k: f64,
    m: f64,
}

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

fn fun(x: f64, conserved: &Conserved) -> f64 {
    let Conserved { f: _, l, k: _, m: _ } = conserved;
    let c = SPEED_OF_LIGHT;
    let gm = GAMMA_LAW_INDEX;
    let e = EPSILON;
    let j = c * c + x;
    let h = c * c + x * gm;
    let u = ((l / h) * (l / h) - 1.0).sqrt();
    let k_prime = h * u + x * (gm - 1.0) / u;
    let l_prime = (1.0 + e) * (k_prime * k_prime + c.powi(4)).sqrt();

    l_prime * l_prime - j * h - k_prime * (l_prime * l_prime - h * h).sqrt()
}

fn secant_solver(primitive: Primitive) -> Primitive {
    let Primitive { r, u: _, d: _, h: _ } = primitive;
    let gm = GAMMA_LAW_INDEX;
    let c = SPEED_OF_LIGHT;
    // let _p = (h - c * c) * d * (gm - 1.0) / gm;
    let conserved = primitive_to_conserved(primitive);
    let Conserved { f, l, k: _, m:_ } = conserved;
    let n: i32 = 100;

    // root finding
    let mut x1: f64 = 0.01 * c * c;
    let mut x2: f64 = 2.0 * x1;
    let f_x1: f64 = fun(x1, &conserved);
    if f_x1 > 0.0 {
        while fun(x2, &conserved) > 0.0 {
            x1 = x2;
            x2 += 2.0 * x1;
        }
    } else if f_x1 < 0.0 {
        while fun(x2, &conserved) < 0.0 {
            x1 = x2;
            x2 += 2.0 * x1;
        }
    } else {
        println!("Exact root found at: {}", x1)
    }
    // println!("e1={}c^2 e2={}c^2", x1/c/c, x2/c/c);
    //root is between x1 and x2 if not exactly x1

    let mut a = x1;
    let mut b = x2;
    for _x in 0..n {
        let f1: f64 = fun(a, &conserved);
        let f2: f64 = fun(b, &conserved);
        let e = b - f2 * (b - a) / (f2 - f1);
        let f_new: f64 = fun(e, &conserved);
        if f1 * f_new < 0.0 {
            b = e
        } else if f2 * f_new < 0.0 {
            a = e
        } else {
            break;
        }
        // println!("fe1 = {}, fe2= {} and fe = {}",f1/c/c, f2/c/c, f_new/c/c);
    }

    let h_new = c * c + a * gm;
    let u_new = (l * l / f / f / h_new / h_new -1.0).sqrt();
    let r_new = r;
    let d_new = f / r_new /r_new / u_new;

    let primitive_new = Primitive{r: (r_new) ,u: (u_new) ,d: (d_new) ,h: (h_new)};

    primitive_new
}

fn primitive_to_conserved(primitive: Primitive) -> Conserved {
    let Primitive {r, u, d, h} = primitive;
    let c = SPEED_OF_LIGHT;
    let gm = GAMMA_LAW_INDEX;
    let p = (h - c * c) * d * (gm - 1.0) / gm;
    let f = d * r * r * u;
    let l = d * h * u * (1.0 + u * u).sqrt();
    let k = d * h * u * u + p;
    let m = 0.0;

    Conserved{f:(f), l:(l), k:(k), m:(m)}
}
fn main() {

    let mut r: f64 = 1e8; // cm
    let mut u: f64 = 10.0; // gamma-beta (dimensionless)
    let c = SPEED_OF_LIGHT;
    let h = 10.0 * c * c;
    let f = 1e20; // g / s / Sr
    let l = f * h * (1.0 + u * u).sqrt();

    let rmax: f64 = 1e10; // cm

    while r < rmax {
        let d = f / r / r / u;
        let h = l / f / (1.0 + u * u).sqrt();
        let prim = Primitive {r, u, d, h};
        let dr = 1e-3 * r;

        r += dr;
        u += du_dr(prim) * dr;
        let cnsrv: Conserved = primitive_to_conserved(prim);
        println!("Preshock{:?}", prim);
        println!("Preshock{:?}", cnsrv);
    }
    while (r >= rmax) & (r < 1.002 * rmax) {
        let d = f / r / r / u;
        let h = l / f / (1.0 + u * u).sqrt();
        let prim = Primitive {r, u, d, h};
        let dr = 1e-3 * r;
        r +=dr;
        let new_prim: Primitive = secant_solver(prim);
        let new_cnsrv: Conserved = primitive_to_conserved(new_prim);
        println!("Postshock{:?}", new_prim);
        println!("Postshock{:?}", new_cnsrv);
    }
}
