use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

#[cfg(test)]
use roots::{find_root_brent, SimpleConvergency};

// Substantial portions of this code are based on Brown and Caldwell's TidyWater package.
// Copyright (c) 2025 Brown and Caldwell
// MIT License

/// Carries validation and convergence failures back to R as readable messages.
type SolverResult<T> = std::result::Result<T, String>;

/// Bounds stack storage because the public schema supports at most three steps.
const MAX_SPECIES: usize = 4;

/// Bounds dissociation-constant storage to one fewer entry than species storage.
const MAX_CONSTANTS: usize = MAX_SPECIES - 1;

/// Lower endpoint of the existing hydrogen-concentration search interval.
const LOG_H_MIN: f64 = -14.0 * std::f64::consts::LN_10;

/// Upper endpoint of the existing hydrogen-concentration search interval.
const LOG_H_MAX: f64 = 0.0;

/// Gives the root finder approximately 1e-12 pH internal precision.
const LOG_H_TOLERANCE: f64 = std::f64::consts::LN_10 * 1e-12;

/// Caps safeguarded iterations well above the bisection worst case.
const MAX_ROOT_ITERATIONS: usize = 100;

/// Stores one dissociation step after its R definition has been validated.
#[derive(Debug, Clone, Copy)]
struct Ion {
    k: f64,
    delta_h: f64,
}

/// Distinguishes systems whose charged species lie on either side of neutral.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChargeDirection {
    Negative,
    Positive,
}

impl ChargeDirection {
    /// Converts the direction into the sign used when summing charge balance.
    fn sign(self) -> i8 {
        match self {
            Self::Negative => -1,
            Self::Positive => 1,
        }
    }
}

/// Caches the immutable chemistry supplied when a solver is constructed.
#[derive(Debug, Clone)]
struct DependentDefinition {
    ions: Vec<Ion>,
    charge_direction: ChargeDirection,
}

/// Combines cached chemistry with concentrations supplied for one solve call.
#[derive(Debug, Clone)]
struct DependentCompound<'a> {
    ions: &'a [Ion],
    charge_direction: ChargeDirection,
    total: f64,
    initial: [f64; MAX_CONSTANTS],
}

/// Holds one fixed-charge dose in the form needed by the charge balance.
#[derive(Debug, Clone)]
struct IndependentCompound {
    charge: i8,
    dose: f64,
}

/// Owns the parsed catalogs behind the external pointer returned to R.
#[derive(Debug)]
struct SolverConfig {
    dependent_definitions: HashMap<String, DependentDefinition>,
    independent_definitions: HashMap<String, i8>,
}

/// Keeps the small set of species fractions on the stack during root finding.
#[cfg(test)]
#[derive(Debug, Clone, Copy)]
struct AlphaFractions {
    values: [f64; MAX_SPECIES],
    len: usize,
}

#[cfg(test)]
impl AlphaFractions {
    /// Exposes only the entries populated for the compound's actual valence.
    fn as_slice(&self) -> &[f64] {
        &self.values[..self.len]
    }
}

/// Stores cumulative corrected log constants for repeated root evaluations.
#[derive(Debug, Clone, Copy)]
struct PreparedCompound {
    cumulative_log_k: [f64; MAX_SPECIES],
    len: usize,
    charge_direction: ChargeDirection,
    total: f64,
}

/// Carries a compound's charge and derivative contribution from one softmax.
#[derive(Debug, Clone, Copy)]
struct ChargeMoments {
    charge: f64,
    derivative: f64,
}

/// Carries charge balance and its derivative with respect to log hydrogen.
#[derive(Debug, Clone, Copy)]
struct BalancePoint {
    balance: f64,
    derivative: f64,
}

/// Holds values shared by every numerical evaluation in one solve.
struct ChargeBalance<'a> {
    compounds: &'a [PreparedCompound],
    kw: f64,
    gamma_h: f64,
    dose_charge: f64,
    starting_charge: f64,
}

impl ChargeBalance<'_> {
    /// Evaluates charge balance and its positive analytic derivative together.
    fn evaluate(&self, log_h: f64) -> BalancePoint {
        let h = log_h.exp();
        let oh = self.kw / (h * self.gamma_h * self.gamma_h);
        let mut dependent_charge = 0.0;
        let mut dependent_derivative = 0.0;
        for compound in self.compounds {
            let moments = charge_moments(compound, log_h);
            dependent_charge += moments.charge;
            dependent_derivative += moments.derivative;
        }
        BalancePoint {
            balance: h - oh + dependent_charge + self.dose_charge - self.starting_charge,
            derivative: h + oh + dependent_derivative,
        }
    }
}

/// Contains a root and test-only evaluation instrumentation.
struct RootSolution {
    log_h: f64,
    #[cfg(test)]
    evaluations: usize,
}

/// Computes a Davies activity coefficient so non-ideal solutions can be solved.
fn calculate_activity(charge: i8, ionic_strength: Option<f64>, temp: f64) -> f64 {
    match ionic_strength {
        Some(value) => {
            let temp_abs = temp + 273.15;
            let de = 78.54
                * (1.0 - 0.004579 * (temp_abs - 298.0)
                    + 11.9e-6 * (temp_abs - 298.0).powi(2)
                    + 28e-9 * (temp_abs - 298.0).powi(3));
            let a = 1.29e6 * (2_f64.sqrt() / (de * temp_abs).powf(1.5));
            let sqrt_i = value.sqrt();
            10_f64.powf(-a * f64::from(charge).powi(2) * (sqrt_i / (1.0 + sqrt_i) - 0.3 * value))
        }
        None => 1.0,
    }
}

/// Precomputes every supported activity coefficient once per solve.
fn activity_coefficients(ionic_strength: Option<f64>, temp: f64) -> [f64; MAX_SPECIES] {
    let gamma_one = calculate_activity(1, ionic_strength, temp);
    [1.0, gamma_one, gamma_one.powi(4), gamma_one.powi(9)]
}

/// Moves a 25 °C equilibrium constant to the requested temperature.
fn k_temp_adjust(delta_h: f64, k: f64, temp: f64) -> f64 {
    let gas_constant = 8.314;
    let temp_abs = temp + 273.15;
    ((delta_h / gas_constant * (1.0 / 298.15 - 1.0 / temp_abs)) + k.ln()).exp()
}

/// Identifies each dissociation transition so its activity correction uses the
/// correct reactant and product charges in Benjamin's distance-from-neutral order.
fn ion_charges(direction: ChargeDirection, ion_index: usize) -> (i8, i8) {
    let charge_magnitude = (ion_index + 1) as i8;
    match direction {
        ChargeDirection::Negative => (-(charge_magnitude - 1), -charge_magnitude),
        ChargeDirection::Positive => (charge_magnitude, charge_magnitude - 1),
    }
}

/// Applies temperature and activity corrections once before root iteration.
fn correct_k(
    compound: &DependentCompound<'_>,
    temp: f64,
    gammas: &[f64; MAX_SPECIES],
) -> [f64; MAX_CONSTANTS] {
    let gamma_h = gammas[1];
    let mut corrected = [0.0; MAX_CONSTANTS];
    for (ion_index, ion) in compound.ions.iter().enumerate() {
        let (reactant_charge, product_charge) = ion_charges(compound.charge_direction, ion_index);
        let temperature_corrected_k = k_temp_adjust(ion.delta_h, ion.k, temp);
        let reactant_activity = gammas[usize::from(reactant_charge.unsigned_abs())];
        let product_activity = gammas[usize::from(product_charge.unsigned_abs())];
        corrected[ion_index] =
            temperature_corrected_k * reactant_activity / (gamma_h * product_activity);
    }
    corrected
}

/// Converts corrected constants into cumulative logs used by every root step.
fn prepare_compound(
    compound: &DependentCompound<'_>,
    temp: f64,
    gammas: &[f64; MAX_SPECIES],
) -> PreparedCompound {
    let corrected = correct_k(compound, temp, gammas);
    let mut cumulative_log_k = [0.0; MAX_SPECIES];
    for index in 0..compound.ions.len() {
        cumulative_log_k[index + 1] = cumulative_log_k[index] + corrected[index].ln();
    }
    PreparedCompound {
        cumulative_log_k,
        len: compound.ions.len() + 1,
        charge_direction: compound.charge_direction,
        total: compound.total,
    }
}

/// Fuses stable alpha normalization with charge mean and variance calculation.
fn charge_moments(compound: &PreparedCompound, log_h: f64) -> ChargeMoments {
    let mut log_weights = [0.0; MAX_SPECIES];
    for (index, weight) in log_weights[..compound.len].iter_mut().enumerate().skip(1) {
        let index = index as f64;
        *weight = match compound.charge_direction {
            ChargeDirection::Negative => compound.cumulative_log_k[index as usize] - index * log_h,
            ChargeDirection::Positive => index * log_h - compound.cumulative_log_k[index as usize],
        };
    }
    let maximum = log_weights[..compound.len]
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut total_weight = 0.0;
    let mut first_moment = 0.0;
    let mut second_moment = 0.0;
    for (index, log_weight) in log_weights[..compound.len].iter().enumerate() {
        let weight = (*log_weight - maximum).exp();
        let index = index as f64;
        total_weight += weight;
        first_moment += index * weight;
        second_moment += index * index * weight;
    }
    let mean = first_moment / total_weight;
    let variance = (second_moment / total_weight - mean * mean).max(0.0);
    ChargeMoments {
        charge: f64::from(compound.charge_direction.sign()) * compound.total * mean,
        derivative: compound.total * variance,
    }
}

/// Alpha fractions indexed by distance from neutral charge.
///
/// For negative compounds each step ratio is K_i / [H+]. For positive
/// compounds it is [H+] / K_i because K_i describes dissociation from charge
/// +i toward charge +(i-1). This exists to provide all species fractions needed
/// by charge balance while log-space normalization keeps extreme cases stable.
#[cfg(test)]
fn calculate_alphas(h: f64, ks: &[f64], direction: ChargeDirection) -> AlphaFractions {
    let len = ks.len() + 1;
    debug_assert!(len <= MAX_SPECIES);
    let log_h = h.ln();
    let mut cumulative_log_weight = 0.0;
    let mut weights = [0.0; MAX_SPECIES];
    for (index, k) in ks.iter().enumerate() {
        cumulative_log_weight += match direction {
            ChargeDirection::Negative => k.ln() - log_h,
            ChargeDirection::Positive => log_h - k.ln(),
        };
        weights[index + 1] = cumulative_log_weight;
    }
    let maximum = weights[..len]
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut total_weight = 0.0;
    for weight in &mut weights[..len] {
        *weight = (*weight - maximum).exp();
        total_weight += *weight;
    }
    for weight in &mut weights[..len] {
        *weight /= total_weight;
    }
    AlphaFractions {
        values: weights,
        len,
    }
}

/// Direct-arithmetic equivalent of [`calculate_alphas`].
///
/// This version is retained as a readable statement of the equilibrium
/// equations, but is deliberately not used by the solver because cumulative
/// products can overflow or underflow for extreme constants or pH values.
#[allow(dead_code)]
#[cfg(test)]
fn calculate_alphas_without_logs(h: f64, ks: &[f64], direction: ChargeDirection) -> AlphaFractions {
    let len = ks.len() + 1;
    debug_assert!(len <= MAX_SPECIES);
    let mut cumulative_weight = 1.0;
    let mut weights = [0.0; MAX_SPECIES];
    weights[0] = cumulative_weight;
    for (index, k) in ks.iter().enumerate() {
        cumulative_weight *= match direction {
            ChargeDirection::Negative => k / h,
            ChargeDirection::Positive => h / k,
        };
        weights[index + 1] = cumulative_weight;
    }
    let total_weight: f64 = weights[..len].iter().sum();
    for weight in &mut weights[..len] {
        *weight /= total_weight;
    }
    AlphaFractions {
        values: weights,
        len,
    }
}

/// Converts species fractions into the compound's equilibrium charge contribution.
#[cfg(test)]
fn equilibrium_charge(compound: &DependentCompound<'_>, alphas: &AlphaFractions) -> f64 {
    let signed_fraction: f64 = alphas
        .as_slice()
        .iter()
        .enumerate()
        .map(|(charge_magnitude, alpha)| {
            f64::from(compound.charge_direction.sign()) * charge_magnitude as f64 * alpha
        })
        .sum();
    compound.total * signed_fraction
}

/// Reconstructs charge present before equilibration from the supplied ion states.
fn initial_charge(compound: &DependentCompound<'_>) -> f64 {
    compound.initial[..compound.ions.len()]
        .iter()
        .enumerate()
        .map(|(index, concentration)| {
            concentration * f64::from(compound.charge_direction.sign()) * (index + 1) as f64
        })
        .sum()
}

/// Finds the unique root with Newton steps guarded by a persistent bracket.
fn find_root_safeguarded(balance: &ChargeBalance<'_>) -> SolverResult<RootSolution> {
    #[cfg(test)]
    let mut evaluations = 0;
    #[allow(unused_mut)]
    let mut evaluate = |log_h| {
        #[cfg(test)]
        {
            evaluations += 1;
        }
        balance.evaluate(log_h)
    };

    let mut lower = LOG_H_MIN;
    let mut upper = LOG_H_MAX;
    let lower_point = evaluate(lower);
    let upper_point = evaluate(upper);
    if !lower_point.balance.is_finite()
        || !upper_point.balance.is_finite()
        || lower_point.balance > 0.0
        || upper_point.balance < 0.0
    {
        return Err("The pH solver failed to converge to a solution on [0, 14].".to_string());
    }
    if lower_point.balance == 0.0 {
        return Ok(RootSolution {
            log_h: lower,
            #[cfg(test)]
            evaluations,
        });
    }
    if upper_point.balance == 0.0 {
        return Ok(RootSolution {
            log_h: upper,
            #[cfg(test)]
            evaluations,
        });
    }

    let mut current = -7.0 * std::f64::consts::LN_10;
    for _ in 0..MAX_ROOT_ITERATIONS {
        let point = evaluate(current);
        if !point.balance.is_finite() || !point.derivative.is_finite() || point.derivative <= 0.0 {
            return Err("The pH solver failed to converge to a solution on [0, 14].".to_string());
        }
        if point.balance < 0.0 {
            lower = current;
        } else if point.balance > 0.0 {
            upper = current;
        } else {
            return Ok(RootSolution {
                log_h: current,
                #[cfg(test)]
                evaluations,
            });
        }

        let bracket_width = upper - lower;
        if bracket_width <= LOG_H_TOLERANCE {
            return Ok(RootSolution {
                log_h: (lower + upper) / 2.0,
                #[cfg(test)]
                evaluations,
            });
        }

        let newton = current - point.balance / point.derivative;
        let next = if newton.is_finite()
            && newton > lower
            && newton < upper
            && (newton - current).abs() <= bracket_width / 2.0
        {
            newton
        } else {
            (lower + upper) / 2.0
        };
        if next == current {
            return Ok(RootSolution {
                log_h: current,
                #[cfg(test)]
                evaluations,
            });
        }
        current = next;
    }
    Err("The pH solver failed to converge to a solution on [0, 14].".to_string())
}

/// Orders Brent endpoints so the better residual is the second point.
#[cfg(test)]
fn arrange_points(a: f64, fa: f64, b: f64, fb: f64) -> (f64, f64, f64, f64) {
    if fa.abs() > fb.abs() {
        (a, fa, b, fb)
    } else {
        (b, fb, a, fa)
    }
}

/// Holds a cached Brent root and test-only function-evaluation count.
#[cfg(test)]
struct CachedBrentSolution {
    root: f64,
    #[cfg(test)]
    evaluations: usize,
}

/// Holds a completed pH solve and test-only numerical instrumentation.
#[cfg(test)]
struct SolveOutcome {
    ph: f64,
    #[cfg(test)]
    evaluations: usize,
}

/// Reproduces the existing Brent algorithm without re-evaluating cached endpoints.
#[cfg(test)]
fn find_root_brent_cached<Func>(
    lower: f64,
    upper: f64,
    mut function: Func,
    tolerance: f64,
    max_iterations: usize,
) -> SolverResult<CachedBrentSolution>
where
    Func: FnMut(f64) -> f64,
{
    #[cfg(test)]
    let mut evaluations = 0;
    let mut evaluate = |value| {
        #[cfg(test)]
        {
            evaluations += 1;
        }
        function(value)
    };
    let lower_value = evaluate(lower);
    let upper_value = evaluate(upper);
    let (mut a, mut fa, mut b, mut fb) = arrange_points(lower, lower_value, upper, upper_value);
    if !fa.is_finite() || !fb.is_finite() || fa * fb > 0.0 {
        return Err("The pH solver failed to converge to a solution on [0, 14].".to_string());
    }
    let (mut c, mut fc, mut d) = (a, fa, a);
    let mut bisected = true;

    for _ in 0..max_iterations {
        if fa.abs() < tolerance {
            return Ok(CachedBrentSolution {
                root: a,
                #[cfg(test)]
                evaluations,
            });
        }
        if fb.abs() < tolerance {
            return Ok(CachedBrentSolution {
                root: b,
                #[cfg(test)]
                evaluations,
            });
        }
        if (a - b).abs() < tolerance {
            return Ok(CachedBrentSolution {
                root: c,
                #[cfg(test)]
                evaluations,
            });
        }

        let mut candidate = if fa != fc && fb != fc {
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };
        let outside_safe_region = (candidate - b) * (candidate - (3.0 * a + b) / 4.0) > 0.0;
        let insufficient_progress = bisected && (candidate - b).abs() >= (b - c).abs() / 2.0;
        let prior_insufficient_progress = !bisected && (candidate - b).abs() >= (c - d).abs() / 2.0;
        let stale_bisection = bisected && (b - c).abs() < tolerance;
        let stale_interpolation = !bisected && (c - d).abs() < tolerance;
        if outside_safe_region
            || insufficient_progress
            || prior_insufficient_progress
            || stale_bisection
            || stale_interpolation
        {
            candidate = (a + b) / 2.0;
            bisected = true;
        } else {
            bisected = false;
        }

        let candidate_value = evaluate(candidate);
        d = c;
        c = b;
        fc = fb;
        if fa * candidate_value < 0.0 {
            (a, fa, b, fb) = arrange_points(a, fa, candidate, candidate_value);
        } else {
            (a, fa, b, fb) = arrange_points(candidate, candidate_value, b, fb);
        }
    }
    Err("The pH solver failed to converge to a solution on [0, 14].".to_string())
}

/// Provides a cached-evaluation log-space Brent candidate for comparisons.
#[cfg(test)]
fn find_root_log_brent(balance: &ChargeBalance<'_>) -> SolverResult<RootSolution> {
    let mut evaluations = 0;
    let mut evaluate = |log_h| {
        evaluations += 1;
        balance.evaluate(log_h).balance
    };
    let fa = evaluate(LOG_H_MIN);
    let fb = evaluate(LOG_H_MAX);
    let (mut a, mut fa, mut b, mut fb) = arrange_points(LOG_H_MIN, fa, LOG_H_MAX, fb);
    if !fa.is_finite() || !fb.is_finite() || fa * fb > 0.0 {
        return Err("The pH solver failed to converge to a solution on [0, 14].".to_string());
    }
    let (mut c, mut fc, mut d) = (a, fa, a);
    let mut bisected = true;

    for _ in 0..MAX_ROOT_ITERATIONS {
        if fa == 0.0 {
            return Ok(RootSolution {
                log_h: a,
                evaluations,
            });
        }
        if fb == 0.0 || (a - b).abs() <= LOG_H_TOLERANCE {
            return Ok(RootSolution {
                log_h: b,
                evaluations,
            });
        }
        let mut candidate = if fa != fc && fb != fc {
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };
        let outside_safe_region = (candidate - b) * (candidate - (3.0 * a + b) / 4.0) > 0.0;
        let insufficient_progress = bisected && (candidate - b).abs() >= (b - c).abs() / 2.0;
        let prior_insufficient_progress = !bisected && (candidate - b).abs() >= (c - d).abs() / 2.0;
        let stale_bisection = bisected && (b - c).abs() <= LOG_H_TOLERANCE;
        let stale_interpolation = !bisected && (c - d).abs() <= LOG_H_TOLERANCE;
        if outside_safe_region
            || insufficient_progress
            || prior_insufficient_progress
            || stale_bisection
            || stale_interpolation
        {
            candidate = (a + b) / 2.0;
            bisected = true;
        } else {
            bisected = false;
        }

        let candidate_balance = evaluate(candidate);
        d = c;
        c = b;
        fc = fb;
        if fa * candidate_balance < 0.0 {
            (a, fa, b, fb) = arrange_points(a, fa, candidate, candidate_balance);
        } else {
            (a, fa, b, fb) = arrange_points(candidate, candidate_balance, b, fb);
        }
    }
    Err("The pH solver failed to converge to a solution on [0, 14].".to_string())
}

/// Prepares a charge balance shared by production and reference root solvers.
fn prepare_balance<'a>(
    temp: f64,
    ionic_strength: Option<f64>,
    dependent_compounds: &[DependentCompound<'a>],
    independent_compounds: &[IndependentCompound],
    h_i: f64,
    oh_i: f64,
) -> (Vec<PreparedCompound>, f64, f64, f64) {
    let gammas = activity_coefficients(ionic_strength, temp);
    let prepared: Vec<PreparedCompound> = dependent_compounds
        .iter()
        .map(|compound| prepare_compound(compound, temp, &gammas))
        .collect();
    let gamma_h = gammas[1];
    let dose_charge: f64 = independent_compounds
        .iter()
        .map(|compound| f64::from(compound.charge) * compound.dose)
        .sum();
    let starting_charge = h_i - oh_i + dependent_compounds.iter().map(initial_charge).sum::<f64>();
    (prepared, gamma_h, dose_charge, starting_charge)
}

/// Solves pH with safeguarded Newton iteration in log-hydrogen space.
fn solve_internal(
    temp: f64,
    ionic_strength: Option<f64>,
    kw: f64,
    dependent_compounds: &[DependentCompound<'_>],
    independent_compounds: &[IndependentCompound],
    h_i: f64,
    oh_i: f64,
) -> SolverResult<f64> {
    let (prepared, gamma_h, dose_charge, starting_charge) = prepare_balance(
        temp,
        ionic_strength,
        dependent_compounds,
        independent_compounds,
        h_i,
        oh_i,
    );
    let balance = ChargeBalance {
        compounds: &prepared,
        kw,
        gamma_h,
        dose_charge,
        starting_charge,
    };
    let root = find_root_safeguarded(&balance)?;
    Ok(-(root.log_h + gamma_h.ln()) / std::f64::consts::LN_10)
}

/// Runs concentration-space Brent while reusing all previously evaluated points.
#[cfg(test)]
fn solve_internal_cached_brent(
    temp: f64,
    ionic_strength: Option<f64>,
    kw: f64,
    dependent_compounds: &[DependentCompound<'_>],
    independent_compounds: &[IndependentCompound],
    h_i: f64,
    oh_i: f64,
) -> SolverResult<SolveOutcome> {
    let gammas = activity_coefficients(ionic_strength, temp);
    let corrected_ks: Vec<[f64; MAX_CONSTANTS]> = dependent_compounds
        .iter()
        .map(|compound| correct_k(compound, temp, &gammas))
        .collect();
    let gamma_h = gammas[1];
    let dose_charge: f64 = independent_compounds
        .iter()
        .map(|compound| f64::from(compound.charge) * compound.dose)
        .sum();
    let starting_charge = h_i - oh_i + dependent_compounds.iter().map(initial_charge).sum::<f64>();
    let charge_balance = |h: f64| {
        let oh = kw / (h * gamma_h * gamma_h);
        let dependent_charge: f64 = dependent_compounds
            .iter()
            .zip(&corrected_ks)
            .map(|(compound, corrected)| {
                let alphas = calculate_alphas(
                    h,
                    &corrected[..compound.ions.len()],
                    compound.charge_direction,
                );
                equilibrium_charge(compound, &alphas)
            })
            .sum();
        h - oh + dependent_charge + dose_charge - starting_charge
    };
    let root = find_root_brent_cached(1e-14, 1.0, charge_balance, 1e-14, 1000)?;
    Ok(SolveOutcome {
        ph: -(root.root * gamma_h).log10(),
        #[cfg(test)]
        evaluations: root.evaluations,
    })
}

/// Retains the previous concentration-space Brent solver as a test oracle.
#[cfg(test)]
fn solve_internal_brent_reference(
    temp: f64,
    ionic_strength: Option<f64>,
    kw: f64,
    dependent_compounds: &[DependentCompound<'_>],
    independent_compounds: &[IndependentCompound],
    h_i: f64,
    oh_i: f64,
) -> SolverResult<(f64, usize)> {
    let gammas = activity_coefficients(ionic_strength, temp);
    let corrected_ks: Vec<[f64; MAX_CONSTANTS]> = dependent_compounds
        .iter()
        .map(|compound| correct_k(compound, temp, &gammas))
        .collect();
    let gamma_h = gammas[1];
    let dose_charge: f64 = independent_compounds
        .iter()
        .map(|compound| f64::from(compound.charge) * compound.dose)
        .sum();
    let starting_charge = h_i - oh_i + dependent_compounds.iter().map(initial_charge).sum::<f64>();
    let mut evaluations = 0;
    let charge_balance = |h: f64| {
        evaluations += 1;
        let oh = kw / (h * gamma_h * gamma_h);
        let dependent_charge: f64 = dependent_compounds
            .iter()
            .zip(&corrected_ks)
            .map(|(compound, corrected)| {
                let alphas = calculate_alphas(
                    h,
                    &corrected[..compound.ions.len()],
                    compound.charge_direction,
                );
                equilibrium_charge(compound, &alphas)
            })
            .sum();
        h - oh + dependent_charge + dose_charge - starting_charge
    };
    let mut convergency = SimpleConvergency {
        eps: 1e-14,
        max_iter: 1000,
    };
    let h = find_root_brent(1e-14, 1.0, charge_balance, &mut convergency)
        .map_err(|_| "The pH solver failed to converge to a solution on [0, 14].".to_string())?;
    Ok((-(h * gamma_h).log10(), evaluations))
}

/// Requires an R value to be a list while preserving its location in errors.
fn as_list(value: &Robj, path: &str) -> SolverResult<List> {
    value
        .as_list()
        .ok_or_else(|| format!("`{path}` must be a list."))
}

/// Extracts named list entries and rejects names that make lookup ambiguous.
fn named_entries(list: &List, path: &str) -> SolverResult<Vec<(String, Robj)>> {
    if !list.is_empty() && list.names().is_none() {
        return Err(format!(
            "Every entry in `{path}` must have a non-empty name."
        ));
    }
    let entries: Vec<(String, Robj)> = list
        .iter()
        .map(|(name, value)| (name.to_string(), value))
        .collect();
    if entries.iter().any(|(name, _)| name.is_empty()) {
        return Err(format!(
            "Every entry in `{path}` must have a non-empty name."
        ));
    }
    let mut seen = HashSet::new();
    if entries.iter().any(|(name, _)| !seen.insert(name.clone())) {
        return Err(format!("Names in `{path}` must be unique."));
    }
    Ok(entries)
}

/// Enforces exact nested-list schemas so misspelled fields cannot be ignored.
fn exact_fields(value: &Robj, path: &str, expected: &[&str]) -> SolverResult<List> {
    let list = as_list(value, path)?;
    let entries = named_entries(&list, path)?;
    let actual: HashSet<&str> = entries.iter().map(|(name, _)| name.as_str()).collect();
    let wanted: HashSet<&str> = expected.iter().copied().collect();
    if actual != wanted || entries.len() != expected.len() {
        return Err(format!(
            "`{path}` must contain exactly {}.",
            expected
                .iter()
                .map(|field| format!("`{field}`"))
                .collect::<Vec<_>>()
                .join(" and ")
        ));
    }
    Ok(list)
}

/// Retrieves a field already proven to exist by [`exact_fields`].
fn field(list: &List, name: &str) -> Robj {
    list.iter()
        .find_map(|(field_name, value)| (field_name == name).then_some(value))
        .expect("exact_fields guarantees required fields")
}

/// Accepts R integer or double scalars and rejects missing or infinite values.
fn scalar_number(value: &Robj, path: &str) -> SolverResult<f64> {
    let number = value
        .as_real()
        .or_else(|| value.as_integer().map(f64::from))
        .ok_or_else(|| format!("`{path}` must be one finite number."))?;
    if !number.is_finite() {
        return Err(format!("`{path}` must be one finite number."));
    }
    Ok(number)
}

/// Validates concentrations and ionic strengths that may be zero.
fn nonnegative_number(value: &Robj, path: &str) -> SolverResult<f64> {
    let number = scalar_number(value, path)?;
    if number < 0.0 {
        return Err(format!("`{path}` must be one finite number >= 0."));
    }
    Ok(number)
}

/// Validates constants such as K values that must be strictly positive.
fn positive_number(value: &Robj, path: &str) -> SolverResult<f64> {
    let number = scalar_number(value, path)?;
    if number <= 0.0 {
        return Err(format!("`{path}` must be one finite number > 0."));
    }
    Ok(number)
}

/// Represents R `NULL` as an omitted nonnegative numerical correction.
fn optional_nonnegative_number(value: &Robj, path: &str) -> SolverResult<Option<f64>> {
    if value.is_null() {
        Ok(None)
    } else {
        nonnegative_number(value, path).map(Some)
    }
}

/// Turns the public -1/1 direction convention into a type-safe Rust enum.
fn parse_charge_direction(value: &Robj, path: &str) -> SolverResult<ChargeDirection> {
    match scalar_number(value, path)? {
        -1.0 => Ok(ChargeDirection::Negative),
        1.0 => Ok(ChargeDirection::Positive),
        _ => Err(format!("`{path}` must be either -1 or 1.")),
    }
}

/// Validates a fixed ion's signed integral valence before narrowing it to `i8`.
fn parse_independent_charge(value: &Robj, path: &str) -> SolverResult<i8> {
    let charge = scalar_number(value, path)?;
    if charge == 0.0 || charge.fract() != 0.0 || charge < i8::MIN as f64 || charge > i8::MAX as f64
    {
        return Err(format!(
            "`{path}` must be a nonzero integer between -128 and 127."
        ));
    }
    Ok(charge as i8)
}

/// Parses pH-dependent catalog entries once for storage in [`SolverConfig`].
fn parse_dependent_definitions(list: &List) -> SolverResult<HashMap<String, DependentDefinition>> {
    named_entries(list, "ph_dependent")?
        .into_iter()
        .map(|(name, value)| {
            let path = format!("ph_dependent${name}");
            let definition = exact_fields(&value, &path, &["constants", "charge"])?;
            let constants_path = format!("{path}$constants");
            let constants = as_list(&field(&definition, "constants"), &constants_path)?;
            if constants.is_empty() {
                return Err(format!("`{constants_path}` must be a non-empty list."));
            }
            if constants.len() >= MAX_SPECIES {
                return Err(format!(
                    "`{constants_path}` must contain one to three dissociation constants (two to four ion states)."
                ));
            }
            let ions = constants
                .values()
                .enumerate()
                .map(|(index, value)| {
                    let ion_path = format!("{constants_path}[[{}]]", index + 1);
                    let ion = exact_fields(&value, &ion_path, &["k", "delta_h"])?;
                    Ok(Ion {
                        k: positive_number(&field(&ion, "k"), &format!("{ion_path}$k"))?,
                        delta_h: scalar_number(
                            &field(&ion, "delta_h"),
                            &format!("{ion_path}$delta_h"),
                        )?,
                    })
                })
                .collect::<SolverResult<Vec<_>>>()?;
            let charge_direction = parse_charge_direction(
                &field(&definition, "charge"),
                &format!("{path}$charge"),
            )?;
            Ok((
                name,
                DependentDefinition {
                    ions,
                    charge_direction,
                },
            ))
        })
        .collect()
}

/// Parses the flat fixed-charge catalog once for storage in [`SolverConfig`].
fn parse_independent_definitions(list: &List) -> SolverResult<HashMap<String, i8>> {
    named_entries(list, "ph_independent_charges")?
        .into_iter()
        .map(|(name, value)| {
            let path = format!("ph_independent_charges${name}");
            let charge = parse_independent_charge(&value, &path)?;
            Ok((name, charge))
        })
        .collect()
}

/// Builds both catalogs together and prevents a name from having two meanings.
fn validate_catalog(
    ph_dependent: &List,
    ph_independent_charges: &List,
) -> SolverResult<(HashMap<String, DependentDefinition>, HashMap<String, i8>)> {
    let dependent = parse_dependent_definitions(ph_dependent)?;
    let independent = parse_independent_definitions(ph_independent_charges)?;
    let mut overlap: Vec<&String> = dependent
        .keys()
        .filter(|name| independent.contains_key(*name))
        .collect();
    overlap.sort();
    if !overlap.is_empty() {
        return Err(format!(
            "Compounds cannot be both pH-dependent and pH-independent: {}.",
            overlap
                .iter()
                .map(|name| name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    Ok((dependent, independent))
}

/// Parse changing pH-dependent concentrations against the cached catalog.
///
/// The expected R value is a named list with one entry per active compound:
/// `list(co3 = list(total = 0.0025, initial = list(0.0024, 0)))`.
/// `initial` is ordered by absolute charge magnitude and must contain one value
/// for every dissociation constant in that compound's cached definition.
fn parse_dependent_compounds<'a>(
    values: &List,
    definitions: &'a HashMap<String, DependentDefinition>,
) -> SolverResult<Vec<DependentCompound<'a>>> {
    if values.is_empty() {
        return Err("`ph_dependent` must contain at least one compound.".to_string());
    }
    named_entries(values, "ph_dependent")?
        .into_iter()
        .map(|(name, value)| {
            let definition = definitions
                .get(&name)
                .ok_or_else(|| format!("Unknown compound in `ph_dependent`: {name}."))?;
            let path = format!("ph_dependent${name}");
            let runtime = exact_fields(&value, &path, &["total", "initial"])?;
            let initial_path = format!("{path}$initial");
            let initial_values = as_list(&field(&runtime, "initial"), &initial_path)?;
            if initial_values.len() != definition.ions.len() {
                return Err(format!(
                    "`{initial_path}` must be a list of length {}.",
                    definition.ions.len()
                ));
            }
            let mut initial = [0.0; MAX_CONSTANTS];
            for (index, value) in initial_values.values().enumerate() {
                initial[index] =
                    nonnegative_number(&value, &format!("{initial_path}[[{}]]", index + 1))?;
            }
            Ok(DependentCompound {
                ions: &definition.ions,
                charge_direction: definition.charge_direction,
                total: nonnegative_number(&field(&runtime, "total"), &format!("{path}$total"))?,
                initial,
            })
        })
        .collect()
}

/// Parse changing fixed-charge doses against the cached catalog.
///
/// The expected R value is a named list with one entry per active compound:
/// `list(na = 0.001, so4 = 0.0005)`.
/// Omitted configured compounds contribute a zero dose.
fn parse_independent_compounds(
    values: &List,
    definitions: &HashMap<String, i8>,
) -> SolverResult<Vec<IndependentCompound>> {
    named_entries(values, "ph_independent")?
        .into_iter()
        .map(|(name, value)| {
            let charge = definitions
                .get(&name)
                .ok_or_else(|| format!("Unknown compound in `ph_independent`: {name}."))?;
            let path = format!("ph_independent${name}");
            Ok(IndependentCompound {
                charge: *charge,
                dose: nonnegative_number(&value, &path)?,
            })
        })
        .collect()
}

/// Packages solver construction success or failure into the R wrapper protocol.
fn solver_response(result: SolverResult<SolverConfig>) -> List {
    match result {
        Ok(config) => list!(value = ExternalPtr::new(config), error = Robj::from(())),
        Err(error) => list!(value = Robj::from(()), error = error),
    }
}

/// Packages a numerical solve result or error into the R wrapper protocol.
fn solve_response(result: SolverResult<f64>) -> List {
    match result {
        Ok(value) => list!(value = value, error = Robj::from(())),
        Err(error) => list!(value = Robj::from(()), error = error),
    }
}

/// Creates an owning pointer so catalog validation and parsing happen only once.
///
/// @keywords internal
/// @usage NULL
#[extendr]
fn create_solver(ph_dependent: List, ph_independent_charges: List) -> List {
    solver_response(
        validate_catalog(&ph_dependent, &ph_independent_charges).map(
            |(dependent_definitions, independent_definitions)| SolverConfig {
                dependent_definitions,
                independent_definitions,
            },
        ),
    )
}

/// Bridges `solve.solvephrust_solver()` inputs into the validated numerical core.
///
/// @keywords internal
/// @usage NULL
#[extendr]
#[allow(clippy::too_many_arguments)]
fn solve_generic(
    solver: ExternalPtr<SolverConfig>,
    temp: Robj,
    ionic_strength: Robj,
    kw: Robj,
    dependent_compounds: List,
    independent_compounds: List,
    h_i: Robj,
    oh_i: Robj,
) -> List {
    solve_response((|| {
        let temp = scalar_number(&temp, "temp")?;
        if temp <= -273.15 {
            return Err("`temp` must be one finite number > -273.15.".to_string());
        }
        let ionic_strength = optional_nonnegative_number(&ionic_strength, "ionic_strength")?;
        let kw = positive_number(&kw, "kw")?;
        let h_i = nonnegative_number(&h_i, "h_i")?;
        let oh_i = nonnegative_number(&oh_i, "oh_i")?;
        let dependent_compounds =
            parse_dependent_compounds(&dependent_compounds, &solver.dependent_definitions)?;
        let independent_compounds =
            parse_independent_compounds(&independent_compounds, &solver.independent_definitions)?;
        solve_internal(
            temp,
            ionic_strength,
            kw,
            &dependent_compounds,
            &independent_compounds,
            h_i,
            oh_i,
        )
    })())
}

extendr_module! {
    mod solvephrust;
    fn create_solver;
    fn solve_generic;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Supplies reproducible pseudo-random values without adding a dependency.
    struct TestRng(u64);

    impl TestRng {
        /// Advances the generator and returns a value in [0, 1).
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1_u64 << 53) as f64
        }
    }

    /// Returns a nearest-rank percentile from a non-empty sample.
    fn percentile(values: &mut [usize], probability: f64) -> usize {
        values.sort_unstable();
        let index = ((probability * values.len() as f64).ceil() as usize)
            .saturating_sub(1)
            .min(values.len() - 1);
        values[index]
    }

    #[test]
    fn alpha_fractions_sum_to_one_for_every_supported_state_count() {
        for ks in [vec![1e-7], vec![1e-6, 1e-10], vec![1e-2, 1e-7, 1e-12]] {
            for direction in [ChargeDirection::Negative, ChargeDirection::Positive] {
                let alphas = calculate_alphas(1e-8, &ks, direction);
                assert!((alphas.as_slice().iter().sum::<f64>() - 1.0).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn alphas_follow_distance_from_neutral_k_ordering() {
        let ks = [2.0, 2.0, 2.0];
        let negative = calculate_alphas(4.0, &ks, ChargeDirection::Negative);
        let positive = calculate_alphas(4.0, &ks, ChargeDirection::Positive);
        let negative_weights = [1.0, 0.5, 0.25, 0.125];
        let positive_weights = [1.0, 2.0, 4.0, 8.0];
        let negative_total: f64 = negative_weights.iter().sum();
        let positive_total: f64 = positive_weights.iter().sum();
        for index in 0..4 {
            assert!(
                (negative.as_slice()[index] - negative_weights[index] / negative_total).abs()
                    < 1e-14
            );
            assert!(
                (positive.as_slice()[index] - positive_weights[index] / positive_total).abs()
                    < 1e-14
            );
        }
    }

    #[test]
    fn log_and_direct_alpha_calculations_are_equivalent() {
        let cases = [
            (1e-8, vec![1e-7]),
            (1e-8, vec![1e-6, 1e-10]),
            (1e-8, vec![1e-2, 1e-7, 1e-12]),
        ];
        for (h, ks) in cases {
            for direction in [ChargeDirection::Negative, ChargeDirection::Positive] {
                let log_alphas = calculate_alphas(h, &ks, direction);
                let direct_alphas = calculate_alphas_without_logs(h, &ks, direction);
                for (log_alpha, direct_alpha) in
                    log_alphas.as_slice().iter().zip(direct_alphas.as_slice())
                {
                    assert!((log_alpha - direct_alpha).abs() < 1e-14);
                }
            }
        }
    }

    #[test]
    fn positive_k_charge_transitions_start_next_to_neutral() {
        assert_eq!(ion_charges(ChargeDirection::Positive, 0), (1, 0));
        assert_eq!(ion_charges(ChargeDirection::Positive, 1), (2, 1));
        assert_eq!(ion_charges(ChargeDirection::Positive, 2), (3, 2));
    }

    #[test]
    fn charge_direction_selects_opposite_limiting_states() {
        let ions = [Ion {
            k: 1e-7,
            delta_h: 0.0,
        }];
        let negative = DependentCompound {
            ions: &ions,
            charge_direction: ChargeDirection::Negative,
            total: 1.0,
            initial: [0.0; MAX_CONSTANTS],
        };
        let positive = DependentCompound {
            charge_direction: ChargeDirection::Positive,
            ..negative.clone()
        };
        let acidic_negative = calculate_alphas(1.0, &[1e-7], ChargeDirection::Negative);
        let basic_negative = calculate_alphas(1e-14, &[1e-7], ChargeDirection::Negative);
        let acidic_positive = calculate_alphas(1.0, &[1e-7], ChargeDirection::Positive);
        let basic_positive = calculate_alphas(1e-14, &[1e-7], ChargeDirection::Positive);
        assert!(equilibrium_charge(&negative, &acidic_negative).abs() < 1e-6);
        assert!((equilibrium_charge(&negative, &basic_negative) + 1.0).abs() < 1e-6);
        assert!((equilibrium_charge(&positive, &acidic_positive) - 1.0).abs() < 1e-6);
        assert!(equilibrium_charge(&positive, &basic_positive).abs() < 1e-6);
    }

    #[test]
    fn activity_correction_uses_transition_charge_states() {
        let ions = [Ion {
            k: 1e-7,
            delta_h: 0.0,
        }];
        let negative = DependentCompound {
            ions: &ions,
            charge_direction: ChargeDirection::Negative,
            total: 0.0,
            initial: [0.0; MAX_CONSTANTS],
        };
        let positive = DependentCompound {
            charge_direction: ChargeDirection::Positive,
            ..negative.clone()
        };
        let ionic_strength = Some(0.1);
        let gamma_one = calculate_activity(1, ionic_strength, 25.0);
        let coefficients = activity_coefficients(ionic_strength, 25.0);
        assert!((coefficients[2] - calculate_activity(2, ionic_strength, 25.0)).abs() < 1e-15);
        assert!((coefficients[3] - calculate_activity(3, ionic_strength, 25.0)).abs() < 1e-15);
        let negative_k = correct_k(&negative, 25.0, &coefficients)[0];
        let positive_k = correct_k(&positive, 25.0, &coefficients)[0];
        assert!((negative_k - 1e-7 / gamma_one.powi(2)).abs() < 1e-20);
        assert!((positive_k - 1e-7).abs() < 1e-20);

        let dipositive_ions = [
            Ion {
                k: 1e-7,
                delta_h: 0.0,
            },
            Ion {
                k: 1e-9,
                delta_h: 0.0,
            },
        ];
        let dipositive = DependentCompound {
            ions: &dipositive_ions,
            charge_direction: ChargeDirection::Positive,
            total: 0.0,
            initial: [0.0; MAX_CONSTANTS],
        };
        let corrected = correct_k(&dipositive, 25.0, &coefficients);
        let gamma_two = calculate_activity(2, ionic_strength, 25.0);
        assert!((corrected[0] - 1e-7).abs() < 1e-20);
        assert!((corrected[1] - 1e-9 * gamma_two / gamma_one.powi(2)).abs() < 1e-20);
    }

    #[test]
    fn charge_moment_derivative_matches_finite_differences() {
        let constant_sets: [Vec<f64>; 3] = [vec![1e-6], vec![1e-4, 1e-9], vec![1e-3, 1e-7, 1e-12]];
        for ks in constant_sets {
            for direction in [ChargeDirection::Negative, ChargeDirection::Positive] {
                let mut cumulative_log_k = [0.0; MAX_SPECIES];
                for (index, k) in ks.iter().enumerate() {
                    cumulative_log_k[index + 1] = cumulative_log_k[index] + k.ln();
                }
                let compound = PreparedCompound {
                    cumulative_log_k,
                    len: ks.len() + 1,
                    charge_direction: direction,
                    total: 0.013,
                };
                let compounds = [compound];
                let balance = ChargeBalance {
                    compounds: &compounds,
                    kw: 1e-14,
                    gamma_h: 0.81,
                    dose_charge: 0.002,
                    starting_charge: -0.003,
                };
                for ph in [1.0, 5.5, 9.5, 13.0] {
                    let log_h = -ph * std::f64::consts::LN_10;
                    let epsilon = 1e-6;
                    let point = balance.evaluate(log_h);
                    let numerical = (balance.evaluate(log_h + epsilon).balance
                        - balance.evaluate(log_h - epsilon).balance)
                        / (2.0 * epsilon);
                    let scale = point.derivative.abs().max(1e-12);
                    assert!(
                        (point.derivative - numerical).abs() / scale < 1e-6,
                        "direction={direction:?} ks={ks:?} ph={ph} analytic={} numerical={numerical}",
                        point.derivative
                    );
                    assert!(point.derivative > 0.0);
                }
            }
        }
    }

    #[test]
    fn safeguarded_solver_handles_domain_boundaries_and_no_bracket() {
        let ions = [Ion {
            k: 1e-7,
            delta_h: 0.0,
        }];
        let compound = DependentCompound {
            ions: &ions,
            charge_direction: ChargeDirection::Negative,
            total: 0.0,
            initial: [0.0; MAX_CONSTANTS],
        };
        let compounds = [compound];
        let lower_h = 1e-14;
        let lower_oh = 1e-14 / lower_h;
        let lower =
            solve_internal(25.0, None, 1e-14, &compounds, &[], 0.0, lower_oh - lower_h).unwrap();
        assert!((lower - 14.0).abs() < 1e-11);

        let upper_h = 1.0;
        let upper_oh = 1e-14 / upper_h;
        let upper =
            solve_internal(25.0, None, 1e-14, &compounds, &[], upper_h - upper_oh, 0.0).unwrap();
        assert!(upper.abs() < 1e-11);
        assert!(solve_internal(25.0, None, 1e-14, &compounds, &[], 2.0, 0.0).is_err());
    }

    #[test]
    fn root_solvers_agree_for_deterministic_random_chemistry() {
        let mut rng = TestRng(0x5eed_fade_cafe_beef);
        let mut newton_evaluations = Vec::new();
        let mut log_brent_evaluations = Vec::new();
        let mut cached_brent_evaluations = Vec::new();
        let mut legacy_evaluations = Vec::new();
        let mut legacy_accurate_cases = 0;
        let mut largest_legacy_error = 0.0_f64;

        for case_index in 0..2_000 {
            let valence = case_index % MAX_CONSTANTS + 1;
            let direction = if case_index % 2 == 0 {
                ChargeDirection::Negative
            } else {
                ChargeDirection::Positive
            };
            let ions: Vec<Ion> = (0..valence)
                .map(|_| Ion {
                    k: 10_f64.powf(-14.0 + 13.0 * rng.next()),
                    delta_h: -50_000.0 + 100_000.0 * rng.next(),
                })
                .collect();
            let temp = -10.0 + 90.0 * rng.next();
            let ionic_strength = if case_index % 4 == 0 {
                None
            } else {
                Some(0.3 * rng.next())
            };
            let total = 10_f64.powf(-8.0 + 7.0 * rng.next());
            let target_log_h = -(0.05 + 13.9 * rng.next()) * std::f64::consts::LN_10;
            let target_h = target_log_h.exp();
            let gammas = activity_coefficients(ionic_strength, temp);
            let empty_initial = [0.0; MAX_CONSTANTS];
            let uninitialized = DependentCompound {
                ions: &ions,
                charge_direction: direction,
                total,
                initial: empty_initial,
            };
            let corrected = correct_k(&uninitialized, temp, &gammas);
            let alphas = calculate_alphas(target_h, &corrected[..valence], direction);
            let mut initial = [0.0; MAX_CONSTANTS];
            for (index, alpha) in alphas.as_slice()[1..].iter().enumerate() {
                initial[index] = total * alpha;
            }
            let compound = DependentCompound {
                initial,
                ..uninitialized
            };
            let compounds = [compound];
            let target_oh = 1e-14 / (target_h * gammas[1] * gammas[1]);

            let newton_result = solve_internal(
                temp,
                ionic_strength,
                1e-14,
                &compounds,
                &[],
                target_h,
                target_oh,
            )
            .unwrap();
            let expected = -(target_log_h + gammas[1].ln()) / std::f64::consts::LN_10;
            assert!(
                (newton_result - expected).abs() <= 1e-10,
                "case={case_index} newton={newton_result:.16} expected={expected:.16} difference={:.3e}",
                (newton_result - expected).abs()
            );
            let cached = solve_internal_cached_brent(
                temp,
                ionic_strength,
                1e-14,
                &compounds,
                &[],
                target_h,
                target_oh,
            )
            .unwrap();
            let (legacy, legacy_count) = solve_internal_brent_reference(
                temp,
                ionic_strength,
                1e-14,
                &compounds,
                &[],
                target_h,
                target_oh,
            )
            .unwrap();
            assert!((cached.ph - legacy).abs() <= 1e-14);
            let legacy_error = (legacy - expected).abs();
            largest_legacy_error = largest_legacy_error.max(legacy_error);
            if legacy_error <= 1e-10 {
                legacy_accurate_cases += 1;
                assert!((newton_result - legacy).abs() <= 1e-10);
            }

            let (prepared, gamma_h, dose_charge, starting_charge) =
                prepare_balance(temp, ionic_strength, &compounds, &[], target_h, target_oh);
            let balance = ChargeBalance {
                compounds: &prepared,
                kw: 1e-14,
                gamma_h,
                dose_charge,
                starting_charge,
            };
            let newton = find_root_safeguarded(&balance).unwrap();
            let log_brent = find_root_log_brent(&balance).unwrap();
            let newton_ph = -(newton.log_h + gamma_h.ln()) / std::f64::consts::LN_10;
            let log_brent_ph = -(log_brent.log_h + gamma_h.ln()) / std::f64::consts::LN_10;
            assert!((log_brent_ph - expected).abs() <= 1e-10);
            assert!((newton_ph - log_brent_ph).abs() <= 1e-10);
            newton_evaluations.push(newton.evaluations);
            log_brent_evaluations.push(log_brent.evaluations);
            cached_brent_evaluations.push(cached.evaluations);
            legacy_evaluations.push(legacy_count);
        }

        let newton_median = percentile(&mut newton_evaluations, 0.5);
        let newton_p95 = percentile(&mut newton_evaluations, 0.95);
        let newton_max = *newton_evaluations.iter().max().unwrap();
        let log_brent_median = percentile(&mut log_brent_evaluations, 0.5);
        let cached_brent_median = percentile(&mut cached_brent_evaluations, 0.5);
        let legacy_median = percentile(&mut legacy_evaluations, 0.5);
        eprintln!(
            "evaluations: newton median={newton_median} p95={newton_p95} max={newton_max}; log-brent median={log_brent_median}; cached-brent median={cached_brent_median}; legacy median={legacy_median}; legacy accurate={legacy_accurate_cases}/2000 largest error={largest_legacy_error:.3e} pH"
        );
        assert!(legacy_accurate_cases > 1_000);
        assert!(newton_median < log_brent_median);
        assert!(newton_median < legacy_median);
        assert!(cached_brent_median < legacy_median);
    }
}
