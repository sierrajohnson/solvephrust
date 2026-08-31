use extendr_api::prelude::*;
use std::cell::RefCell;
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
#[cfg(test)]
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
    #[cfg(test)]
    k: f64,
    #[cfg(test)]
    delta_h: f64,
    ln_k: f64,
    delta_h_over_r: f64,
}

impl Ion {
    /// Caches log-space temperature terms once when catalog chemistry is parsed.
    fn new(k: f64, delta_h: f64) -> Self {
        Self {
            #[cfg(test)]
            k,
            #[cfg(test)]
            delta_h,
            ln_k: k.ln(),
            delta_h_over_r: delta_h / 8.314,
        }
    }
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

/// Stores one validated fixed-charge catalog entry in indexed solver order.
#[derive(Debug, Clone, Copy)]
struct IndependentDefinition {
    charge: i8,
}

/// Combines cached chemistry with concentrations supplied for one solve call.
#[cfg(test)]
#[derive(Debug, Clone)]
struct DependentCompound<'a> {
    ions: &'a [Ion],
    charge_direction: ChargeDirection,
    total: f64,
    initial: [f64; MAX_CONSTANTS],
}

/// Holds one fixed-charge dose in the form needed by the charge balance.
#[cfg(test)]
#[derive(Debug, Clone)]
struct IndependentCompound {
    charge: i8,
    dose: f64,
}

/// Reuses invocation-sized buffers so successful solves do not allocate.
#[derive(Debug)]
struct SolveScratch {
    prepared: Vec<PreparedCompound>,
    definition_cache: Vec<Option<CachedDefinition>>,
    activity_cache: Option<CachedActivity>,
    seen_dependent: Vec<bool>,
    seen_independent: Vec<bool>,
    dependent_order: Vec<usize>,
    independent_order: Vec<usize>,
}

impl SolveScratch {
    /// Allocates buffers once at solver construction using catalog capacities.
    fn new(dependent_len: usize, independent_len: usize) -> Self {
        Self {
            prepared: Vec::with_capacity(dependent_len),
            definition_cache: vec![None; dependent_len],
            activity_cache: None,
            seen_dependent: vec![false; dependent_len],
            seen_independent: vec![false; independent_len],
            dependent_order: Vec::with_capacity(dependent_len),
            independent_order: Vec::with_capacity(independent_len),
        }
    }

    /// Clears per-call state while retaining every backing allocation.
    fn reset(&mut self) {
        self.prepared.clear();
        self.seen_dependent.fill(false);
        self.seen_independent.fill(false);
    }
}

/// Owns indexed catalogs and reusable scratch behind the external pointer.
#[derive(Debug)]
struct SolverConfig {
    dependent_names: Vec<Box<str>>,
    dependent_definitions: Vec<DependentDefinition>,
    dependent_index: HashMap<Box<str>, usize>,
    independent_names: Vec<Box<str>>,
    independent_definitions: Vec<IndependentDefinition>,
    independent_index: HashMap<Box<str>, usize>,
    scratch: RefCell<SolveScratch>,
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

/// Identifies exact temperature and ionic-strength inputs for chemistry caches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ChemistryKey {
    temp_bits: u64,
    ionic_strength_bits: Option<u64>,
}

impl ChemistryKey {
    /// Preserves exact floating-point inputs without approximate cache matching.
    fn new(temp: f64, ionic_strength: Option<f64>) -> Self {
        Self {
            temp_bits: temp.to_bits(),
            ionic_strength_bits: ionic_strength.map(f64::to_bits),
        }
    }
}

/// Caches activity terms shared by every compound at one solution condition.
#[derive(Debug, Clone, Copy)]
struct CachedActivity {
    key: ChemistryKey,
    gamma_h: f64,
    log_gamma_h: f64,
    log_gammas: [f64; MAX_SPECIES],
}

/// Caches one definition's corrected constants for its most recent condition.
#[derive(Debug, Clone, Copy)]
struct CachedDefinition {
    key: ChemistryKey,
    prepared: PreparedCompound,
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
#[cfg(test)]
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
#[cfg(test)]
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
#[cfg(test)]
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

/// Prepares cached chemistry directly from one validated runtime total.
fn prepare_definition(
    definition: &DependentDefinition,
    temp: f64,
    log_gammas: &[f64; MAX_SPECIES],
) -> PreparedCompound {
    let inverse_temperature_delta = 1.0 / 298.15 - 1.0 / (temp + 273.15);
    let mut cumulative_log_k = [0.0; MAX_SPECIES];
    for (ion_index, ion) in definition.ions.iter().enumerate() {
        let (reactant_charge, product_charge) = ion_charges(definition.charge_direction, ion_index);
        let corrected_log_k = ion.ln_k
            + ion.delta_h_over_r * inverse_temperature_delta
            + log_gammas[usize::from(reactant_charge.unsigned_abs())]
            - log_gammas[1]
            - log_gammas[usize::from(product_charge.unsigned_abs())];
        cumulative_log_k[ion_index + 1] = cumulative_log_k[ion_index] + corrected_log_k;
    }
    PreparedCompound {
        cumulative_log_k,
        len: definition.ions.len() + 1,
        charge_direction: definition.charge_direction,
        total: 0.0,
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
#[cfg(test)]
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
#[cfg(test)]
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

/// Runs Newton from an already prepared charge balance without allocating.
fn solve_prepared(
    kw: f64,
    compounds: &[PreparedCompound],
    gamma_h: f64,
    log_gamma_h: f64,
    dose_charge: f64,
    starting_charge: f64,
) -> SolverResult<f64> {
    let balance = ChargeBalance {
        compounds,
        kw,
        gamma_h,
        dose_charge,
        starting_charge,
    };
    let root = find_root_safeguarded(&balance)?;
    Ok(-(root.log_h + log_gamma_h) / std::f64::consts::LN_10)
}

/// Solves pH with safeguarded Newton iteration in log-hydrogen space.
#[cfg(test)]
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
    solve_prepared(
        kw,
        &prepared,
        gamma_h,
        gamma_h.ln(),
        dose_charge,
        starting_charge,
    )
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
fn finite_scalar(value: &Robj) -> Option<f64> {
    value
        .as_real()
        .or_else(|| value.as_integer().map(f64::from))
        .filter(|number| number.is_finite())
}

/// Accepts R integer or double scalars and adds a path only on failure.
fn scalar_number(value: &Robj, path: &str) -> SolverResult<f64> {
    finite_scalar(value).ok_or_else(|| format!("`{path}` must be one finite number."))
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
                    Ok(Ion::new(
                        positive_number(&field(&ion, "k"), &format!("{ion_path}$k"))?,
                        scalar_number(
                            &field(&ion, "delta_h"),
                            &format!("{ion_path}$delta_h"),
                        )?,
                    ))
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

impl SolverConfig {
    /// Converts validated named catalogs into stable indexed execution tables.
    fn from_catalogs(
        dependent: HashMap<String, DependentDefinition>,
        independent: HashMap<String, i8>,
    ) -> Self {
        let dependent_len = dependent.len();
        let independent_len = independent.len();
        let mut dependent_names = Vec::with_capacity(dependent_len);
        let mut dependent_definitions = Vec::with_capacity(dependent_len);
        let mut dependent_index = HashMap::with_capacity(dependent_len);
        for (index, (name, definition)) in dependent.into_iter().enumerate() {
            let name = name.into_boxed_str();
            dependent_index.insert(name.clone(), index);
            dependent_names.push(name);
            dependent_definitions.push(definition);
        }

        let mut independent_names = Vec::with_capacity(independent_len);
        let mut independent_definitions = Vec::with_capacity(independent_len);
        let mut independent_index = HashMap::with_capacity(independent_len);
        for (index, (name, charge)) in independent.into_iter().enumerate() {
            let name = name.into_boxed_str();
            independent_index.insert(name.clone(), index);
            independent_names.push(name);
            independent_definitions.push(IndependentDefinition { charge });
        }

        Self {
            dependent_names,
            dependent_definitions,
            dependent_index,
            independent_names,
            independent_definitions,
            independent_index,
            scratch: RefCell::new(SolveScratch::new(dependent_len, independent_len)),
        }
    }
}

/// Resolves a runtime name through a cached positional fast path or borrowed lookup.
fn resolve_runtime_index(
    name: &str,
    position: usize,
    names: &[Box<str>],
    index: &HashMap<Box<str>, usize>,
    cached_order: &mut Vec<usize>,
) -> Option<usize> {
    if let Some(&cached) = cached_order.get(position) {
        if names[cached].as_ref() == name {
            return Some(cached);
        }
    }
    let resolved = *index.get(name)?;
    if position < cached_order.len() {
        cached_order[position] = resolved;
    } else {
        cached_order.push(resolved);
    }
    Some(resolved)
}

/// Produces the fixed runtime-schema error only when a malformed entry is found.
fn dependent_schema_error(name: &str) -> String {
    format!("`ph_dependent${name}` must contain exactly `total` and `initial`.")
}

/// Validates one runtime dependent entry without allocating successful values.
fn parse_dependent_values(
    name: &str,
    value: &Robj,
    definition: &DependentDefinition,
) -> SolverResult<(f64, f64)> {
    let runtime = value
        .as_list()
        .ok_or_else(|| format!("`ph_dependent${name}` must be a list."))?;
    if runtime.len() != 2 || runtime.names().is_none() {
        return Err(dependent_schema_error(name));
    }

    let mut total_value = None;
    let mut initial_value = None;
    for (field_name, field_value) in runtime.iter() {
        match field_name {
            "total" if total_value.is_none() => total_value = Some(field_value),
            "initial" if initial_value.is_none() => initial_value = Some(field_value),
            _ => return Err(dependent_schema_error(name)),
        }
    }
    let total_value = total_value.ok_or_else(|| dependent_schema_error(name))?;
    let initial_value = initial_value.ok_or_else(|| dependent_schema_error(name))?;
    let total = finite_scalar(&total_value)
        .filter(|number| *number >= 0.0)
        .ok_or_else(|| format!("`ph_dependent${name}$total` must be one finite number >= 0."))?;
    let initial_values = initial_value
        .as_list()
        .ok_or_else(|| format!("`ph_dependent${name}$initial` must be a list."))?;
    if initial_values.len() != definition.ions.len() {
        return Err(format!(
            "`ph_dependent${name}$initial` must be a list of length {}.",
            definition.ions.len()
        ));
    }

    let sign = f64::from(definition.charge_direction.sign());
    let mut initial_charge = 0.0;
    for (ion_index, initial_value) in initial_values.values().enumerate() {
        let concentration = finite_scalar(&initial_value)
            .filter(|number| *number >= 0.0)
            .ok_or_else(|| {
                format!(
                    "`ph_dependent${name}$initial[[{}]]` must be one finite number >= 0.",
                    ion_index + 1
                )
            })?;
        initial_charge += concentration * sign * (ion_index + 1) as f64;
    }

    Ok((total, initial_charge))
}

/// Reuses activity terms while temperature and ionic strength remain unchanged.
fn activity_for_condition(
    scratch: &mut SolveScratch,
    key: ChemistryKey,
    temp: f64,
    ionic_strength: Option<f64>,
) -> CachedActivity {
    if let Some(cached) = scratch.activity_cache {
        if cached.key == key {
            return cached;
        }
    }
    let gammas = activity_coefficients(ionic_strength, temp);
    let log_gamma_h = gammas[1].ln();
    let cached = CachedActivity {
        key,
        gamma_h: gammas[1],
        log_gamma_h,
        log_gammas: [0.0, log_gamma_h, 4.0 * log_gamma_h, 9.0 * log_gamma_h],
    };
    scratch.activity_cache = Some(cached);
    cached
}

/// Parses changing compounds directly into reusable prepared solver storage.
#[allow(clippy::too_many_arguments)]
fn prepare_runtime_inputs(
    solver: &SolverConfig,
    scratch: &mut SolveScratch,
    temp: f64,
    ionic_strength: Option<f64>,
    dependent_values: &List,
    independent_values: &List,
    h_i: f64,
    oh_i: f64,
) -> SolverResult<(f64, f64, f64, f64)> {
    if dependent_values.is_empty() {
        return Err("`ph_dependent` must contain at least one compound.".to_string());
    }
    if dependent_values.names().is_none() {
        return Err("Every entry in `ph_dependent` must have a non-empty name.".to_string());
    }

    scratch.reset();
    let chemistry_key = ChemistryKey::new(temp, ionic_strength);
    let activity = activity_for_condition(scratch, chemistry_key, temp, ionic_strength);
    let mut starting_charge = h_i - oh_i;
    for (position, (name, value)) in dependent_values.iter().enumerate() {
        if name.is_empty() {
            return Err("Every entry in `ph_dependent` must have a non-empty name.".to_string());
        }
        let definition_index = resolve_runtime_index(
            name,
            position,
            &solver.dependent_names,
            &solver.dependent_index,
            &mut scratch.dependent_order,
        )
        .ok_or_else(|| format!("Unknown compound in `ph_dependent`: {name}."))?;
        if scratch.seen_dependent[definition_index] {
            return Err("Names in `ph_dependent` must be unique.".to_string());
        }
        scratch.seen_dependent[definition_index] = true;
        let (total, initial_charge) = parse_dependent_values(
            name,
            &value,
            &solver.dependent_definitions[definition_index],
        )?;
        let mut prepared = match scratch.definition_cache[definition_index] {
            Some(cached) if cached.key == chemistry_key => cached.prepared,
            _ => {
                let prepared = prepare_definition(
                    &solver.dependent_definitions[definition_index],
                    temp,
                    &activity.log_gammas,
                );
                scratch.definition_cache[definition_index] = Some(CachedDefinition {
                    key: chemistry_key,
                    prepared,
                });
                prepared
            }
        };
        prepared.total = total;
        scratch.prepared.push(prepared);
        starting_charge += initial_charge;
    }
    scratch.dependent_order.truncate(dependent_values.len());

    if !independent_values.is_empty() && independent_values.names().is_none() {
        return Err("Every entry in `ph_independent` must have a non-empty name.".to_string());
    }
    let mut dose_charge = 0.0;
    for (position, (name, value)) in independent_values.iter().enumerate() {
        if name.is_empty() {
            return Err("Every entry in `ph_independent` must have a non-empty name.".to_string());
        }
        let definition_index = resolve_runtime_index(
            name,
            position,
            &solver.independent_names,
            &solver.independent_index,
            &mut scratch.independent_order,
        )
        .ok_or_else(|| format!("Unknown compound in `ph_independent`: {name}."))?;
        if scratch.seen_independent[definition_index] {
            return Err("Names in `ph_independent` must be unique.".to_string());
        }
        scratch.seen_independent[definition_index] = true;
        let dose = finite_scalar(&value)
            .filter(|number| *number >= 0.0)
            .ok_or_else(|| format!("`ph_independent${name}` must be one finite number >= 0."))?;
        dose_charge += f64::from(solver.independent_definitions[definition_index].charge) * dose;
    }
    scratch.independent_order.truncate(independent_values.len());

    Ok((
        activity.gamma_h,
        activity.log_gamma_h,
        dose_charge,
        starting_charge,
    ))
}

/// Creates an owning pointer so catalog validation and parsing happen only once.
///
/// @keywords internal
/// @usage NULL
#[extendr]
fn create_solver(ph_dependent: List, ph_independent_charges: List) -> ExternalPtr<SolverConfig> {
    let result =
        validate_catalog(&ph_dependent, &ph_independent_charges).map(|(dependent, independent)| {
            ExternalPtr::new(SolverConfig::from_catalogs(dependent, independent))
        });
    match result {
        Ok(solver) => solver,
        Err(error) => throw_r_error(error),
    }
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
) -> f64 {
    let result = (|| {
        let temp = scalar_number(&temp, "temp")?;
        if temp <= -273.15 {
            return Err("`temp` must be one finite number > -273.15.".to_string());
        }
        let ionic_strength = optional_nonnegative_number(&ionic_strength, "ionic_strength")?;
        let kw = positive_number(&kw, "kw")?;
        let h_i = nonnegative_number(&h_i, "h_i")?;
        let oh_i = nonnegative_number(&oh_i, "oh_i")?;
        let mut scratch = solver.scratch.try_borrow_mut().map_err(|_| {
            "A solvephrust solver cannot be used reentrantly from the same R call.".to_string()
        })?;
        let (gamma_h, log_gamma_h, dose_charge, starting_charge) = prepare_runtime_inputs(
            &solver,
            &mut scratch,
            temp,
            ionic_strength,
            &dependent_compounds,
            &independent_compounds,
            h_i,
            oh_i,
        )?;
        solve_prepared(
            kw,
            &scratch.prepared,
            gamma_h,
            log_gamma_h,
            dose_charge,
            starting_charge,
        )
    })();
    match result {
        Ok(value) => value,
        Err(error) => throw_r_error(error),
    }
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
        let ions = [Ion::new(1e-7, 0.0)];
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
        let ions = [Ion::new(1e-7, 0.0)];
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

        let dipositive_ions = [Ion::new(1e-7, 0.0), Ion::new(1e-9, 0.0)];
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
        let ions = [Ion::new(1e-7, 0.0)];
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
                .map(|_| {
                    Ion::new(
                        10_f64.powf(-14.0 + 13.0 * rng.next()),
                        -50_000.0 + 100_000.0 * rng.next(),
                    )
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
            legacy_evaluations.push(legacy_count);
        }

        let newton_median = percentile(&mut newton_evaluations, 0.5);
        let newton_p95 = percentile(&mut newton_evaluations, 0.95);
        let newton_max = *newton_evaluations.iter().max().unwrap();
        let log_brent_median = percentile(&mut log_brent_evaluations, 0.5);
        let legacy_median = percentile(&mut legacy_evaluations, 0.5);
        eprintln!(
            "evaluations: newton median={newton_median} p95={newton_p95} max={newton_max}; log-brent median={log_brent_median}; legacy median={legacy_median}; legacy accurate={legacy_accurate_cases}/2000 largest error={largest_legacy_error:.3e} pH"
        );
        assert!(legacy_accurate_cases > 1_000);
        assert!(newton_median < log_brent_median);
        assert!(newton_median < legacy_median);
    }
}
