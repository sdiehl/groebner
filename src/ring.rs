//! Polynomial ring helpers for parsing and formatting user-facing polynomials.
//!
//! A `PolynomialRing` stores variable names, their order, and the monomial order used by parsed
//! polynomials. It can format results in plain text or LaTeX using those variable names.
//!
//! # Example
//! ```
//! use groebner::{MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x0", "y"], MonomialOrder::Lex)?;
//! let polynomial = ring.parse("3/2*x0^2*y - y")?;
//!
//! assert_eq!(ring.format(&polynomial)?, "3/2*x0^2*y - y");
//! assert_eq!(
//!     ring.format_latex(&polynomial)?,
//!     "\\frac{3}{2} x_{0}^{2} y - y"
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::Field;
use crate::finite_field::{PrimeField, Zp};
use crate::monomial::{Monomial, MonomialOrder};
use crate::polynomial::Polynomial;
use crate::rational_function::RationalFunction;
use num_rational::BigRational;
use std::collections::HashMap;
use std::fmt;
use std::marker::PhantomData;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParsePolynomialError {
    EmptyVariableList,
    DuplicateVariable(String),
    UnknownVariable(String),
    ExpectedTerm,
    ExpectedFactor,
    ExpectedExponent,
    ExpectedDenominator,
    InvalidNumber(String),
    InvalidExponent(String),
    DivisionByZero,
    TrailingInput(String),
    WrongVariableCount { expected: usize, actual: usize },
    MissingModulus,
}

impl fmt::Display for ParsePolynomialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParsePolynomialError::EmptyVariableList => write!(f, "variable list cannot be empty"),
            ParsePolynomialError::DuplicateVariable(var) => {
                write!(f, "duplicate variable name `{var}`")
            }
            ParsePolynomialError::UnknownVariable(var) => write!(f, "unknown variable `{var}`"),
            ParsePolynomialError::ExpectedTerm => write!(f, "expected polynomial term"),
            ParsePolynomialError::ExpectedFactor => write!(f, "expected term factor"),
            ParsePolynomialError::ExpectedExponent => write!(f, "expected exponent after `^`"),
            ParsePolynomialError::ExpectedDenominator => {
                write!(f, "expected denominator after `/`")
            }
            ParsePolynomialError::InvalidNumber(number) => write!(f, "invalid number `{number}`"),
            ParsePolynomialError::InvalidExponent(exponent) => {
                write!(f, "invalid exponent `{exponent}`")
            }
            ParsePolynomialError::DivisionByZero => {
                write!(f, "rational denominator cannot be zero")
            }
            ParsePolynomialError::TrailingInput(input) => {
                write!(f, "unexpected input after polynomial: `{input}`")
            }
            ParsePolynomialError::WrongVariableCount { expected, actual } => write!(
                f,
                "polynomial has {actual} variables, but this ring has {expected}"
            ),
            ParsePolynomialError::MissingModulus => {
                write!(f, "ring has no modulus; use PolynomialRing::with_modulus")
            }
        }
    }
}

impl std::error::Error for ParsePolynomialError {}

#[derive(Debug, Clone)]
pub struct PolynomialRing<F> {
    variables: Vec<String>,
    variable_indices: HashMap<String, usize>,
    order: MonomialOrder,
    modulus: Option<u64>,
    parameters: Vec<String>,
    _field: PhantomData<F>,
}

impl<F> PolynomialRing<F> {
    pub fn new<I, S>(variables: I, order: MonomialOrder) -> Result<Self, ParsePolynomialError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self::build(variables, order, None, Vec::new())
    }

    /// A ring over a runtime prime field, needed to parse [`Zp`] coefficients.
    pub fn with_modulus<I, S>(
        variables: I,
        order: MonomialOrder,
        modulus: u64,
    ) -> Result<Self, ParsePolynomialError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self::build(variables, order, Some(modulus), Vec::new())
    }

    /// A ring over [`RationalFunction`] whose coefficients may mention `parameter`.
    pub fn with_parameter<I, S>(
        variables: I,
        order: MonomialOrder,
        parameter: impl Into<String>,
    ) -> Result<Self, ParsePolynomialError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self::build(variables, order, None, vec![parameter.into()])
    }

    /// A ring whose coefficients may mention several `parameters`, for `Frac` coefficients
    /// (feature `parameters`).
    pub fn with_parameters<I, S, P, T>(
        variables: I,
        order: MonomialOrder,
        parameters: P,
    ) -> Result<Self, ParsePolynomialError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
        P: IntoIterator<Item = T>,
        T: Into<String>,
    {
        Self::build(
            variables,
            order,
            None,
            parameters.into_iter().map(Into::into).collect(),
        )
    }

    fn build<I, S>(
        variables: I,
        order: MonomialOrder,
        modulus: Option<u64>,
        parameters: Vec<String>,
    ) -> Result<Self, ParsePolynomialError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let variables: Vec<String> = variables.into_iter().map(Into::into).collect();
        if variables.is_empty() {
            return Err(ParsePolynomialError::EmptyVariableList);
        }

        let mut variable_indices = HashMap::with_capacity(variables.len());
        for (index, variable) in variables.iter().enumerate() {
            if variable_indices.insert(variable.clone(), index).is_some() {
                return Err(ParsePolynomialError::DuplicateVariable(variable.clone()));
            }
        }
        let mut seen = std::collections::HashSet::new();
        if let Some(name) = parameters
            .iter()
            .find(|p| variable_indices.contains_key(*p) || !seen.insert(*p))
        {
            return Err(ParsePolynomialError::DuplicateVariable(name.clone()));
        }

        Ok(Self {
            variables,
            variable_indices,
            order,
            modulus,
            parameters,
            _field: PhantomData,
        })
    }

    pub fn variables(&self) -> &[String] {
        &self.variables
    }

    pub fn order(&self) -> MonomialOrder {
        self.order.clone()
    }

    pub fn modulus(&self) -> Option<u64> {
        self.modulus
    }

    /// The first parameter.
    pub fn parameter(&self) -> Option<&str> {
        self.parameters.first().map(String::as_str)
    }

    pub fn parameters(&self) -> &[String] {
        &self.parameters
    }

    /// The same ring under a different monomial order.
    pub fn with_order(&self, order: MonomialOrder) -> Self {
        Self {
            variables: self.variables.clone(),
            variable_indices: self.variable_indices.clone(),
            order,
            modulus: self.modulus,
            parameters: self.parameters.clone(),
            _field: PhantomData,
        }
    }
}

impl<F: ParseCoefficient> PolynomialRing<F> {
    /// Each parameter the field can represent, with its value.
    fn params<G: ParseCoefficient>(&self) -> Vec<(String, G)> {
        let n = self.parameters.len();
        self.parameters
            .iter()
            .enumerate()
            .filter_map(|(i, p)| G::parameter(i, n).map(|v| (p.clone(), v)))
            .collect()
    }

    pub fn parse(&self, input: &str) -> Result<Polynomial<F>, ParsePolynomialError> {
        let tokens = Lexer::new(input).tokenize()?;
        // Retain the legacy grouped-number and implicit-product syntax while delegating
        // expression parsing and expansion to polycore.
        let mut source = String::new();
        let mut previous: Option<&Token> = None;
        for token in &tokens {
            if let Token::Ident(name) = token
                && !self.variable_indices.contains_key(name)
                && !self.params::<F>().iter().any(|(p, _)| p == name)
            {
                return Err(ParsePolynomialError::UnknownVariable(name.clone()));
            }
            let ends = matches!(
                previous,
                Some(Token::Number(_) | Token::Ident(_) | Token::Close)
            );
            let starts = matches!(token, Token::Number(_) | Token::Ident(_) | Token::Open);
            if ends && starts {
                source.push('*');
            }
            if matches!(token, Token::Plus) && matches!(previous, None | Some(Token::Open)) {
                previous = Some(token);
                continue;
            }
            source.push_str(&token.to_string());
            previous = Some(token);
        }
        let core = polycore::Ring::try_new(self.variables.clone(), self.order.clone())
            .map_err(|e| ParsePolynomialError::TrailingInput(e.to_string()))?;
        // Validate the coefficient context before the infallible integer-lifting callback.
        F::parse_coefficient("1", None, self.modulus)?;
        let error = std::cell::RefCell::new(None);
        let lift = |n: num_bigint::BigInt| {
            F::parse_coefficient(&n.to_string(), None, self.modulus).unwrap_or_else(|e| {
                *error.borrow_mut() = Some(e);
                F::zero()
            })
        };
        let params = self.params::<F>();
        let params: Vec<(&str, F)> = params
            .iter()
            .map(|(p, v)| (p.as_str(), v.clone()))
            .collect();
        let result = core.parse_with(&source, &lift, &params);
        if let Some(e) = error.into_inner() {
            return Err(e);
        }
        result
            .map(|p| p.map(|c| c.clone().bind_modulus(self.modulus)))
            .map_err(|e| ParsePolynomialError::TrailingInput(e.to_string()))
    }

    pub fn parse_many(&self, input: &str) -> Result<Vec<Polynomial<F>>, ParsePolynomialError> {
        input
            .split([';', ','])
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(|line| self.parse(line))
            .collect()
    }

    pub fn format(&self, polynomial: &Polynomial<F>) -> Result<String, ParsePolynomialError> {
        if polynomial.nvars != self.variables.len() {
            return Err(ParsePolynomialError::WrongVariableCount {
                expected: self.variables.len(),
                actual: polynomial.nvars,
            });
        }

        if polynomial.is_zero() {
            return Ok("0".to_string());
        }

        let mut output = String::new();
        for term in &polynomial.terms {
            let coeff = term.1.format_coefficient_in(&self.parameters);
            let is_negative = coeff.starts_with('-');
            let abs_coeff = if is_negative { &coeff[1..] } else { &coeff };
            let monomial = self.format_monomial(&term.0)?;
            let term_body = if monomial == "1" {
                abs_coeff.to_string()
            } else if abs_coeff == "1" {
                monomial
            } else {
                format!("{abs_coeff}*{monomial}")
            };

            if output.is_empty() {
                if is_negative {
                    output.push('-');
                }
                output.push_str(&term_body);
            } else if is_negative {
                output.push_str(" - ");
                output.push_str(&term_body);
            } else {
                output.push_str(" + ");
                output.push_str(&term_body);
            }
        }
        Ok(output)
    }

    pub fn format_latex(&self, polynomial: &Polynomial<F>) -> Result<String, ParsePolynomialError> {
        if polynomial.nvars != self.variables.len() {
            return Err(ParsePolynomialError::WrongVariableCount {
                expected: self.variables.len(),
                actual: polynomial.nvars,
            });
        }

        if polynomial.is_zero() {
            return Ok("0".to_string());
        }

        let mut output = String::new();
        for term in &polynomial.terms {
            let names: Vec<String> = self.parameters.iter().map(|p| latex_variable(p)).collect();
            let coeff = term.1.format_coefficient_latex_in(&names);
            let is_negative = coeff.starts_with('-');
            let coeff_latex = if is_negative { &coeff[1..] } else { &coeff }.to_string();
            let abs_coeff = coeff_latex.as_str();
            let monomial = self.format_monomial_latex(&term.0)?;
            let term_body = if monomial == "1" {
                coeff_latex
            } else if abs_coeff == "1" {
                monomial
            } else {
                format!("{coeff_latex} {monomial}")
            };

            if output.is_empty() {
                if is_negative {
                    output.push('-');
                }
                output.push_str(&term_body);
            } else if is_negative {
                output.push_str(" - ");
                output.push_str(&term_body);
            } else {
                output.push_str(" + ");
                output.push_str(&term_body);
            }
        }

        Ok(output)
    }

    pub fn format_monomial(&self, monomial: &Monomial) -> Result<String, ParsePolynomialError> {
        if monomial.nvars() != self.variables.len() {
            return Err(ParsePolynomialError::WrongVariableCount {
                expected: self.variables.len(),
                actual: monomial.nvars(),
            });
        }

        let mut factors = Vec::new();
        for (variable, exponent) in self.variables.iter().zip(monomial.exps().iter()) {
            match exponent {
                0 => {}
                1 => factors.push(variable.clone()),
                exp => factors.push(format!("{variable}^{exp}")),
            }
        }

        if factors.is_empty() {
            Ok("1".to_string())
        } else {
            Ok(factors.join("*"))
        }
    }

    pub fn format_monomial_latex(
        &self,
        monomial: &Monomial,
    ) -> Result<String, ParsePolynomialError> {
        if monomial.nvars() != self.variables.len() {
            return Err(ParsePolynomialError::WrongVariableCount {
                expected: self.variables.len(),
                actual: monomial.nvars(),
            });
        }

        let mut factors = Vec::new();
        for (variable, exponent) in self.variables.iter().zip(monomial.exps().iter()) {
            match exponent {
                0 => {}
                1 => factors.push(latex_variable(variable)),
                exp => factors.push(format!("{}^{{{exp}}}", latex_variable(variable))),
            }
        }

        if factors.is_empty() {
            Ok("1".to_string())
        } else {
            Ok(factors.join(" "))
        }
    }
}

fn latex_coefficient(coefficient: &str) -> String {
    if let Some((numerator, denominator)) = coefficient.split_once('/') {
        format!("\\frac{{{numerator}}}{{{denominator}}}")
    } else {
        coefficient.to_string()
    }
}

fn latex_variable(variable: &str) -> String {
    let split = variable
        .char_indices()
        .rev()
        .find(|(_, ch)| !ch.is_ascii_digit())
        .map(|(index, ch)| index + ch.len_utf8())
        .unwrap_or(0);
    let (base, digits) = variable.split_at(split);
    let escaped_base = base.replace('_', "\\_");
    if digits.is_empty() || base.is_empty() {
        escaped_base
    } else {
        format!("{escaped_base}_{{{digits}}}")
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    Open,
    Close,
    Number(String),
    Ident(String),
}

impl fmt::Display for Token {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Token::Plus => write!(f, "+"),
            Token::Minus => write!(f, "-"),
            Token::Star => write!(f, "*"),
            Token::Slash => write!(f, "/"),
            Token::Caret => write!(f, "^"),
            Token::Open => write!(f, "("),
            Token::Close => write!(f, ")"),
            Token::Number(number) => write!(f, "{number}"),
            Token::Ident(ident) => write!(f, "{ident}"),
        }
    }
}

struct Lexer<'a> {
    input: &'a str,
    chars: Vec<char>,
    index: usize,
}

impl<'a> Lexer<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            input,
            chars: input.chars().collect(),
            index: 0,
        }
    }

    fn tokenize(mut self) -> Result<Vec<Token>, ParsePolynomialError> {
        let mut tokens = Vec::new();
        while self.index < self.chars.len() {
            match self.chars[self.index] {
                c if c.is_whitespace() => self.index += 1,
                '+' => {
                    self.index += 1;
                    tokens.push(Token::Plus);
                }
                '-' => {
                    self.index += 1;
                    tokens.push(Token::Minus);
                }
                '*' => {
                    self.index += 1;
                    tokens.push(Token::Star);
                }
                '/' => {
                    self.index += 1;
                    tokens.push(Token::Slash);
                }
                '^' => {
                    self.index += 1;
                    tokens.push(Token::Caret);
                }
                '(' => {
                    self.index += 1;
                    tokens.push(Token::Open);
                }
                ')' => {
                    self.index += 1;
                    tokens.push(Token::Close);
                }
                ',' | ';' => {
                    return Err(ParsePolynomialError::TrailingInput(
                        self.chars[self.index].to_string(),
                    ));
                }
                c if c.is_ascii_digit() => tokens.push(Token::Number(self.read_number())),
                c if is_ident_start(c) => tokens.push(Token::Ident(self.read_ident())),
                _ => {
                    return Err(ParsePolynomialError::TrailingInput(
                        self.input.chars().skip(self.index).collect(),
                    ));
                }
            }
        }
        Ok(tokens)
    }

    fn read_number(&mut self) -> String {
        let mut number = self.read_ascii_digits();

        loop {
            let after_digits = self.index;
            self.skip_inline_whitespace();
            if self.index >= self.chars.len() || !self.chars[self.index].is_ascii_digit() {
                self.index = after_digits;
                break;
            }

            let group = self.read_ascii_digits();
            if group.len() == 3 {
                number.push_str(&group);
            } else {
                self.index = after_digits;
                break;
            }
        }

        number
    }

    fn read_ascii_digits(&mut self) -> String {
        let start = self.index;
        while self.index < self.chars.len() && self.chars[self.index].is_ascii_digit() {
            self.index += 1;
        }
        self.chars[start..self.index].iter().collect()
    }

    fn skip_inline_whitespace(&mut self) {
        while self.index < self.chars.len() && self.chars[self.index].is_whitespace() {
            self.index += 1;
        }
    }

    fn read_ident(&mut self) -> String {
        let start = self.index;
        self.index += 1;
        while self.index < self.chars.len() && is_ident_continue(self.chars[self.index]) {
            self.index += 1;
        }
        self.chars[start..self.index].iter().collect()
    }
}

/// Coefficient parsing and printing for [`PolynomialRing`]; implement for custom fields.
pub trait ParseCoefficient: Field + fmt::Display {
    /// The ring parameter raised to `exponent`, for fields with a parameter.
    fn parameter_power(_exponent: u32) -> Option<Self> {
        None
    }

    /// Parameter `index` of `count`, for fields with parameters. The default is
    /// [`Self::parameter_power`] for a single parameter.
    fn parameter(index: usize, count: usize) -> Option<Self> {
        (index == 0 && count == 1)
            .then(|| Self::parameter_power(1))
            .flatten()
    }

    /// Text for this coefficient as a factor of a term, with the parameters named `parameters`.
    fn format_coefficient_in(&self, parameters: &[String]) -> String {
        self.format_coefficient(parameters.first().map_or("a", String::as_str))
    }

    /// LaTeX for this coefficient, with the parameters already in LaTeX.
    fn format_coefficient_latex_in(&self, parameters: &[String]) -> String {
        self.format_coefficient_latex(parameters.first().map_or("a", String::as_str))
    }

    /// Text for this coefficient as a factor of a term, with the parameter named `parameter`.
    fn format_coefficient(&self, _parameter: &str) -> String {
        self.to_string()
    }

    fn format_coefficient_latex(&self, parameter: &str) -> String {
        let text = self.format_coefficient(parameter);
        match text.strip_prefix('-') {
            Some(rest) => format!("-{}", latex_coefficient(rest)),
            None => latex_coefficient(&text),
        }
    }

    fn minus_one() -> Self {
        -Self::one()
    }

    fn bind_modulus(self, _modulus: Option<u64>) -> Self {
        self
    }

    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
        modulus: Option<u64>,
    ) -> Result<Self, ParsePolynomialError>;
}

impl ParseCoefficient for BigRational {
    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
        _modulus: Option<u64>,
    ) -> Result<Self, ParsePolynomialError> {
        let text = if let Some(denominator) = denominator {
            format!("{numerator}/{denominator}")
        } else {
            numerator.to_string()
        };
        BigRational::from_str(&text).map_err(|_| ParsePolynomialError::InvalidNumber(text))
    }
}

fn modular_quotient<F: Field>(
    numerator: F,
    denominator: Option<Result<F, ParsePolynomialError>>,
) -> Result<F, ParsePolynomialError> {
    let Some(denominator) = denominator else {
        return Ok(numerator);
    };
    let inverse = denominator?
        .inverse()
        .ok_or(ParsePolynomialError::DivisionByZero)?;
    Ok(numerator * inverse)
}

impl<const P: u64> ParseCoefficient for PrimeField<P> {
    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
        _modulus: Option<u64>,
    ) -> Result<Self, ParsePolynomialError> {
        let parse = |digits: &str| {
            parse_digits(digits, P)
                .map(PrimeField::<P>::new)
                .map_err(|_| ParsePolynomialError::InvalidNumber(digits.to_string()))
        };
        modular_quotient(parse(numerator)?, denominator.map(parse))
    }
}

impl ParseCoefficient for Zp {
    fn bind_modulus(self, modulus: Option<u64>) -> Self {
        modulus.map_or(self, |m| self.bind(m))
    }

    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
        modulus: Option<u64>,
    ) -> Result<Self, ParsePolynomialError> {
        let modulus = modulus.ok_or(ParsePolynomialError::MissingModulus)?;
        let parse = |digits: &str| {
            parse_digits(digits, modulus)
                .map(|v| Zp::new(v, modulus))
                .map_err(|_| ParsePolynomialError::InvalidNumber(digits.to_string()))
        };
        modular_quotient(parse(numerator)?, denominator.map(parse))
    }
}

impl ParseCoefficient for RationalFunction {
    fn parameter_power(exponent: u32) -> Option<Self> {
        Some(Self::parameter_power(exponent))
    }

    fn format_coefficient(&self, parameter: &str) -> String {
        self.render(parameter, true)
    }

    fn format_coefficient_latex(&self, parameter: &str) -> String {
        self.render_latex(parameter)
    }

    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
        modulus: Option<u64>,
    ) -> Result<Self, ParsePolynomialError> {
        BigRational::parse_coefficient(numerator, denominator, modulus).map(Self::constant)
    }
}

fn is_ident_start(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '_'
}

fn is_ident_continue(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

fn parse_digits(digits: &str, p: u64) -> Result<u64, ()> {
    if p < 2 || digits.is_empty() {
        return Err(());
    }
    digits.bytes().try_fold(0u64, |v, c| {
        if !c.is_ascii_digit() {
            return Err(());
        }
        Ok(((u128::from(v) * 10 + u128::from(c - b'0')) % u128::from(p)) as u64)
    })
}
