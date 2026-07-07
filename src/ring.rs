//! Polynomial ring helpers for parsing and formatting user-facing polynomials.

use crate::field::Field;
use crate::finite_field::PrimeField;
use crate::monomial::{Monomial, MonomialOrder};
use crate::polynomial::{Polynomial, Term};
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
        }
    }
}

impl std::error::Error for ParsePolynomialError {}

#[derive(Debug, Clone)]
pub struct PolynomialRing<F> {
    variables: Vec<String>,
    variable_indices: HashMap<String, usize>,
    order: MonomialOrder,
    _field: PhantomData<F>,
}

impl<F> PolynomialRing<F> {
    pub fn new<I, S>(variables: I, order: MonomialOrder) -> Result<Self, ParsePolynomialError>
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

        Ok(Self {
            variables,
            variable_indices,
            order,
            _field: PhantomData,
        })
    }

    pub fn variables(&self) -> &[String] {
        &self.variables
    }

    pub fn order(&self) -> MonomialOrder {
        self.order
    }
}

impl<F: ParseCoefficient> PolynomialRing<F> {
    pub fn parse(&self, input: &str) -> Result<Polynomial<F>, ParsePolynomialError> {
        let tokens = Lexer::new(input).tokenize()?;
        let mut parser = ParserFor::new(tokens, self);
        let polynomial = parser.parse_polynomial()?;
        if let Some(token) = parser.peek() {
            return Err(ParsePolynomialError::TrailingInput(token.to_string()));
        }
        Ok(polynomial)
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
            let coeff = term.coefficient.to_string();
            let is_negative = coeff.starts_with('-');
            let abs_coeff = if is_negative { &coeff[1..] } else { &coeff };
            let monomial = self.format_monomial(&term.monomial)?;
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

    pub fn format_monomial(&self, monomial: &Monomial) -> Result<String, ParsePolynomialError> {
        if monomial.nvars() != self.variables.len() {
            return Err(ParsePolynomialError::WrongVariableCount {
                expected: self.variables.len(),
                actual: monomial.nvars(),
            });
        }

        let mut factors = Vec::new();
        for (variable, exponent) in self.variables.iter().zip(monomial.exponents.iter()) {
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
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
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

struct ParserFor<'a, F> {
    tokens: Vec<Token>,
    index: usize,
    ring: &'a PolynomialRing<F>,
}

impl<'a, F: ParseCoefficient> ParserFor<'a, F> {
    fn new(tokens: Vec<Token>, ring: &'a PolynomialRing<F>) -> Self {
        Self {
            tokens,
            index: 0,
            ring,
        }
    }

    fn parse_polynomial(&mut self) -> Result<Polynomial<F>, ParsePolynomialError> {
        let mut terms = Vec::new();
        let mut saw_term = false;

        while self.peek().is_some() {
            let sign = self.parse_sign();
            if self.peek().is_none() {
                return Err(ParsePolynomialError::ExpectedTerm);
            }
            terms.push(self.parse_term(sign)?);
            saw_term = true;

            if !matches!(self.peek(), Some(Token::Plus | Token::Minus) | None) {
                return Err(ParsePolynomialError::TrailingInput(
                    self.peek()
                        .map(ToString::to_string)
                        .unwrap_or_else(|| String::from("<end>")),
                ));
            }
        }

        if !saw_term {
            return Err(ParsePolynomialError::ExpectedTerm);
        }

        Ok(Polynomial::new(
            terms,
            self.ring.variables.len(),
            self.ring.order,
        ))
    }

    fn parse_sign(&mut self) -> i32 {
        match self.peek() {
            Some(Token::Plus) => {
                self.index += 1;
                1
            }
            Some(Token::Minus) => {
                self.index += 1;
                -1
            }
            _ => 1,
        }
    }

    fn parse_term(&mut self, sign: i32) -> Result<Term<F>, ParsePolynomialError> {
        let mut coefficient = if sign < 0 { F::minus_one() } else { F::one() };
        let mut exponents = vec![0u32; self.ring.variables.len()];
        let mut saw_factor = false;
        let mut require_factor = false;

        while let Some(token) = self.peek().cloned() {
            match token {
                Token::Star => {
                    if !saw_factor || require_factor {
                        return Err(ParsePolynomialError::ExpectedFactor);
                    }
                    self.index += 1;
                    require_factor = true;
                }
                Token::Number(number) => {
                    self.index += 1;
                    let factor = self.parse_coefficient(&number)?;
                    coefficient = coefficient.multiply(&factor);
                    saw_factor = true;
                    require_factor = false;
                }
                Token::Ident(variable) => {
                    self.index += 1;
                    let index = self
                        .ring
                        .variable_indices
                        .get(&variable)
                        .copied()
                        .ok_or(ParsePolynomialError::UnknownVariable(variable))?;
                    let exponent = self.parse_optional_exponent()?;
                    exponents[index] = exponents[index].checked_add(exponent).ok_or_else(|| {
                        ParsePolynomialError::InvalidExponent(exponent.to_string())
                    })?;
                    saw_factor = true;
                    require_factor = false;
                }
                Token::Plus | Token::Minus | Token::Slash | Token::Caret => break,
            }
        }

        if require_factor {
            return Err(ParsePolynomialError::ExpectedFactor);
        }
        if !saw_factor {
            return Err(ParsePolynomialError::ExpectedTerm);
        }

        Ok(Term::new(coefficient, Monomial::new(exponents)))
    }

    fn parse_coefficient(&mut self, numerator: &str) -> Result<F, ParsePolynomialError> {
        let denominator = if matches!(self.peek(), Some(Token::Slash)) {
            self.index += 1;
            let denominator = match self.peek().cloned() {
                Some(Token::Number(number)) => {
                    self.index += 1;
                    number
                }
                _ => return Err(ParsePolynomialError::ExpectedDenominator),
            };
            if denominator == "0" {
                return Err(ParsePolynomialError::DivisionByZero);
            }
            Some(denominator)
        } else {
            None
        };

        F::parse_coefficient(numerator, denominator.as_deref())
    }

    fn parse_optional_exponent(&mut self) -> Result<u32, ParsePolynomialError> {
        if !matches!(self.peek(), Some(Token::Caret)) {
            return Ok(1);
        }

        self.index += 1;
        let exponent = match self.peek().cloned() {
            Some(Token::Number(number)) => {
                self.index += 1;
                number
            }
            _ => return Err(ParsePolynomialError::ExpectedExponent),
        };

        exponent
            .parse::<u32>()
            .map_err(|_| ParsePolynomialError::InvalidExponent(exponent))
    }

    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.index)
    }
}

pub trait ParseCoefficient: Field {
    fn minus_one() -> Self {
        Self::one().negate()
    }

    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
    ) -> Result<Self, ParsePolynomialError>;
}

impl ParseCoefficient for BigRational {
    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
    ) -> Result<Self, ParsePolynomialError> {
        let text = if let Some(denominator) = denominator {
            format!("{numerator}/{denominator}")
        } else {
            numerator.to_string()
        };

        BigRational::from_str(&text).map_err(|_| ParsePolynomialError::InvalidNumber(text))
    }
}

impl<const P: u32> ParseCoefficient for PrimeField<P> {
    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
    ) -> Result<Self, ParsePolynomialError> {
        let numerator = PrimeField::<P>::parse_digits_mod(numerator)
            .map_err(|_| ParsePolynomialError::InvalidNumber(numerator.to_string()))?;
        let Some(denominator) = denominator else {
            return Ok(numerator);
        };
        let denominator = PrimeField::<P>::parse_digits_mod(denominator)
            .map_err(|_| ParsePolynomialError::InvalidNumber(denominator.to_string()))?;
        let inverse = denominator
            .inverse()
            .ok_or(ParsePolynomialError::DivisionByZero)?;
        Ok(numerator.multiply(&inverse))
    }
}

fn is_ident_start(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '_'
}

fn is_ident_continue(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}
