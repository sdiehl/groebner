//! Stress runner for the public benchmark corpus in `tests/corpus`.
//!
//! Every system runs in a child process so that timeouts, crashes and runaway memory are
//! isolated, and its basis is checked against its reference fingerprint.
//!
//! ```text
//! cargo run --release --example corpus -- [--filter katsura] [--field zp,qq,gf32003]
//!     [--algo f4,buchberger] [--timeout 60] [--jobs 1] [--out results.tsv] [--dir tests/corpus]
//! ```

#[path = "../tests/common/corpus.rs"]
mod corpus;

use corpus::{Algorithm, Field, System};
use std::env;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::{Duration, Instant};

struct Options {
    dir: PathBuf,
    filters: Vec<String>,
    fields: Vec<Field>,
    algorithms: Vec<Algorithm>,
    timeout: Duration,
    jobs: usize,
    out: Option<PathBuf>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Status {
    Ok,
    Mismatch,
    Error,
    Crash,
    Timeout,
}

impl Status {
    fn label(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::Mismatch => "MISMATCH",
            Self::Error => "ERROR",
            Self::Crash => "CRASH",
            Self::Timeout => "timeout",
        }
    }
}

struct Outcome {
    name: String,
    field: Field,
    algorithm: Algorithm,
    status: Status,
    seconds: f64,
    reference: f64,
    size: usize,
    detail: String,
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.first().is_some_and(|a| a == "--one") {
        return child(&args[1..]);
    }
    match parse_options(&args) {
        Ok(options) => parent(&options),
        Err(message) => {
            eprintln!("{message}");
            ExitCode::from(2)
        }
    }
}

fn parse_options(args: &[String]) -> Result<Options, String> {
    let mut options = Options {
        dir: PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/corpus")),
        filters: Vec::new(),
        fields: vec![Field::Zp],
        algorithms: vec![Algorithm::F4],
        timeout: Duration::from_secs(60),
        jobs: 1,
        out: None,
    };
    let mut iter = args.iter();
    while let Some(flag) = iter.next() {
        let value = iter.next().ok_or_else(|| format!("{flag} needs a value"))?;
        match flag.as_str() {
            "--dir" => options.dir = PathBuf::from(value),
            "--filter" => options.filters.extend(value.split(',').map(String::from)),
            "--field" => {
                options.fields = value
                    .split(',')
                    .map(|f| Field::parse(f).ok_or_else(|| format!("unknown field {f}")))
                    .collect::<Result<_, _>>()?;
            }
            "--algo" => {
                options.algorithms = value
                    .split(',')
                    .map(|a| Algorithm::parse(a).ok_or_else(|| format!("unknown algorithm {a}")))
                    .collect::<Result<_, _>>()?;
            }
            "--timeout" => {
                let secs: f64 = value.parse().map_err(|_| "bad --timeout".to_string())?;
                options.timeout = Duration::from_secs_f64(secs);
            }
            "--jobs" => options.jobs = value.parse().map_err(|_| "bad --jobs".to_string())?,
            "--out" => options.out = Some(PathBuf::from(value)),
            _ => return Err(format!("unknown flag {flag}")),
        }
    }
    Ok(options)
}

fn child(args: &[String]) -> ExitCode {
    let [path, field, algorithm] = args else {
        return ExitCode::from(2);
    };
    let (Some(field), Some(algorithm)) = (Field::parse(field), Algorithm::parse(algorithm)) else {
        return ExitCode::from(2);
    };
    let system = match corpus::load(Path::new(path)) {
        Ok(system) => system,
        Err(e) => {
            println!("error\t0\t0\t{e}");
            return ExitCode::from(1);
        }
    };
    let Some(reference) = &system.reference else {
        println!("error\t0\t0\tno reference");
        return ExitCode::from(1);
    };
    let start = Instant::now();
    let result = corpus::run(&system, field, algorithm);
    let seconds = start.elapsed().as_secs_f64();
    match result {
        Ok(rows) => match corpus::compare(&rows, &reference.rows) {
            Ok(()) => println!("ok\t{seconds}\t{}\t", rows.len()),
            Err(e) => println!("mismatch\t{seconds}\t{}\t{e}", rows.len()),
        },
        Err(e) => println!("error\t{seconds}\t0\t{e}"),
    }
    ExitCode::SUCCESS
}

fn parent(options: &Options) -> ExitCode {
    let systems: Vec<System> = corpus::load_dir(&options.dir)
        .into_iter()
        .filter(|s| {
            options.filters.is_empty()
                || options.filters.iter().any(|f| s.name.contains(f.as_str()))
        })
        .collect();
    let jobs: Vec<(&System, Field, Algorithm)> = systems
        .iter()
        .flat_map(|s| {
            options
                .fields
                .iter()
                .flat_map(move |&f| options.algorithms.iter().map(move |&a| (s, f, a)))
        })
        .filter(|(s, f, _)| f.applies(s))
        .collect();
    if jobs.is_empty() {
        eprintln!("no systems matched in {}", options.dir.display());
        return ExitCode::from(2);
    }

    let exe = match env::current_exe() {
        Ok(exe) => exe,
        Err(e) => {
            eprintln!("{e}");
            return ExitCode::from(2);
        }
    };
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let outcomes = Mutex::new(Vec::new());
    let started = Instant::now();
    thread::scope(|scope| {
        for _ in 0..options.jobs.max(1) {
            scope.spawn(|| {
                while let Some(&(system, field, algorithm)) = jobs.get(next.fetch_add(1, Ordering::SeqCst)) {
                    let outcome = execute(&exe, system, field, algorithm, options.timeout);
                    let count = done.fetch_add(1, Ordering::SeqCst) + 1;
                    println!(
                        "[{count:>4}/{}] {:<44} {:<8} {:<10} {:<8} {:>9.3}s  ref {:>8.3}s  {:>5} polys  {}",
                        jobs.len(),
                        outcome.name,
                        format!("{:?}", outcome.field).to_lowercase(),
                        format!("{:?}", outcome.algorithm).to_lowercase(),
                        outcome.status.label(),
                        outcome.seconds,
                        outcome.reference,
                        outcome.size,
                        outcome.detail,
                    );
                    if let Ok(mut all) = outcomes.lock() {
                        all.push(outcome);
                    }
                }
            });
        }
    });
    let outcomes = outcomes.into_inner().unwrap_or_default();
    summarize(&outcomes, started.elapsed(), options)
}

fn execute(
    exe: &Path,
    system: &System,
    field: Field,
    algorithm: Algorithm,
    timeout: Duration,
) -> Outcome {
    let mut outcome = Outcome {
        name: system.name.clone(),
        field,
        algorithm,
        status: Status::Crash,
        seconds: 0.0,
        reference: system.reference.as_ref().map_or(0.0, |r| r.time),
        size: 0,
        detail: String::new(),
    };
    let spawned = Command::new(exe)
        .arg("--one")
        .arg(&system.path)
        .arg(format!("{field:?}").to_lowercase())
        .arg(format!("{algorithm:?}").to_lowercase())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn();
    let mut process = match spawned {
        Ok(process) => process,
        Err(e) => {
            outcome.detail = e.to_string();
            return outcome;
        }
    };
    let start = Instant::now();
    let exit = loop {
        match process.try_wait() {
            Ok(Some(status)) => break Some(status),
            Ok(None) if start.elapsed() > timeout => {
                let _ = process.kill();
                let _ = process.wait();
                break None;
            }
            Ok(None) => thread::sleep(Duration::from_millis(10)),
            Err(_) => break None,
        }
    };
    outcome.seconds = start.elapsed().as_secs_f64();
    let Some(exit) = exit else {
        outcome.status = Status::Timeout;
        return outcome;
    };
    let mut stdout = String::new();
    let mut stderr = String::new();
    if let Some(mut pipe) = process.stdout.take() {
        let _ = pipe.read_to_string(&mut stdout);
    }
    if let Some(mut pipe) = process.stderr.take() {
        let _ = pipe.read_to_string(&mut stderr);
    }
    let fields: Vec<&str> = stdout
        .trim_end_matches(['\n', '\r'])
        .splitn(4, '\t')
        .collect();
    match (exit.success(), fields.as_slice()) {
        (true, [status, seconds, size, detail]) => {
            outcome.status = match *status {
                "ok" => Status::Ok,
                "mismatch" => Status::Mismatch,
                _ => Status::Error,
            };
            outcome.seconds = seconds.parse().unwrap_or(outcome.seconds);
            outcome.size = size.parse().unwrap_or(0);
            outcome.detail = (*detail).to_string();
        }
        _ => {
            outcome.detail = stderr
                .lines()
                .find(|l| l.contains("panicked") || l.contains("error"))
                .or_else(|| stderr.lines().next())
                .map_or_else(|| format!("exit {exit}"), str::to_string);
        }
    }
    outcome
}

fn summarize(outcomes: &[Outcome], elapsed: Duration, options: &Options) -> ExitCode {
    let count = |s: Status| outcomes.iter().filter(|o| o.status == s).count();
    println!(
        "\n{} runs in {:.1}s: {} ok, {} mismatch, {} error, {} crash, {} timeout",
        outcomes.len(),
        elapsed.as_secs_f64(),
        count(Status::Ok),
        count(Status::Mismatch),
        count(Status::Error),
        count(Status::Crash),
        count(Status::Timeout),
    );
    let mut failures: Vec<&Outcome> = outcomes
        .iter()
        .filter(|o| !matches!(o.status, Status::Ok | Status::Timeout))
        .collect();
    failures.sort_by(|a, b| a.name.cmp(&b.name));
    for o in &failures {
        println!(
            "  {:<8} {} {:?} {:?}: {}",
            o.status.label(),
            o.name,
            o.field,
            o.algorithm,
            o.detail
        );
    }
    if let Some(path) = &options.out {
        let mut sorted: Vec<&Outcome> = outcomes.iter().collect();
        sorted.sort_by(|a, b| a.name.cmp(&b.name));
        let mut tsv =
            String::from("name\tfield\talgorithm\tstatus\tseconds\treference\tsize\tdetail\n");
        for o in sorted {
            tsv.push_str(&format!(
                "{}\t{:?}\t{:?}\t{}\t{:.4}\t{:.4}\t{}\t{}\n",
                o.name,
                o.field,
                o.algorithm,
                o.status.label(),
                o.seconds,
                o.reference,
                o.size,
                o.detail
            ));
        }
        if let Err(e) = fs::write(path, tsv) {
            eprintln!("could not write {}: {e}", path.display());
        }
    }
    if failures.is_empty() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
