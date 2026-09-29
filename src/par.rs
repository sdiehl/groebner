//! Serial and Rayon versions of the parallel patterns the crate uses, chosen by the
//! `parallel` feature, so callers need no feature gates of their own.

#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cmp::Ordering;

pub(crate) fn threads() -> usize {
    #[cfg(feature = "parallel")]
    return rayon::current_num_threads();
    #[cfg(not(feature = "parallel"))]
    1
}

// Runs `f` on a pool thread, so parallel steps inside it split work by stealing instead
// of each waking the pool from outside and sleeping until it finishes.
pub(crate) fn on_pool<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    #[cfg(feature = "parallel")]
    if rayon::current_thread_index().is_none() {
        return rayon::scope(|_| f());
    }
    f()
}

pub(crate) fn join<A: Send, B: Send>(
    a: impl FnOnce() -> A + Send,
    b: impl FnOnce() -> B + Send,
) -> (A, B) {
    #[cfg(feature = "parallel")]
    return rayon::join(a, b);
    #[cfg(not(feature = "parallel"))]
    (a(), b())
}

pub(crate) fn map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    #[cfg(feature = "parallel")]
    return items.par_iter().map(&f).collect();
    #[cfg(not(feature = "parallel"))]
    items.iter().map(f).collect()
}

pub(crate) fn sort_unstable_by<T: Send>(items: &mut [T], cmp: impl Fn(&T, &T) -> Ordering + Sync) {
    #[cfg(feature = "parallel")]
    items.par_sort_unstable_by(cmp);
    #[cfg(not(feature = "parallel"))]
    items.sort_unstable_by(cmp);
}

// The outputs of `f` on consecutive chunks of `size` items, concatenated.
pub(crate) fn flat_map_chunks<T: Sync, R: Send, I: IntoIterator<Item = R>>(
    items: &[T],
    size: usize,
    f: impl Fn(&[T]) -> I + Sync,
) -> Vec<R> {
    #[cfg(feature = "parallel")]
    return items.par_chunks(size).flat_map_iter(&f).collect();
    #[cfg(not(feature = "parallel"))]
    items.chunks(size).flat_map(f).collect()
}

// `f` on each item and its index, each worker reusing one scratch value.
pub(crate) fn for_each_init<T: Sync, S>(
    items: &[T],
    scratch: impl Fn() -> S + Sync + Send,
    f: impl Fn(&mut S, usize, &T) + Sync + Send,
) {
    #[cfg(feature = "parallel")]
    items
        .par_iter()
        .enumerate()
        .for_each_init(scratch, |s, (i, item)| f(s, i, item));
    #[cfg(not(feature = "parallel"))]
    {
        let mut s = scratch();
        items
            .iter()
            .enumerate()
            .for_each(|(i, item)| f(&mut s, i, item));
    }
}

// Map rows in parallel blocks, each worker reusing one scratch value.
pub(crate) fn map_rows<T: Sync, R: Send, S>(
    rows: &[T],
    scratch: impl Fn() -> S + Sync,
    f: impl Fn(&T, &mut S) -> R + Sync,
) -> Vec<R> {
    map_blocks(rows, 32, scratch, f)
}

// `map_rows` in runs of about equal total `weight`, as rows of similar cost cluster
// and equal counts would leave one worker a heavy run.
pub(crate) fn map_weighted<T: Sync, R: Send, S>(
    rows: &[T],
    weight: impl Fn(&T) -> usize,
    scratch: impl Fn() -> S + Sync,
    f: impl Fn(&T, &mut S) -> R + Sync,
) -> Vec<R> {
    #[cfg(feature = "parallel")]
    if rows.len() >= 64 {
        let total: usize = rows.iter().map(&weight).sum();
        let target = total.div_ceil(16 * threads()).max(1);
        let (mut runs, mut start, mut sum) = (Vec::new(), 0, 0);
        for (i, r) in rows.iter().enumerate() {
            sum += weight(r);
            if sum >= target {
                runs.push(&rows[start..=i]);
                (start, sum) = (i + 1, 0);
            }
        }
        runs.push(&rows[start..]);
        return runs
            .par_iter()
            .with_max_len(1)
            .flat_map_iter(|run| {
                let mut s = scratch();
                run.iter().map(|r| f(r, &mut s)).collect::<Vec<_>>()
            })
            .collect();
    }
    #[cfg(not(feature = "parallel"))]
    let _ = weight;
    map_rows(rows, scratch, f)
}

// `map_rows` in blocks of `size` rows, split only when parallel and there are at
// least two.
pub(crate) fn map_blocks<T: Sync, R: Send, S>(
    rows: &[T],
    size: usize,
    scratch: impl Fn() -> S + Sync,
    f: impl Fn(&T, &mut S) -> R + Sync,
) -> Vec<R> {
    let block = |block: &[T]| {
        let mut s = scratch();
        block.iter().map(|r| f(r, &mut s)).collect::<Vec<_>>()
    };
    if cfg!(feature = "parallel") && rows.len() >= 2 * size {
        return flat_map_chunks(rows, size, block);
    }
    block(rows)
}

pub(crate) fn for_rows_mut<T: Send>(rows: &mut [T], f: impl Fn(&mut T) + Sync) {
    #[cfg(feature = "parallel")]
    if rows.len() >= 64 {
        return rows
            .par_chunks_mut(32)
            .for_each(|block| block.iter_mut().for_each(&f));
    }
    rows.iter_mut().for_each(f);
}

// `map_rows` for cheap per-row work, which only pays to split when there is a lot.
pub(crate) fn map_big<T: Sync, R: Send>(rows: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    if rows.len() < 1024 {
        return rows.iter().map(f).collect();
    }
    map_rows(rows, || (), |r, ()| f(r))
}
