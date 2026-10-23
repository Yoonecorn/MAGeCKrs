/// Opt-in scheduler policy. A fixed limit always bypasses load-based selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThreadMode {
    Fixed(usize),
    Auto { maximum: Option<usize> },
}
impl ThreadMode {
    pub fn is_auto(self) -> bool {
        matches!(self, Self::Auto { .. })
    }
}
impl std::str::FromStr for ThreadMode {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value == "auto" {
            return Ok(Self::Auto { maximum: None });
        }
        if let Some(maximum) = value.strip_prefix("auto:") {
            return maximum
                .parse::<usize>()
                .ok()
                .filter(|n| *n > 0)
                .map(|n| Self::Auto { maximum: Some(n) })
                .ok_or_else(|| "--threads auto:MAX requires a positive integer maximum".into());
        }
        value
            .parse::<usize>()
            .map(Self::Fixed)
            .map_err(|_| "--threads must be auto, auto:MAX, or an integer worker limit".into())
    }
}

// Opt-in bounded stage timings; ordinary execution emits no profiling metadata.
use std::time::Instant;

pub(crate) struct Profile {
    start: Instant,
    previous: Instant,
    enabled: bool,
    workflow: &'static str,
    events: usize,
}

impl Profile {
    pub(crate) fn new(workflow: &'static str) -> Self {
        let start = Instant::now();
        Self {
            start,
            previous: start,
            enabled: std::env::var_os("MAGECKRS_PROFILE").is_some(),
            workflow,
            events: 0,
        }
    }
    pub(crate) fn stage(&mut self, stage: &str) {
        if self.enabled {
            if self.events >= 64 {
                return;
            }
            self.events += 1;
            let now = Instant::now();
            eprintln!(
                "PROFILE\t{}\t{}\t{:.9}",
                self.workflow,
                stage,
                now.duration_since(self.previous).as_secs_f64()
            );
            self.previous = now;
        }
    }
    pub(crate) fn finish(&self) {
        if self.enabled {
            eprintln!(
                "PROFILE\t{}\ttotal\t{:.9}",
                self.workflow,
                self.start.elapsed().as_secs_f64()
            );
        }
    }
}

use crate::error::MageckError;
use rayon::prelude::*;

// Native macOS APIs from SDK mach/host_info.h and sys/sysctl.h. No affinity,
// core-class queries, subprocesses, dependencies or global pool configuration.
unsafe extern "C" {
    fn mach_host_self() -> u32;
    static mach_task_self_: u32;
    fn mach_port_deallocate(task: u32, port: u32) -> i32;
    fn host_statistics(host: u32, flavor: i32, info: *mut i32, count: *mut u32) -> i32;
    fn sysctlbyname(
        name: *const std::ffi::c_char,
        old: *mut std::ffi::c_void,
        size: *mut usize,
        new: *mut std::ffi::c_void,
        new_size: usize,
    ) -> i32;
}

fn sysctl(name: &std::ffi::CStr, buffer: &mut [u8]) -> Option<usize> {
    let mut size = buffer.len();
    // SAFETY: buffer is writable for size bytes; name is NUL terminated. Read only.
    let status = unsafe {
        sysctlbyname(
            name.as_ptr(),
            buffer.as_mut_ptr().cast(),
            &mut size,
            std::ptr::null_mut(),
            0,
        )
    };
    (status == 0 && size <= buffer.len()).then_some(size)
}
fn sys_u64(name: &std::ffi::CStr) -> Option<u64> {
    let mut bytes = [0; 8];
    (sysctl(name, &mut bytes)? == 8).then(|| u64::from_ne_bytes(bytes))
}
fn sys_u32(name: &std::ffi::CStr) -> Option<usize> {
    let mut bytes = [0; 4];
    (sysctl(name, &mut bytes)? == 4).then(|| u32::from_ne_bytes(bytes) as usize)
}
fn sys_text(name: &std::ffi::CStr) -> String {
    let mut bytes = [0; 256];
    sysctl(name, &mut bytes).map_or_else(
        || "unavailable".into(),
        |size| {
            String::from_utf8_lossy(&bytes[..size])
                .trim_end_matches('\0')
                .replace(['\t', '\r', '\n'], " ")
        },
    )
}
fn ticks() -> Result<[u32; 4], &'static str> {
    let mut ticks = [0i32; 4];
    let mut count = 4;
    // HOST_CPU_LOAD_INFO=3, CPU_STATE_IDLE=2, count=4 (32-bit tick counters).
    // SAFETY: exact SDK layout and count; host send right is released every call.
    let status = unsafe {
        let host = mach_host_self();
        let status = host_statistics(host, 3, ticks.as_mut_ptr(), &mut count);
        mach_port_deallocate(mach_task_self_, host);
        status
    };
    if status != 0 || count != 4 {
        return Err("host_statistics_unavailable");
    }
    Ok(ticks.map(|tick| tick as u32))
}
fn busy_basis_points(before: [u32; 4], after: [u32; 4]) -> Result<u32, &'static str> {
    let delta = std::array::from_fn::<_, 4, _>(|i| after[i].wrapping_sub(before[i]) as u64);
    let total: u64 = delta.iter().sum();
    if total == 0 {
        return Err("no_cpu_tick_progress");
    }
    Ok(((total - delta[2]) * 10000 / total) as u32)
}
fn sample_load() -> Result<(u32, u128, usize), &'static str> {
    let start = Instant::now();
    let mut before = ticks()?;
    let mut busy = 0;
    let mut valid_windows = 0;
    // Two 250ms windows; use the busier window with tick progress. macOS
    // can cache snapshots; an unchanged window is skipped, never treated as idle.
    for _ in 0..2 {
        std::thread::sleep(std::time::Duration::from_millis(250));
        let after = ticks()?;
        if let Ok(sample) = busy_basis_points(before, after) {
            busy = busy.max(sample);
            valid_windows += 1;
        }
        before = after;
    }
    let elapsed = start.elapsed().as_millis();
    if elapsed > 1000 {
        return Err("sample_delayed_over_1000ms");
    }
    if valid_windows == 0 {
        return Err("no_cpu_tick_progress");
    }
    Ok((busy, elapsed, valid_windows))
}

const STACK_BYTES: usize = 8 * 1024 * 1024;
const MAX_INCREMENTAL_BYTES: usize = 256 * 1024 * 1024;
const MAX_AUTO_WORKERS: usize = 8;

/// Conservative incremental fit allowance; includes IRLS's n*n hat matrix.
/// Saturation means one serial fit, never an overflowing memory worker budget.
pub(crate) fn mle_worker_bytes(samples: usize, guides: usize, conditions: usize) -> usize {
    let n = samples
        .saturating_mul(2)
        .saturating_sub(1)
        .saturating_mul(guides);
    let p = guides.saturating_add(conditions);
    n.saturating_mul(n)
        .saturating_mul(2)
        .saturating_add(n.saturating_mul(p).saturating_mul(4))
        .saturating_add(p.saturating_mul(p).saturating_mul(16))
        .saturating_add(n.saturating_mul(32))
        .saturating_add(p.saturating_mul(16))
        .saturating_mul(8)
        .saturating_add(STACK_BYTES)
}

#[derive(Clone, Debug)]
struct Observation {
    capacity: usize,
    busy: Option<u32>,
    elapsed_ms: u128,
    valid_windows: usize,
    memory_budget: usize,
    fallback: String,
}
impl Observation {
    fn detect(auto: bool) -> Self {
        let capacity = std::thread::available_parallelism().ok().map(usize::from);
        let active = sys_u32(
            std::ffi::CStr::from_bytes_with_nul(b"hw.activecpu\0").expect("static sysctl name"),
        )
        .filter(|n| *n > 0);
        let memory = sys_u64(
            std::ffi::CStr::from_bytes_with_nul(b"hw.memsize\0").expect("static sysctl name"),
        );
        let mut fallback = Vec::new();
        if capacity.is_none() {
            fallback.push("available_parallelism_unavailable");
        }
        if active.is_none() {
            fallback.push("activecpu_unavailable");
        }
        if memory.is_none() {
            fallback.push("memory_capacity_unavailable");
        }
        let (busy, elapsed_ms, valid_windows) = if auto {
            match sample_load() {
                Ok((busy, elapsed, windows)) => (Some(busy), elapsed, windows),
                Err(reason) => {
                    fallback.push(reason);
                    (None, 0, 0)
                }
            }
        } else {
            (None, 0, 0)
        };
        Self {
            capacity: capacity.unwrap_or(1).min(active.unwrap_or(1)),
            busy,
            elapsed_ms,
            valid_windows,
            memory_budget: memory.map_or(0, |bytes| {
                (bytes / 32).min(MAX_INCREMENTAL_BYTES as u64) as usize
            }),
            fallback: if fallback.is_empty() {
                "none".into()
            } else {
                fallback.join(",")
            },
        }
    }
}

fn select(
    mode: Option<ThreadMode>,
    jobs: usize,
    scalability: usize,
    per_worker: usize,
    observation: &Observation,
) -> usize {
    let budget = match mode {
        Some(ThreadMode::Auto { maximum }) => {
            if observation.fallback != "none" {
                return 1;
            }
            // Floor 75% of observed idle capacity; reserve a logical CPU for
            // the serial controller/system. Always allow a scalar forward path.
            let idle = 10000 - observation.busy.unwrap_or(10000).min(10000);
            let idle_workers = observation
                .capacity
                .saturating_mul(idle as usize)
                .saturating_mul(3)
                / 40000;
            idle_workers
                .min(observation.capacity.saturating_sub(1).max(1))
                .min(MAX_AUTO_WORKERS)
                .min(maximum.unwrap_or(usize::MAX))
        }
        Some(ThreadMode::Fixed(n)) => n,
        None => 1,
    };
    budget
        .min(observation.capacity)
        .min(jobs.max(1))
        .min(scalability.max(1))
        .min(observation.memory_budget / per_worker.max(1))
        .max(1)
}

/// Workflow-local pools only. Auto MLE maps create/drop one pool per major stage.
/// Fixed mode bypasses load observations. No BLAS or environment settings change.
pub(crate) struct Execution {
    pool: Option<rayon::ThreadPool>,
    auto: Option<ThreadMode>,
    records: std::sync::Mutex<Vec<String>>,
    threads: Vec<std::thread::JoinHandle<()>>,
}
impl Drop for Execution {
    fn drop(&mut self) {
        // Rayon pool drop signals shutdown; join actual OS threads before a
        // subsequent stage can create workers, preventing transient overlap.
        drop(self.pool.take());
        for thread in self.threads.drain(..) {
            let _ = thread.join();
        }
    }
}
impl Execution {
    pub(crate) fn serial() -> Self {
        Self {
            pool: None,
            auto: None,
            records: std::sync::Mutex::new(Vec::new()),
            threads: Vec::new(),
        }
    }
    pub(crate) fn new(threads: Option<usize>, jobs: usize) -> Result<Self, MageckError> {
        if threads == Some(0) {
            return Err(MageckError::Config(
                "--threads must be greater than zero".into(),
            ));
        }
        let cpus = std::thread::available_parallelism().map_or(1, usize::from);
        Self::pool(threads.unwrap_or(1).min(jobs.max(1)).min(cpus), None)
    }
    fn pool(workers: usize, stack_bytes: Option<usize>) -> Result<Self, MageckError> {
        let mut result = Self::serial();
        if workers > 1 && rayon::current_thread_index().is_none() {
            let handles = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
            let spawned = handles.clone();
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .spawn_handler(move |thread| {
                    let builder = std::thread::Builder::new()
                        .name(format!("mageck-worker-{}", thread.index()));
                    // Fixed workers retain the platform's default stack setting.
                    let builder = match stack_bytes {
                        Some(bytes) => builder.stack_size(bytes),
                        None => builder,
                    };
                    let handle = builder.spawn(move || thread.run())?;
                    spawned.lock().expect("unpoisoned handles").push(handle);
                    Ok(())
                })
                .build();
            result.threads = std::mem::take(&mut *handles.lock().expect("unpoisoned handles"));
            result.pool = Some(pool.map_err(|error| {
                MageckError::Config(format!("Unable to create CPU worker pool: {error}"))
            })?);
        }
        Ok(result)
    }
    pub(crate) fn workflow(mode: Option<ThreadMode>, jobs: usize) -> Result<Self, MageckError> {
        match mode {
            Some(mode @ ThreadMode::Auto { .. }) => {
                let mut execution = Self::serial();
                execution.auto = Some(mode);
                Ok(execution)
            }
            Some(ThreadMode::Fixed(n)) => Self::new(Some(n), jobs),
            None => Self::new(None, jobs),
        }
    }
    pub(crate) fn stage(
        mode: Option<ThreadMode>,
        jobs: usize,
        scalability: usize,
        per_worker: usize,
        stage: &str,
    ) -> Result<Self, MageckError> {
        if !mode.is_some_and(ThreadMode::is_auto) {
            return Self::workflow(mode, jobs);
        }
        Self::stage_observed(
            mode,
            jobs,
            scalability,
            per_worker,
            stage,
            Observation::detect(true),
        )
    }
    fn stage_observed(
        mode: Option<ThreadMode>,
        jobs: usize,
        scalability: usize,
        per_worker: usize,
        stage: &str,
        mut observation: Observation,
    ) -> Result<Self, MageckError> {
        let mut workers = select(mode, jobs, scalability, per_worker, &observation);
        if rayon::current_thread_index().is_some() {
            workers = 1;
            observation.fallback = "nested_rayon_serial".into();
        }
        let mut execution = match Self::pool(workers, Some(STACK_BYTES)) {
            Ok(execution) => execution,
            Err(_) => {
                workers = 1;
                observation.fallback = "worker_pool_unavailable".into();
                Self::serial()
            }
        };
        execution.records.get_mut().expect("unpoisoned metadata").push(format!(
            "threads.{stage}.policy\tauto-v1\nthreads.{stage}.capacity\t{}\nthreads.{stage}.busy_basis_points\t{}\nthreads.{stage}.sample_ms\t{}\nthreads.{stage}.valid_windows\t{}\nthreads.{stage}.workers\t{workers}\nthreads.{stage}.jobs\t{jobs}\nthreads.{stage}.scalability_cap\t{scalability}\nthreads.{stage}.maximum\t{:?}\nthreads.{stage}.worker_bytes\t{per_worker}\nthreads.{stage}.memory_budget_bytes\t{}\nthreads.{stage}.fallback\t{}\nthreads.{stage}.host_model\t{}\nthreads.{stage}.os_release\t{}\nthreads.{stage}.hostname\t{}\nthreads.{stage}.architecture\taarch64\nthreads.{stage}.placement\tmacos_scheduler\nthreads.{stage}.sampler\thost_statistics_two_250ms_windows_max_busy_with_tick_progress\n",
            observation.capacity, observation.busy.map_or_else(|| "unavailable".into(), |n| n.to_string()),
            observation.elapsed_ms, observation.valid_windows, mode, observation.memory_budget, observation.fallback,
            sys_text(std::ffi::CStr::from_bytes_with_nul(b"hw.model\0").expect("static sysctl name")), sys_text(std::ffi::CStr::from_bytes_with_nul(b"kern.osrelease\0").expect("static sysctl name")),
            sys_text(std::ffi::CStr::from_bytes_with_nul(b"kern.hostname\0").expect("static sysctl name"))));
        Ok(execution)
    }
    pub(crate) fn metadata(&self) -> String {
        self.records.lock().expect("unpoisoned metadata").concat()
    }
    pub(crate) fn is_parallel(&self) -> bool {
        self.pool.is_some()
    }
    pub(crate) fn map_stage<T: Send, U: Send, F: Fn(T) -> Result<U, MageckError> + Sync + Send>(
        &self,
        stage: &str,
        per_worker: usize,
        items: Vec<T>,
        operation: F,
    ) -> Result<Vec<U>, MageckError> {
        if let Some(mode) = self.auto {
            let execution =
                Self::stage(Some(mode), items.len(), MAX_AUTO_WORKERS, per_worker, stage)?;
            self.records
                .lock()
                .expect("unpoisoned metadata")
                .push(execution.metadata());
            // Local pool is dropped before the next stage can sample or create one.
            execution.map(items, operation)
        } else {
            self.map(items, operation)
        }
    }
    /// Collect every outcome in indexed source order, then select the first error.
    /// Rayon short-circuiting would select errors according to scheduling.
    pub(crate) fn map<T: Send, U: Send, F: Fn(T) -> Result<U, MageckError> + Sync + Send>(
        &self,
        items: Vec<T>,
        operation: F,
    ) -> Result<Vec<U>, MageckError> {
        match &self.pool {
            Some(pool) => {
                let outcomes: Vec<_> =
                    pool.install(|| items.into_par_iter().map(operation).collect());
                outcomes.into_iter().collect()
            }
            None => items.into_iter().map(operation).collect(),
        }
    }
}

#[cfg(test)]
#[allow(dead_code)]
impl Execution {
    pub(crate) fn corrected_probe_budget(jobs: usize, busy: Option<u32>, reason: &str) -> Self {
        Self::stage_observed(
            Some(ThreadMode::Auto { maximum: Some(8) }),
            jobs,
            8,
            STACK_BYTES,
            "test.corrected",
            Observation {
                capacity: 16,
                busy,
                elapsed_ms: 500,
                valid_windows: usize::from(busy.is_some()),
                memory_budget: MAX_INCREMENTAL_BYTES,
                fallback: reason.into(),
            },
        )
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn observation(capacity: usize, busy: Option<u32>, fallback: &str) -> Observation {
        Observation {
            capacity,
            busy,
            elapsed_ms: 500,
            valid_windows: usize::from(busy.is_some()),
            memory_budget: MAX_INCREMENTAL_BYTES,
            fallback: fallback.into(),
        }
    }
    #[test]
    fn native_stage_and_workflow_metadata_are_available() {
        let auto = Some(ThreadMode::Auto { maximum: Some(2) });
        let stage = Execution::stage(auto, 1, 1, STACK_BYTES, "test.native").unwrap();
        assert!(stage.metadata().contains("test.native.capacity"));
        assert!(stage.metadata().contains("workers\t1\n"));
        let workflow = Execution::workflow(auto, 4).unwrap();
        assert_eq!(
            workflow
                .map_stage("test.workflow", STACK_BYTES, (0..4).collect(), Ok)
                .unwrap(),
            (0..4).collect::<Vec<_>>()
        );
        assert!(workflow
            .metadata()
            .contains("test.workflow.policy\tauto-v1"));
        assert!(Execution::workflow(Some(ThreadMode::Fixed(0)), 4).is_err());
        assert!(Execution::workflow(None, 4).unwrap().metadata().is_empty());
    }
    #[test]
    fn auto_policy_is_conservative_monotone_and_bounded() {
        let auto = Some(ThreadMode::Auto { maximum: None });
        for capacity in [1, 2, 4, 8, 16, 64] {
            let mut previous = usize::MAX;
            for busy in (0..=10000).step_by(100) {
                let obs = observation(capacity, Some(busy), "none");
                let workers = select(auto, 100, 64, STACK_BYTES, &obs);
                assert!(workers >= 1 && workers <= capacity && workers <= MAX_AUTO_WORKERS);
                assert!(workers <= previous);
                previous = workers;
                assert_eq!(workers, select(auto, 100, 64, STACK_BYTES, &obs));
            }
        }
        let idle = observation(16, Some(0), "none");
        assert_eq!(select(auto, 100, 64, STACK_BYTES, &idle), 8);
        assert_eq!(select(auto, 2, 64, STACK_BYTES, &idle), 2);
        assert_eq!(select(auto, 100, 3, STACK_BYTES, &idle), 3);
        assert_eq!(
            select(
                Some(ThreadMode::Auto { maximum: Some(2) }),
                100,
                64,
                STACK_BYTES,
                &idle
            ),
            2
        );
        assert_eq!(select(auto, 100, 64, MAX_INCREMENTAL_BYTES / 2, &idle), 2);
        assert_eq!(select(auto, 100, 64, usize::MAX, &idle), 1);
        let loaded = observation(16, Some(10000), "none");
        assert_eq!(select(auto, 100, 64, STACK_BYTES, &loaded), 1);
        assert_eq!(
            select(Some(ThreadMode::Fixed(6)), 100, 64, STACK_BYTES, &loaded),
            6
        );
        assert_eq!(
            select(
                Some(ThreadMode::Fixed(usize::MAX)),
                100,
                64,
                STACK_BYTES,
                &idle
            ),
            16
        );
        assert_eq!(select(None, 100, 64, STACK_BYTES, &idle), 1);
    }
    #[test]
    fn unavailable_detection_records_serial_fallback_and_correct_errors() {
        let auto = Some(ThreadMode::Auto { maximum: None });
        for reason in [
            "host_statistics_unavailable",
            "no_cpu_tick_progress",
            "sample_delayed_over_1000ms",
            "available_parallelism_unavailable",
            "memory_capacity_unavailable",
            "activecpu_unavailable",
        ] {
            let execution = Execution::stage_observed(
                auto,
                20,
                8,
                STACK_BYTES,
                "test.fallback",
                observation(16, None, reason),
            )
            .unwrap();
            assert!(!execution.is_parallel());
            assert!(execution.metadata().contains(reason));
            assert!(execution.metadata().contains("workers\t1\n"));
            assert_eq!(
                execution.map((0..20).collect(), Ok).unwrap(),
                (0..20).collect::<Vec<_>>()
            );
            assert_eq!(
                execution
                    .map((0..5).collect(), |i| if i >= 2 {
                        Err(MageckError::Math(format!("source {i}")))
                    } else {
                        Ok(i)
                    })
                    .unwrap_err()
                    .to_string(),
                "Math/Modeling error: source 2"
            );
        }
    }
    #[test]
    fn native_tick_math_handles_wrap_and_no_progress() {
        assert_eq!(busy_basis_points([0; 4], [20, 10, 70, 0]).unwrap(), 3000);
        assert!(busy_basis_points([1; 4], [1; 4]).is_err());
        assert_eq!(
            busy_basis_points([u32::MAX - 9, 0, 0, 0], [0, 0, 10, 0]).unwrap(),
            5000
        );
        assert_eq!(mle_worker_bytes(usize::MAX, usize::MAX, 1), usize::MAX);
        assert!(mle_worker_bytes(3, 100, 2) > mle_worker_bytes(3, 10, 2));
    }
    #[test]
    fn live_workers_are_bounded_and_nested_execution_stays_serial() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let active = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let auto = Some(ThreadMode::Auto { maximum: Some(3) });
        let execution = Execution::stage_observed(
            auto,
            30,
            8,
            STACK_BYTES,
            "test.workers",
            observation(16, Some(0), "none"),
        )
        .unwrap();
        let result = execution
            .map((0..30).collect(), |i| {
                let now = active.fetch_add(1, Ordering::SeqCst) + 1;
                peak.fetch_max(now, Ordering::SeqCst);
                let nested = Execution::stage_observed(
                    auto,
                    20,
                    8,
                    STACK_BYTES,
                    "test.nested",
                    observation(16, Some(0), "none"),
                )?;
                assert!(!nested.is_parallel());
                assert!(nested.metadata().contains("nested_rayon_serial"));
                std::thread::sleep(std::time::Duration::from_millis(2));
                active.fetch_sub(1, Ordering::SeqCst);
                Ok(i)
            })
            .unwrap();
        assert_eq!(result, (0..30).collect::<Vec<_>>());
        assert!((2..=3).contains(&peak.load(Ordering::SeqCst)));
        assert_eq!(active.load(Ordering::SeqCst), 0);
        for _ in 0..3 {
            let error = execution
                .map((0..20).collect(), |i| {
                    if i == 2 {
                        std::thread::sleep(std::time::Duration::from_millis(5));
                    }
                    if i == 2 || i == 3 {
                        Err(MageckError::Math(format!("source {i}")))
                    } else {
                        Ok(i)
                    }
                })
                .unwrap_err();
            assert_eq!(error.to_string(), "Math/Modeling error: source 2");
        }
        assert_eq!(execution.threads.len(), 3);
        drop(execution); // all handles joined, no worker survives this boundary
    }
    #[test]
    fn huge_positive_limit_is_bounded_and_zero_is_invalid() {
        assert!(Execution::new(Some(0), 1).is_err());
        let execution = Execution::new(Some(usize::MAX), 32).unwrap();
        assert_eq!(
            execution.map((0..32).collect(), Ok).unwrap(),
            (0..32).collect::<Vec<_>>()
        );
    }
    #[test]
    fn scheduling_never_selects_a_later_error() {
        for threads in [None, Some(1), Some(2), Some(4), Some(8)] {
            let execution = Execution::new(threads, 20).unwrap();
            for _ in 0..3 {
                let error = execution
                    .map((0..20).collect(), |i| {
                        if i == 2 {
                            std::thread::sleep(std::time::Duration::from_millis(5));
                        }
                        if i == 2 || i == 3 {
                            Err(MageckError::Math(format!("source {i}")))
                        } else {
                            Ok(i)
                        }
                    })
                    .unwrap_err();
                assert_eq!(error.to_string(), "Math/Modeling error: source 2");
            }
        }
    }
}
