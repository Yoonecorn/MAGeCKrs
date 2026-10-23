/// build.rs — MAGeCKrs
///
/// Responsibilities:
///   1. Enforce macOS-only target (panic on Linux/Windows).
///   2. Link Apple Accelerate framework; emit `has_accelerate` and `apple_silicon` cfg flags.
///   3. If `metal` feature is active: link Metal + Foundation; emit `has_metal`.
///   4. If `metal` feature is active: attempt to precompile the Metal shader to .metallib;
///      emit `has_metallib` and `MAGECK_METALLIB_PATH` if successful.
///      Fall back to JIT source compilation if xcrun is unavailable.
///
/// Agents: do not modify the platform check or Accelerate linking without explicit instruction.
/// The has_accelerate and apple_silicon cfg flags must always be emitted on macOS aarch64.
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();

    // ── Platform gate ──────────────────────────────────────────────────────────
    if target_os != "macos" {
        panic!(
            "mageckrs targets macOS only.\n\
             Linux and Windows are not supported build targets.\n\
             If you are on a non-macOS machine and only need to check syntax,\n\
             this build script still enforces the platform requirement by design."
        );
    }

    // ── Accelerate (always active on macOS) ────────────────────────────────────
    // Links the Accelerate umbrella framework, which includes vDSP, vForce, and
    // BLAS/LAPACK (routed through the AMX coprocessor on Apple Silicon).
    println!("cargo:rustc-link-lib=framework=Accelerate");
    println!("cargo:rustc-cfg=has_accelerate");

    if target_arch == "aarch64" {
        // Enables Apple Silicon-specific code paths (AMX verification, P-core tuning).
        println!("cargo:rustc-cfg=apple_silicon");
    }

    // ── Metal (optional feature) ───────────────────────────────────────────────
    let metal_feature = std::env::var("CARGO_FEATURE_METAL").is_ok();
    if metal_feature {
        println!("cargo:rustc-link-lib=framework=Metal");
        println!("cargo:rustc-link-lib=framework=Foundation");
        println!("cargo:rustc-cfg=has_metal");

        compile_metal_shader();
    }
}

/// Attempts to AOT-compile shaders/permutation_fdr.metal to a .metallib binary.
/// Sets `has_metallib` cfg and `MAGECK_METALLIB_PATH` env var if successful.
/// Falls back to JIT source compilation (via include_str! in metal_gpu.rs) if xcrun
/// is unavailable or compilation fails.
fn compile_metal_shader() {
    let shader_src = PathBuf::from("shaders/permutation_fdr.metal");

    // Rerun this script if the shader source changes.
    println!("cargo:rerun-if-changed=shaders/permutation_fdr.metal");

    // If the shader file doesn't exist yet (Phase 8 not started), skip AOT compilation.
    // The metal_gpu.rs module will use JIT compilation until the shader file is created.
    if !shader_src.exists() {
        println!(
            "cargo:warning=shaders/permutation_fdr.metal not found; \
             Metal JIT compilation will be used at runtime."
        );
        return;
    }

    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let air_path = out_dir.join("permutation_fdr.air");
    let metallib_path = out_dir.join("permutation_fdr.metallib");

    // Step 1: Compile .metal → .air
    let metal_result = Command::new("xcrun")
        .args([
            "-sdk",
            "macosx",
            "metal",
            "-o",
            air_path.to_str().unwrap(),
            shader_src.to_str().unwrap(),
        ])
        .status();

    match metal_result {
        Ok(status) if status.success() => {}
        Ok(status) => {
            println!(
                "cargo:warning=xcrun metal exited with status {status}; \
                 falling back to JIT shader compilation."
            );
            return;
        }
        Err(e) => {
            println!(
                "cargo:warning=xcrun metal not available ({e}); \
                 falling back to JIT shader compilation. \
                 Install Xcode Command Line Tools to enable AOT compilation."
            );
            return;
        }
    }

    // Step 2: Link .air → .metallib
    let metallib_result = Command::new("xcrun")
        .args([
            "-sdk",
            "macosx",
            "metallib",
            "-o",
            metallib_path.to_str().unwrap(),
            air_path.to_str().unwrap(),
        ])
        .status();

    match metallib_result {
        Ok(status) if status.success() => {
            // AOT compilation succeeded. Embed the .metallib path so metal_gpu.rs
            // can use include_bytes!(env!("MAGECK_METALLIB_PATH")) for a self-contained binary.
            println!(
                "cargo:rustc-env=MAGECK_METALLIB_PATH={}",
                metallib_path.display()
            );
            println!("cargo:rustc-cfg=has_metallib");
        }
        Ok(status) => {
            println!(
                "cargo:warning=xcrun metallib exited with status {status}; \
                 falling back to JIT shader compilation."
            );
        }
        Err(e) => {
            println!(
                "cargo:warning=xcrun metallib failed ({e}); \
                 falling back to JIT shader compilation."
            );
        }
    }
}
