//! Integration tests for AC analysis, porting ngspice-upstream/tests/filters/ tests.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use ferrospice_core::simulate_ac;
use ferrospice_netlist::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Port of ngspice-upstream/tests/filters/lowpass.cir
///
/// Circuit: V1 (AC=1) -> R1(10k) -> node 2, R2(10k) from node 2 to ground,
///          C1(1n) from node 2 to ground.
/// Sweep: .ac DEC 10 1k 1Meg → 31 frequency points.
#[test]
fn test_lowpass_filter() {
    let cir_path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../ngspice-upstream/tests/filters/lowpass.cir"
    );
    let cir_content = std::fs::read_to_string(cir_path).expect("failed to read lowpass.cir");
    let netlist = Netlist::parse(&cir_content).expect("failed to parse lowpass.cir");

    let result = simulate_ac(&netlist).expect("AC simulation failed");
    assert_eq!(result.plots.len(), 1);

    let plot = &result.plots[0];
    let freq_vec = plot.vecs.iter().find(|v| v.name == "frequency").unwrap();
    let v2_vec = plot.vecs.iter().find(|v| v.name == "v(2)").unwrap();

    // Should have 31 frequency points (3 decades × 10 pts/decade + 1).
    assert_eq!(freq_vec.real.len(), 31);
    assert_eq!(v2_vec.complex.len(), 31);

    // Reference values from ngspice lowpass.out (selected points).
    // Format: (index, frequency, v2_real, v2_imag)
    let reference: &[(usize, f64, f64, f64)] = &[
        (0, 1.000000e+03, 4.995070e-01, -1.569248e-02),
        (10, 1.000000e+04, 4.550849e-01, -1.429691e-01),
        (15, 3.162278e+04, 2.516406e-01, -2.499946e-01),
        (20, 1.000000e+05, 4.599983e-02, -1.445127e-01),
        (30, 1.000000e+06, 5.060931e-04, -1.589938e-02),
    ];

    for &(idx, ref_freq, ref_re, ref_im) in reference {
        let freq = freq_vec.real[idx];
        let re = v2_vec.complex[idx].re;
        let im = v2_vec.complex[idx].im;

        // Frequency should match within 0.01%.
        let freq_rel_err = (freq - ref_freq).abs() / ref_freq;
        assert!(
            freq_rel_err < 1e-4,
            "freq[{idx}]: got {freq:.6e}, expected {ref_freq:.6e}, rel_err={freq_rel_err:.2e}"
        );

        // Complex voltage should match within 1% relative to magnitude.
        let ref_mag = (ref_re * ref_re + ref_im * ref_im).sqrt();
        let re_err = (re - ref_re).abs();
        let im_err = (im - ref_im).abs();
        let tol = 0.01 * ref_mag + 1e-10; // 1% of magnitude + small absolute

        assert!(
            re_err < tol,
            "v(2)[{idx}] real: got {re:.6e}, expected {ref_re:.6e}, err={re_err:.2e}, tol={tol:.2e}"
        );
        assert!(
            im_err < tol,
            "v(2)[{idx}] imag: got {im:.6e}, expected {ref_im:.6e}, err={im_err:.2e}, tol={tol:.2e}"
        );
    }
}

/// Verify the -3dB point of a simple RC lowpass filter.
/// f_3dB = 1/(2πRC) for a single-pole RC filter.
#[test]
fn test_rc_lowpass_3db_point() {
    let netlist = Netlist::parse(
        "RC lowpass 3dB test
V1 1 0 DC 0 AC 1
R1 1 2 1k
C1 2 0 1u
.ac DEC 100 10 100k
.end
",
    )
    .unwrap();

    let result = simulate_ac(&netlist).unwrap();
    let plot = &result.plots[0];
    let freq_vec = plot.vecs.iter().find(|v| v.name == "frequency").unwrap();
    let v2_vec = plot.vecs.iter().find(|v| v.name == "v(2)").unwrap();

    let r = 1000.0;
    let c = 1e-6;
    let f_3db = 1.0 / (2.0 * std::f64::consts::PI * r * c); // ~159.15 Hz

    // Find the closest frequency to f_3dB.
    let mut closest_idx = 0;
    let mut min_dist = f64::MAX;
    for (i, &f) in freq_vec.real.iter().enumerate() {
        let dist = ((f / f_3db).ln()).abs();
        if dist < min_dist {
            min_dist = dist;
            closest_idx = i;
        }
    }

    let f_closest = freq_vec.real[closest_idx];
    let re = v2_vec.complex[closest_idx].re;
    let im = v2_vec.complex[closest_idx].im;
    let mag = (re * re + im * im).sqrt();

    // Compute analytical magnitude at f_closest: |H(f)| = 1/√(1 + (f/f_3dB)²)
    let ratio = f_closest / f_3db;
    let expected_mag = 1.0 / (1.0 + ratio * ratio).sqrt();

    assert_abs_diff_eq!(mag, expected_mag, epsilon = 0.01);
}
