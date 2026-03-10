#[cfg(test)]
#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);

use thevenin_types::Expr;

/// Parse a numeric expression value, returning an error if it's not a literal number.
pub(crate) fn expr_val(expr: &Expr, context: &str) -> Result<f64, mna::MnaError> {
    match expr {
        Expr::Num(v) => Ok(*v),
        _ => Err(MnaError::NonNumericValue {
            element: context.to_string(),
        }),
    }
}

/// Parse a numeric expression value, returning a default if it's not a literal number.
pub(crate) fn expr_val_or(expr: &Expr, default: f64) -> f64 {
    match expr {
        Expr::Num(v) => *v,
        _ => default,
    }
}

pub mod ac;
pub mod bjt;
pub mod bsim3;
pub mod bsim3soi_dd;
pub mod bsim3soi_fd;
pub mod bsim3soi_pd;
pub mod bsim4;
pub mod cpl;
pub(crate) mod device_stamp;
pub mod diode;
pub mod expr;
pub mod hfet;
pub mod jfet;
pub mod libproc;
pub mod ltra;
pub mod mesa;
pub mod mesfet;
pub mod mna;
pub mod mos6;
pub mod mosfet;
pub mod newton;
pub mod noise;
pub(crate) mod physics;
pub mod pz;
pub mod sens;
pub mod simulate;
pub mod sparse;
pub mod subckt;
pub mod tf;
pub mod transient;
pub mod txl;
pub mod vbic;
pub mod waveform;

pub use ac::simulate_ac;
pub use bjt::{BjtInstance, BjtModel, BjtType, stamp_bjt};
pub use bsim3::{Bsim3Companion, Bsim3Instance, Bsim3Model, stamp_bsim3};
pub use bsim3soi_dd::{
    Bsim3SoiDdCompanion, Bsim3SoiDdInstance, Bsim3SoiDdModel, stamp_bsim3soi_dd,
};
pub use bsim3soi_fd::{
    Bsim3SoiFdCompanion, Bsim3SoiFdInstance, Bsim3SoiFdModel, stamp_bsim3soi_fd,
};
pub use bsim3soi_pd::{
    Bsim3SoiPdCompanion, Bsim3SoiPdInstance, Bsim3SoiPdModel, stamp_bsim3soi_pd,
};
pub use bsim4::{Bsim4Companion, Bsim4Instance, Bsim4Model, stamp_bsim4};
pub use diode::DiodeModel;
pub use hfet::{HfetCompanion, HfetInstance, HfetModel, HfetType};
pub use jfet::{JfetInstance, JfetModel, JfetType, stamp_jfet};
pub use mesa::{MesaCompanion, MesaInstance, MesaModel, stamp_mesa};
pub use mesfet::{MesfetCompanion, MesfetInstance, MesfetModel, MesfetType};
pub use mna::{
    CapacitorInstance, CurrentSourceInstance, DiodeInstance, InductorInstance, MnaError, MnaSystem,
    ResistorInstance, VoltageSourceInstance, assemble_mna, stamp_conductance,
};
pub use mos6::{Mos6Instance, Mos6Model, stamp_mos6};
pub use mosfet::{MosfetInstance, MosfetModel, MosfetType, stamp_mosfet};
pub use newton::{NrError, NrOptions, NrResult, newton_raphson_solve};
pub use noise::simulate_noise;
pub use pz::simulate_pz;
pub use sens::simulate_sens;
pub use simulate::{nr_options_from_netlist, simulate_dc, simulate_op, solve_nonlinear_op};
pub use sparse::{ComplexLinearSystem, LinearSystem, SparseMatrix, SparseMatrixError};
pub use subckt::{SubcktError, flatten_netlist};
pub use tf::simulate_tf;
pub use transient::simulate_tran;
pub use vbic::{
    VbicCompanion, VbicInstance, VbicModel, VbicType, stamp_vbic, stamp_vbic_with_voltages,
};
