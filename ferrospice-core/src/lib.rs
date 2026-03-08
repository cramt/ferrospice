#[cfg(test)]
#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);

pub mod ac;
pub mod bjt;
pub mod diode;
pub mod jfet;
pub mod mna;
pub mod mosfet;
pub mod newton;
pub mod noise;
pub mod pz;
pub mod sens;
pub mod simulate;
pub mod sparse;
pub mod subckt;
pub mod tf;
pub mod transient;
pub mod waveform;

pub use ac::simulate_ac;
pub use bjt::{BjtInstance, BjtModel, BjtType, stamp_bjt};
pub use diode::DiodeModel;
pub use jfet::{JfetInstance, JfetModel, JfetType, stamp_jfet};
pub use mna::{
    CapacitorInstance, CurrentSourceInstance, DiodeInstance, InductorInstance, MnaError, MnaSystem,
    ResistorInstance, VoltageSourceInstance, assemble_mna, stamp_conductance,
};
pub use mosfet::{MosfetInstance, MosfetModel, MosfetType, stamp_mosfet};
pub use newton::{NrError, NrOptions, NrResult, newton_raphson_solve};
pub use noise::simulate_noise;
pub use pz::simulate_pz;
pub use sens::simulate_sens;
pub use simulate::{simulate_dc, simulate_op};
pub use sparse::{ComplexLinearSystem, LinearSystem, SparseMatrix, SparseMatrixError};
pub use subckt::{SubcktError, flatten_netlist};
pub use tf::simulate_tf;
pub use transient::simulate_tran;
