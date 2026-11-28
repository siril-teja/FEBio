//! FEBio-Rust FE core (ROADMAP Phase 4).
//!
//! 4.2: [`feb::read_file`] reads FEBio `.feb` input files (`febio_spec` 2.5 and 4.0)
//! into a [`model::Model`]: nodes, element blocks, named node sets and surfaces,
//! discrete (spring) sets and materials. Sections the core does not model yet are kept
//! as generic [`model::Node`] trees.
//!
//! 4.3: [`element::ElementTraits`]: hex8, tet4 and tet10 shape functions, FEBio's default
//! integration rules, and kinematics (Jacobian, shape-function gradients, deformation
//! gradient).
//!
//! 4.4: [`material`]: neo-Hookean and orthotropic elastic (Cauchy stress and spatial
//! elasticity tensor), FEBio's formulations.
pub mod element;
pub mod feb;
pub mod material;
pub mod model;

pub use feb::{read_file, read_str, Error};
pub use model::{ElemType, ElementBlock, Material, Model, Node};
