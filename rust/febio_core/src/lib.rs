//! FEBio-Rust FE core (ROADMAP Phase 4).
//!
//! 4.2: [`feb::read_file`] reads FEBio `.feb` input files (`febio_spec` 2.5 and 4.0)
//! into a [`model::Model`]: nodes, element blocks, named node sets and surfaces,
//! discrete (spring) sets and materials. Sections the core does not model yet are kept
//! as generic [`model::Node`] trees.
pub mod feb;
pub mod model;

pub use feb::{read_file, read_str, Error};
pub use model::{ElemType, ElementBlock, Material, Model, Node};
