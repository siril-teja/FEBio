//! Print a summary of a .feb file: febinfo FILE.feb
fn main() {
    let path = std::env::args().nth(1).expect("usage: febinfo FILE.feb");
    let t = std::time::Instant::now();
    let m = febio_core::read_file(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    println!(
        "febio_spec {} module {:?} ({:.2} s)",
        m.version,
        m.module,
        t.elapsed().as_secs_f64()
    );
    println!("nodes {}", m.node_count());
    println!("solid elements {}", m.solid_element_count());
    println!("element blocks {}", m.blocks.len());
    println!(
        "materials {} + discrete {}",
        m.materials.len(),
        m.discrete_materials.len()
    );
    println!(
        "node sets {}, surfaces {}, discrete sets {}",
        m.node_sets.len(),
        m.surfaces.len(),
        m.discrete_sets.len()
    );
    let secs: Vec<_> = m.sections.iter().map(|s| s.name.as_str()).collect();
    println!("other sections {secs:?}");
}
