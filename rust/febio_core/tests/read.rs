//! .feb reader tests: small 2.5 and 4.0 models, and (when FINLEY_FEB_DIR points to the
//! unpacked Finley et al. 2018 models) every model against the counts in its FEBio log.
use febio_core::{read_file, read_str, ElemType};

const SPEC25: &str = r#"<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="2.5">
  <Module type="solid"/>
  <Control><time_steps>10</time_steps></Control>
  <Material>
    <material id="1" name="disc" type="neo-Hookean"><E>3.4</E><v>0.4</v></material>
  </Material>
  <Geometry>
    <Nodes>
      <node id="10">0,0,0</node><node id="11">1,0,0</node><node id="12">0,1,0</node>
      <node id="13">0,0,1</node><node id="14">1,1,1</node>
    </Nodes>
    <Elements type="tet4" mat="1" name="part1">
      <elem id="1">10,11,12,13</elem>
      <elem id="2">11,12,13,14</elem>
    </Elements>
    <NodeSet name="bottom"><node id="10"/><node id="11"/></NodeSet>
    <Surface name="top"><tri3 id="1">11,12,13</tri3></Surface>
    <DiscreteSet name="spring"><delem>10,14</delem></DiscreteSet>
  </Geometry>
  <Discrete><discrete_material id="2" name="lig" type="linear spring"><E>1</E></discrete_material></Discrete>
</febio_spec>"#;

const SPEC40: &str = r#"<?xml version="1.0"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Material><material id="1" name="bone" type="neo-Hookean"><E>1000</E><v>0.3</v></material></Material>
  <Mesh>
    <Nodes name="all">
      <node id="1">0,0,0</node><node id="2">1,0,0</node><node id="3">1,1,0</node><node id="4">0,1,0</node>
      <node id="5">0,0,1</node><node id="6">1,0,1</node><node id="7">1,1,1</node><node id="8">0,1,1</node>
    </Nodes>
    <Elements type="hex8" name="block"><elem id="1">1,2,3,4,5,6,7,8</elem></Elements>
    <NodeSet name="base">1,2,3,4</NodeSet>
  </Mesh>
  <MeshDomains><SolidDomain name="block" mat="bone"/></MeshDomains>
  <Boundary><bc name="fix" type="zero displacement" node_set="base"/></Boundary>
</febio_spec>"#;

#[test]
fn reads_spec_2_5() {
    let m = read_str(SPEC25).unwrap();
    assert_eq!(m.version, "2.5");
    assert_eq!(m.module.as_deref(), Some("solid"));
    assert_eq!(m.node_count(), 5);
    assert_eq!(m.node_ids, vec![10, 11, 12, 13, 14]);
    assert_eq!(m.solid_element_count(), 2);
    let b = &m.blocks[0];
    assert_eq!(b.kind, ElemType::Tet4);
    assert_eq!(b.material.as_deref(), Some("1"));
    assert_eq!(b.element(1), &[1, 2, 3, 4]);
    assert_eq!(m.node_sets["bottom"], vec![0, 1]);
    assert_eq!(m.surfaces["top"], vec![(ElemType::Tri3, vec![1, 2, 3])]);
    assert_eq!(m.discrete_sets["spring"], vec![[0, 4]]);
    assert_eq!(m.materials[0].param("E"), Some(3.4));
    assert_eq!(m.discrete_materials[0].kind, "linear spring");
    assert_eq!(
        m.section("Control")
            .unwrap()
            .child("time_steps")
            .unwrap()
            .value(),
        Some(10.0)
    );
}

#[test]
fn reads_spec_4_0() {
    let m = read_str(SPEC40).unwrap();
    assert_eq!(m.version, "4.0");
    assert_eq!(m.node_count(), 8);
    assert_eq!(m.blocks[0].kind, ElemType::Hex8);
    assert_eq!(m.blocks[0].material.as_deref(), Some("bone"));
    assert_eq!(m.node_sets["base"], vec![0, 1, 2, 3]);
    assert!(m.section("Boundary").is_some());
}

#[test]
fn rejects_bad_input() {
    assert!(read_str("<notfebio/>").is_err());
    let bad = SPEC40.replace("1,2,3,4,5,6,7,8", "1,2,3,4,5,6,7,99");
    assert!(read_str(&bad).is_err());
}

fn log_count(log: &str, label: &str) -> Option<usize> {
    let i = log.find(label)?;
    let rest = &log[i..];
    let colon = rest.find(':')?;
    rest[colon + 1..].split_whitespace().next()?.parse().ok()
}

#[test]
fn finley_models_match_febio_logs() {
    let Ok(dir) = std::env::var("FINLEY_FEB_DIR") else {
        eprintln!("FINLEY_FEB_DIR not set: skipping");
        return;
    };
    let mut checked = 0;
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        if p.extension().and_then(|s| s.to_str()) != Some("feb") {
            continue;
        }
        let logp = p.with_extension("log");
        let Ok(log) = std::fs::read(&logp) else {
            continue;
        };
        let log: String = log.into_iter().map(char::from).collect();
        let m = read_file(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
        let name = p.file_name().unwrap().to_string_lossy();
        assert_eq!(
            Some(m.node_count()),
            log_count(&log, "Number of nodes"),
            "{name}: nodes"
        );
        assert_eq!(
            Some(m.solid_element_count()),
            log_count(&log, "Number of solid elements"),
            "{name}: elements"
        );
        assert_eq!(
            Some(m.materials.len() + m.discrete_materials.len()),
            log_count(&log, "Number of materials"),
            "{name}: materials"
        );
        checked += 1;
    }
    assert!(checked > 0, "no .feb with a .log found in {dir}");
    eprintln!("checked {checked} Finley models");
}
