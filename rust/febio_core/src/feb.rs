//! Reader for FEBio `.feb` input files, `febio_spec` 2.5 and 4.0 (ROADMAP 4.2).
//!
//! * 2.5: mesh in `<Geometry>` (`Nodes`, `Elements type mat name`, `NodeSet`, `Surface`,
//!   `DiscreteSet`); material assigned per element block by id.
//! * 4.0: mesh in `<Mesh>`; materials assigned by name in `<MeshDomains>`.
//!
//! Files are read as UTF-8, or as ISO-8859-1 when they are not valid UTF-8 (FEBio writes
//! `encoding="ISO-8859-1"`).
use crate::model::{ElemType, ElementBlock, Material, Model, Node};
use std::collections::HashMap;
use std::fmt;

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    Xml(String),
    Format(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Xml(e) => write!(f, "XML error: {e}"),
            Error::Format(e) => write!(f, "invalid .feb file: {e}"),
        }
    }
}

impl std::error::Error for Error {}

fn bad<T>(msg: impl Into<String>) -> Result<T, Error> {
    Err(Error::Format(msg.into()))
}

/// Read a `.feb` file.
pub fn read_file(path: impl AsRef<std::path::Path>) -> Result<Model, Error> {
    let bytes = std::fs::read(path).map_err(Error::Io)?;
    let text = match String::from_utf8(bytes) {
        Ok(s) => s,
        // ISO-8859-1: every byte is the code point of the same value
        Err(e) => e.into_bytes().into_iter().map(char::from).collect(),
    };
    read_str(&text)
}

/// Read a `.feb` document from a string.
pub fn read_str(text: &str) -> Result<Model, Error> {
    let opt = roxmltree::ParsingOptions {
        allow_dtd: true,
        ..Default::default()
    };
    let doc = roxmltree::Document::parse_with_options(text, opt)
        .map_err(|e| Error::Xml(e.to_string()))?;
    let root = doc.root_element();
    if root.tag_name().name() != "febio_spec" {
        return bad(format!(
            "root element is <{}>, not <febio_spec>",
            root.tag_name().name()
        ));
    }
    let mut m = Model {
        version: root.attribute("version").unwrap_or("").to_string(),
        ..Default::default()
    };
    let mut ids: HashMap<i64, usize> = HashMap::new();
    let mut domains: Vec<(String, String)> = Vec::new(); // 4.0 MeshDomains: (part, material name)
    for sec in root.children().filter(|n| n.is_element()) {
        match sec.tag_name().name() {
            "Module" => m.module = sec.attribute("type").map(String::from),
            "Geometry" | "Mesh" => read_mesh(sec, &mut m, &mut ids)?,
            "MeshDomains" => {
                for d in sec.children().filter(|n| n.is_element()) {
                    if let (Some(name), Some(mat)) = (d.attribute("name"), d.attribute("mat")) {
                        domains.push((name.to_string(), mat.to_string()));
                    }
                }
            }
            "Material" => {
                for mat in sec
                    .children()
                    .filter(|n| n.is_element() && n.tag_name().name() == "material")
                {
                    m.materials.push(material(mat));
                }
            }
            "Discrete" => {
                for d in sec
                    .children()
                    .filter(|n| n.is_element() && n.tag_name().name() == "discrete_material")
                {
                    m.discrete_materials.push(material(d));
                }
                m.sections.push(tree(sec));
            }
            _ => m.sections.push(tree(sec)),
        }
    }
    for (part, mat) in domains {
        for b in m
            .blocks
            .iter_mut()
            .filter(|b| b.name.as_deref() == Some(part.as_str()))
        {
            b.material = Some(mat.clone());
        }
    }
    Ok(m)
}

fn tree(n: roxmltree::Node) -> Node {
    Node {
        name: n.tag_name().name().to_string(),
        attrs: n
            .attributes()
            .map(|a| (a.name().to_string(), a.value().to_string()))
            .collect(),
        text: n
            .children()
            .filter(|c| c.is_text())
            .filter_map(|c| c.text())
            .collect::<String>()
            .trim()
            .to_string(),
        children: n.children().filter(|c| c.is_element()).map(tree).collect(),
    }
}

fn material(n: roxmltree::Node) -> Material {
    Material {
        id: n.attribute("id").and_then(|s| s.trim().parse().ok()),
        name: n.attribute("name").map(String::from),
        kind: n.attribute("type").unwrap_or("").to_string(),
        tree: tree(n),
    }
}

fn numbers<T: std::str::FromStr>(s: &str, what: &str) -> Result<Vec<T>, Error> {
    s.split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.is_empty())
        .map(|t| {
            t.parse::<T>()
                .map_err(|_| Error::Format(format!("bad number '{t}' in {what}")))
        })
        .collect()
}

fn text_of(n: roxmltree::Node) -> String {
    n.children()
        .filter(|c| c.is_text())
        .filter_map(|c| c.text())
        .collect()
}

fn index(ids: &HashMap<i64, usize>, id: i64, what: &str) -> Result<usize, Error> {
    ids.get(&id)
        .copied()
        .ok_or_else(|| Error::Format(format!("unknown node {id} in {what}")))
}

fn read_mesh(
    sec: roxmltree::Node,
    m: &mut Model,
    ids: &mut HashMap<i64, usize>,
) -> Result<(), Error> {
    // nodes first (a section may list them after other entities only in theory)
    for nodes in sec
        .children()
        .filter(|n| n.is_element() && n.tag_name().name() == "Nodes")
    {
        for nd in nodes
            .children()
            .filter(|n| n.is_element() && n.tag_name().name() == "node")
        {
            let id: i64 = match nd.attribute("id").and_then(|s| s.trim().parse().ok()) {
                Some(v) => v,
                None => return bad("<node> without a numeric id"),
            };
            let x: Vec<f64> = numbers(&text_of(nd), "node coordinates")?;
            if x.len() != 3 {
                return bad(format!("node {id}: expected 3 coordinates"));
            }
            if ids.insert(id, m.nodes.len()).is_some() {
                return bad(format!("duplicate node id {id}"));
            }
            m.nodes.push([x[0], x[1], x[2]]);
            m.node_ids.push(id);
        }
    }
    for n in sec.children().filter(|n| n.is_element()) {
        match n.tag_name().name() {
            "Elements" => {
                let kind = ElemType::parse(n.attribute("type").unwrap_or(""));
                let mut b = ElementBlock {
                    kind: kind.clone(),
                    name: n.attribute("name").map(String::from),
                    material: n.attribute("mat").map(String::from),
                    ids: Vec::new(),
                    conn: Vec::new(),
                };
                for e in n
                    .children()
                    .filter(|c| c.is_element() && c.tag_name().name() == "elem")
                {
                    let id: i64 = e
                        .attribute("id")
                        .and_then(|s| s.trim().parse().ok())
                        .unwrap_or(-1);
                    let nn: Vec<i64> = numbers(&text_of(e), "element connectivity")?;
                    if let Some(k) = kind.nodes() {
                        if nn.len() != k {
                            return bad(format!("element {id}: {} nodes, expected {k}", nn.len()));
                        }
                    }
                    for v in nn {
                        b.conn.push(index(ids, v, "element connectivity")?);
                    }
                    b.ids.push(id);
                }
                m.blocks.push(b);
            }
            "NodeSet" => {
                let name = n.attribute("name").unwrap_or("").to_string();
                let mut set = Vec::new();
                let kids: Vec<_> = n
                    .children()
                    .filter(|c| c.is_element() && c.tag_name().name() == "node")
                    .collect();
                if kids.is_empty() {
                    for v in numbers::<i64>(&text_of(n), "node set")? {
                        set.push(index(ids, v, "node set")?);
                    }
                } else {
                    for k in kids {
                        let v: i64 = k
                            .attribute("id")
                            .and_then(|s| s.trim().parse().ok())
                            .unwrap_or(-1);
                        set.push(index(ids, v, "node set")?);
                    }
                }
                m.node_sets.insert(name, set);
            }
            "Surface" => {
                let name = n.attribute("name").unwrap_or("").to_string();
                let mut facets = Vec::new();
                for f in n.children().filter(|c| c.is_element()) {
                    let kind = ElemType::parse(f.tag_name().name());
                    let nn: Vec<i64> = numbers(&text_of(f), "surface facet")?;
                    let mut v = Vec::with_capacity(nn.len());
                    for x in nn {
                        v.push(index(ids, x, "surface facet")?);
                    }
                    facets.push((kind, v));
                }
                m.surfaces.insert(name, facets);
            }
            "DiscreteSet" => {
                let name = n.attribute("name").unwrap_or("").to_string();
                let mut pairs = Vec::new();
                for d in n
                    .children()
                    .filter(|c| c.is_element() && c.tag_name().name() == "delem")
                {
                    let nn: Vec<i64> = numbers(&text_of(d), "discrete element")?;
                    if nn.len() != 2 {
                        return bad("discrete element without 2 nodes");
                    }
                    pairs.push([
                        index(ids, nn[0], "discrete set")?,
                        index(ids, nn[1], "discrete set")?,
                    ]);
                }
                m.discrete_sets.insert(name, pairs);
            }
            "Nodes" => {}
            _ => {}
        }
    }
    Ok(())
}
