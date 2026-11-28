//! Data structures of a model read from a `.feb` file.
use std::collections::BTreeMap;

/// Element type (FEBio names).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElemType {
    Hex8,
    Hex20,
    Hex27,
    Penta6,
    Penta15,
    Tet4,
    Tet10,
    Tet15,
    Pyra5,
    Quad4,
    Quad8,
    Quad9,
    Tri3,
    Tri6,
    Line2,
    Other(String),
}

impl ElemType {
    pub fn parse(s: &str) -> ElemType {
        match s {
            "hex8" => ElemType::Hex8,
            "hex20" => ElemType::Hex20,
            "hex27" => ElemType::Hex27,
            "penta6" => ElemType::Penta6,
            "penta15" => ElemType::Penta15,
            "tet4" => ElemType::Tet4,
            "tet10" => ElemType::Tet10,
            "tet15" => ElemType::Tet15,
            "pyra5" => ElemType::Pyra5,
            "quad4" => ElemType::Quad4,
            "quad8" => ElemType::Quad8,
            "quad9" => ElemType::Quad9,
            "tri3" => ElemType::Tri3,
            "tri6" => ElemType::Tri6,
            "line2" => ElemType::Line2,
            o => ElemType::Other(o.to_string()),
        }
    }
    /// number of nodes, if known
    pub fn nodes(&self) -> Option<usize> {
        Some(match self {
            ElemType::Hex8 => 8,
            ElemType::Hex20 => 20,
            ElemType::Hex27 => 27,
            ElemType::Penta6 => 6,
            ElemType::Penta15 => 15,
            ElemType::Tet4 => 4,
            ElemType::Tet10 => 10,
            ElemType::Tet15 => 15,
            ElemType::Pyra5 => 5,
            ElemType::Quad4 => 4,
            ElemType::Quad8 => 8,
            ElemType::Quad9 => 9,
            ElemType::Tri3 => 3,
            ElemType::Tri6 => 6,
            ElemType::Line2 => 2,
            ElemType::Other(_) => return None,
        })
    }
    pub fn is_solid(&self) -> bool {
        matches!(
            self,
            ElemType::Hex8
                | ElemType::Hex20
                | ElemType::Hex27
                | ElemType::Penta6
                | ElemType::Penta15
                | ElemType::Tet4
                | ElemType::Tet10
                | ElemType::Tet15
                | ElemType::Pyra5
        )
    }
}

/// A block of elements of one type (an `<Elements>` section).
#[derive(Clone, Debug)]
pub struct ElementBlock {
    pub kind: ElemType,
    /// part / domain name (`name` attribute), if any
    pub name: Option<String>,
    /// material: id (2.5 `mat` attribute) or name (4.0 `MeshDomains`), if assigned
    pub material: Option<String>,
    /// element ids as in the file
    pub ids: Vec<i64>,
    /// connectivity, 0-based node indices, `kind.nodes()` per element
    pub conn: Vec<usize>,
}

impl ElementBlock {
    pub fn len(&self) -> usize {
        self.ids.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
    pub fn nodes_per_element(&self) -> usize {
        if self.ids.is_empty() {
            0
        } else {
            self.conn.len() / self.ids.len()
        }
    }
    /// 0-based node indices of element `e` of this block
    pub fn element(&self, e: usize) -> &[usize] {
        let k = self.nodes_per_element();
        &self.conn[e * k..(e + 1) * k]
    }
}

/// Generic XML element: name, attributes, text and children.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Node {
    pub name: String,
    pub attrs: BTreeMap<String, String>,
    pub text: String,
    pub children: Vec<Node>,
}

impl Node {
    pub fn attr(&self, k: &str) -> Option<&str> {
        self.attrs.get(k).map(|s| s.as_str())
    }
    pub fn child(&self, name: &str) -> Option<&Node> {
        self.children.iter().find(|c| c.name == name)
    }
    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Node> + 'a {
        self.children.iter().filter(move |c| c.name == name)
    }
    /// numeric value of the text, if it parses
    pub fn value(&self) -> Option<f64> {
        self.text.trim().parse().ok()
    }
}

/// A material (`<material>`): id, name, type and its parameter tree.
#[derive(Clone, Debug)]
pub struct Material {
    pub id: Option<i64>,
    pub name: Option<String>,
    pub kind: String,
    /// the full `<material>` element (parameters, nested solids, fibers, axes)
    pub tree: Node,
}

impl Material {
    /// value of a direct scalar parameter, e.g. `E`, `v`, `c1`
    pub fn param(&self, name: &str) -> Option<f64> {
        self.tree.child(name).and_then(|c| c.value())
    }
}

/// A model read from a `.feb` file.
#[derive(Clone, Debug, Default)]
pub struct Model {
    /// `febio_spec` version, e.g. "2.5" or "4.0"
    pub version: String,
    /// `<Module type>` (e.g. "solid", "biphasic")
    pub module: Option<String>,
    /// node coordinates, in file order
    pub nodes: Vec<[f64; 3]>,
    /// node ids as in the file (same order as `nodes`)
    pub node_ids: Vec<i64>,
    pub blocks: Vec<ElementBlock>,
    /// named node sets: 0-based node indices
    pub node_sets: BTreeMap<String, Vec<usize>>,
    /// named surfaces: facets (type and 0-based node indices)
    pub surfaces: BTreeMap<String, Vec<(ElemType, Vec<usize>)>>,
    /// discrete sets (springs): pairs of 0-based node indices
    pub discrete_sets: BTreeMap<String, Vec<[usize; 2]>>,
    /// solid / rigid materials (`Material` section)
    pub materials: Vec<Material>,
    /// discrete (spring) materials (`Discrete` section)
    pub discrete_materials: Vec<Material>,
    /// all other top-level sections (Control, Boundary, Loads, Contact, Rigid, Step, ...)
    pub sections: Vec<Node>,
}

impl Model {
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    pub fn solid_element_count(&self) -> usize {
        self.blocks
            .iter()
            .filter(|b| b.kind.is_solid())
            .map(|b| b.len())
            .sum()
    }
    pub fn section(&self, name: &str) -> Option<&Node> {
        self.sections.iter().find(|s| s.name == name)
    }
}
