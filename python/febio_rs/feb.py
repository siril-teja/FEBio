"""Edit .feb input files for parameter studies."""
import xml.etree.ElementTree as ET


def set_values(src, dst, values):
    """Copy src to dst with values changed. values maps a path to a value:
      "Control/time_steps": 20                      element text
      "Material/material[@name='nucleus']/c1": 0.1  ElementTree XPath (first match)
      "Control/solver/symmetric_stiffness": 0
      "Material/material[@id='3']@name": "disc"     attribute (path@attribute)
    Paths are relative to the root element. Raises KeyError for a path with no match."""
    tree = ET.parse(src)
    root = tree.getroot()
    for path, v in values.items():
        attr = None
        if "@" in path and not path.endswith("]") and path.rsplit("@", 1)[1].isidentifier():
            path, attr = path.rsplit("@", 1)
        el = root.find(path)
        if el is None:
            raise KeyError(f"no element matches {path!r} in {src}")
        if attr:
            el.set(attr, str(v))
        else:
            el.text = str(v)
    tree.write(dst, encoding="ISO-8859-1", xml_declaration=True)
    return dst


def get_value(src, path):
    """Text (or path@attribute) of the first element matching path, or None."""
    root = ET.parse(src).getroot()
    attr = None
    if "@" in path and not path.endswith("]") and path.rsplit("@", 1)[1].isidentifier():
        path, attr = path.rsplit("@", 1)
    el = root.find(path)
    if el is None:
        return None
    return el.get(attr) if attr else el.text
