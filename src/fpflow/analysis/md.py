#region modules
import numpy as np 
from lxml import etree 
import h5py 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class MdXml:
    def __init__(self, xml_path: str='md.xml'):
        self.xml_path = xml_path
        self.data = self._parse_xml()

    def _parse_xml(self):
        # Placeholder for XML parsing logic
        tree = etree.parse(self.xml_path)
        root = tree.getroot()

        # Helper: strip namespace
        def lname(tag):
            return tag.split('}')[-1]

        # Collect all <atomic_structure> nodes (initial + steps)
        structures = []
        for node in root.iter():
            if lname(node.tag) == "atomic_structure":
                structures.append(node)

        frames = []
        for struct in structures:
            apos = [n for n in struct if lname(n.tag) == "atomic_positions"][0]
            coords = []
            for atom in apos:
                if lname(atom.tag) != "atom":
                    continue
                coords.append([float(x) for x in atom.text.split()])
            frames.append(coords)

        arr = np.array(frames)  # shape (nstep, natoms, 3)
        return arr 
    
    def write_h5(self):
        with h5py.File('md.h5', 'w') as f:
            f.create_dataset('/positions', data=self.data)

#endregion